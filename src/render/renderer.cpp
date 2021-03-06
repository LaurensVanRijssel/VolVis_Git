#include "renderer.h"
#include <algorithm>
#include <algorithm> // std::fill
#include <cmath>
#include <functional>
#include <glm/common.hpp>
#include <glm/gtx/component_wise.hpp>
#include <iostream>
#include <tbb/blocked_range2d.h>
#include <tbb/parallel_for.h>
#include <tuple>

namespace render {

// The renderer is passed a pointer to the volume, gradient volume, camera and an initial renderConfig.
// The camera being pointed to may change each frame (when the user interacts). When the renderConfig
// changes the setConfig function is called with the updated render config. This gives the Renderer an
// opportunity to resize the framebuffer.
Renderer::Renderer(
    const volume::Volume* pVolume,
    const volume::GradientVolume* pGradientVolume,
    const render::RayTraceCamera* pCamera,
    const RenderConfig& initialConfig)
    : m_pVolume(pVolume)
    , m_pGradientVolume(pGradientVolume)
    , m_pCamera(pCamera)
    , m_config(initialConfig)
{
    resizeImage(initialConfig.renderResolution);
}

// Set a new render config if the user changed the settings.
void Renderer::setConfig(const RenderConfig& config)
{
    if (config.renderResolution != m_config.renderResolution)
        resizeImage(config.renderResolution);

    m_config = config;
}

// Resize the framebuffer and fill it with black pixels.
void Renderer::resizeImage(const glm::ivec2& resolution)
{
    m_frameBuffer.resize(size_t(resolution.x) * size_t(resolution.y), glm::vec4(0.0f));
}

// Clear the framebuffer by setting all pixels to black.
void Renderer::resetImage()
{
    std::fill(std::begin(m_frameBuffer), std::end(m_frameBuffer), glm::vec4(0.0f));
}

// Return a VIEW into the framebuffer. This view is merely a reference to the m_frameBuffer member variable.
// This does NOT make a copy of the framebuffer.
gsl::span<const glm::vec4> Renderer::frameBuffer() const
{
    return m_frameBuffer;
}

// Main render function. It computes an image according to the current renderMode.
// Multithreading is enabled in Release/RelWithDebInfo modes. In Debug mode multithreading is disabled to make debugging easier.
void Renderer::render()
{
    resetImage();

    //static constexpr float sampleStep =1.0f;
    float sampleStep = m_config.stepSize;
    const glm::vec3 planeNormal = -glm::normalize(m_pCamera->forward());
    const glm::vec3 volumeCenter = glm::vec3(m_pVolume->dims()) / 2.0f;
    const Bounds bounds { glm::vec3(0.0f), glm::vec3(m_pVolume->dims() - glm::ivec3(1)) };

    // 0 = sequential (single-core), 1 = TBB (multi-core)
#ifdef NDEBUG
    // If NOT in debug mode then enable parallelism using the TBB library (Intel Threaded Building Blocks).
#define PARALLELISM 1
#else
    // Disable multithreading in debug mode.
#define PARALLELISM 0
#endif

#if PARALLELISM == 0
    // Regular (single threaded) for loops.
    for (int x = 0; x < m_config.renderResolution.x; x++) {
        for (int y = 0; y < m_config.renderResolution.y; y++) {
#else
    // Parallel for loop (in 2 dimensions) that subdivides the screen into tiles. 
    const tbb::blocked_range2d<int> screenRange { 0, m_config.renderResolution.y, 0, m_config.renderResolution.x };
        tbb::parallel_for(screenRange, [&](tbb::blocked_range2d<int> localRange) {
        // Loop over the pixels in a tile. This function is called on multiple threads at the same time.
        for (int y = std::begin(localRange.rows()); y != std::end(localRange.rows()); y++) {
            for (int x = std::begin(localRange.cols()); x != std::end(localRange.cols()); x++) {
#endif
            // Compute a ray for the current pixel.
            const glm::vec2 pixelPos = glm::vec2(x, y) / glm::vec2(m_config.renderResolution);
            Ray ray = m_pCamera->generateRay(pixelPos * 2.0f - 1.0f);

            // Compute where the ray enters and exists the volume.
            // If the ray misses the volume then we continue to the next pixel.
            if (!instersectRayVolumeBounds(ray, bounds))
                continue;

            // Get a color for the current pixel according to the current render mode.
            glm::vec4 color {};
            switch (m_config.renderMode) {
            case RenderMode::RenderSlicer: {
                color = traceRaySlice(ray, volumeCenter, planeNormal);
                break;
            }
            case RenderMode::RenderMIP: {
                color = traceRayMIP(ray, sampleStep);
                break;
            }
            case RenderMode::RenderComposite: {
                color = traceRayComposite(ray, sampleStep);
                break;
            }
            case RenderMode::RenderIso: {
                color = traceRayISO(ray, sampleStep);
                break;
            }
            case RenderMode::RenderTF2D: {
                color = traceRayTF2D(ray, sampleStep);
                break;
            }
            };
            // Write the resulting color to the screen.
            fillColor(x, y, color);

#if PARALLELISM == 1
        }
    }
});
#else
            }
        }
#endif
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// This function generates a view alongside a plane perpendicular to the camera through the center of the volume
//  using the slicing technique.
glm::vec4 Renderer::traceRaySlice(const Ray& ray, const glm::vec3& volumeCenter, const glm::vec3& planeNormal) const
{
    const float t = glm::dot(volumeCenter - ray.origin, planeNormal) / glm::dot(ray.direction, planeNormal);
    const glm::vec3 samplePos = ray.origin + ray.direction * t;
    const float val = m_pVolume->getVoxelInterpolate(samplePos);
    return glm::vec4(glm::vec3(std::max(val / m_pVolume->maximum(), 0.0f)), 1.f);
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// Function that implements maximum-intensity-projection (MIP) raycasting.
// It returns the color assigned to a ray/pixel given it's origin, direction and the distances
// at which it enters/exits the volume (ray.tmin & ray.tmax respectively).
// The ray must be sampled with a distance defined by the sampleStep
glm::vec4 Renderer::traceRayMIP(const Ray& ray, float sampleStep) const
{
    float maxVal = 0.0f;

    // Incrementing samplePos directly instead of recomputing it each frame gives a measureable speed-up.
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {
        const float val = m_pVolume->getVoxelInterpolate(samplePos);
        maxVal = std::max(val, maxVal);
    }

    // Normalize the result to a range of [0 to mpVolume->maximum()].
    return glm::vec4(glm::vec3(maxVal) / m_pVolume->maximum(), 1.0f);
}

// ======= TODO: IMPLEMENT ========
// This function should find the position where the ray intersects with the volume's isosurface.
// If volume shading is DISABLED then simply return the isoColor.
// If volume shading is ENABLED then return the phong-shaded color at that location using the local gradient (from m_pGradientVolume).
// Use the camera position (m_pCamera->position()) as the light position.
// Use the bisectionAccuracy function (to be implemented) to get a more precise isosurface location between two steps.
glm::vec4 Renderer::traceRayISO(const Ray& ray, float sampleStep) const
{
    static constexpr glm::vec3 isoColor { 0.8f, 0.8f, 0.2f };

    // Incrementing samplePos directly instead of recomputing it each frame gives a measureable speed-up.
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {

        const float val = m_pVolume->getVoxelInterpolate(samplePos);

        if (val >= m_config.isoValue) {
            // when a high enough Iso value has been found, find the right point in between
            float t_eval = bisectionAccuracy(ray, t -= sampleStep, t, m_config.isoValue);
            samplePos = ray.origin + ray.direction * t_eval;
            // Check if volume shading is enabled, if so compute the reflection color, otherwise, return the isocolor
            if (m_config.volumeShading) {
                return glm::vec4(
                    computePhongShading(
                        isoColor,                                           // Color
                        m_pGradientVolume->getGradientVoxel(samplePos),     // gradient
                        -ray.direction,                                     // light source direction (same as camera)
                        -ray.direction),                                    // camera direction
                    1.0f);
            } else {
                return glm::vec4(isoColor, 1.0f);
            }
        }
    }
    // If no value was above the threshold, return transparent black.
    return glm::vec4(0.0f);
}

// ======= TODO: IMPLEMENT ========
// Given that the iso value lies somewhere between t0 and t1, find a t for which the value
// closely matches the iso value (less than 0.01 difference). Add a limit to the number of
// iterations such that it does not get stuck in degerate cases.
float Renderer::bisectionAccuracy(const Ray& ray, float t0, float t1, float isoValue) const
{
    const float tolerance = 0.01f;
    const int max_iterations = 5; // range decreases to 3,1% of start range (t1 - t0)

    float t_mid(0.0);

    for (int i(0); i < max_iterations; ++i) {

        // estimate t_mid, get closer every step.
        t_mid = (t0 + t1) / 2.0f;
        float t_diff = t_mid - t0;
        const glm::vec3 p_mid = ray.origin + ray.direction * t_mid;
        const float v_mid = m_pVolume->getVoxelInterpolate(p_mid);

        // if v_mid is higher than isoValue, the IsoValue lays in between t0 and t_mid, otherwise in between t_mid and t1
        // If the step size is below 1, there are only 2 voxels, in this case, return the right one.
        if (v_mid >= isoValue) {
            if (t_diff < 1.0f) {
                return t_mid;
            }
            t1 = t_mid;
        } else {
            if (t_diff < 1.0f) {
                return t1;
            }
            t0 = t_mid;
        }

        // if v_mid is very close to the isoValue, return t_mid
        if (abs(v_mid - isoValue) < tolerance) {
            return t_mid;
        }
    }
    // if it has been tried more than the max_iter, return t_mid, is will be close to the actual value.
    return t_mid;
}

// ======= TODO: IMPLEMENT ========
// In this function, implement 1D transfer function raycasting.
// Use getTFValue to compute the color for a given volume value according to the 1D transfer function.
glm::vec4 Renderer::traceRayComposite(const Ray& ray, float sampleStep) const
{
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    glm::vec3 ct(0.0f);
    float at = 0.0f;

    // for each step, interpolate the intensity and compute the color and opacity accordingly.
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {

        const float intensityVal = m_pVolume->getVoxelInterpolate(samplePos);
        const glm::vec4 volumeValue = getTFValue(intensityVal);

        glm::vec3 ci;
        // if volume shading is on, update the color based on the gradient.
        if (m_config.volumeShading) {
            ci = computePhongShading(
                glm::vec3(volumeValue.r, volumeValue.g, volumeValue.b), // Color
                m_pGradientVolume->getGradientVoxel(samplePos),         // gradient
                -ray.direction,                                         // Light source is at the camera position.
                -ray.direction);                                        // camera direction
        } else {
            ci = glm::vec3( volumeValue.r, volumeValue.g, volumeValue.b );
        }
 
        float ai = volumeValue.a;
        
        ci *= ai;

        ct += (1 - at) * ci;
        at = at + (1 - at) * ai;

        // if the total opacity is close to 1, the remaining voxel do not matter.
        if (at >= 0.99f) {
            break;
        }
    }

    return glm::vec4(ct, at);
}

// ======= TODO: IMPLEMENT ========
// In this function, implement 2D transfer function raycasting.
// Use the getTF2DOpacity function that you implemented to compute the opacity according to the 2D transfer function.
glm::vec4 Renderer::traceRayTF2D(const Ray& ray, float sampleStep) const
{
    glm::vec3 samplePos = ray.origin + ray.tmin * ray.direction;
    const glm::vec3 increment = sampleStep * ray.direction;

    glm::vec3 composite_color(0.0f);
    float composite_opacity = 0.0f;
    for (float t = ray.tmin; t <= ray.tmax; t += sampleStep, samplePos += increment) {

        //Determine the intensity and gradient of this voxel
        const float intensityVal = m_pVolume->getVoxelInterpolate(samplePos);
        const volume::GradientVoxel& gradient = m_pGradientVolume->getGradientVoxel(samplePos);
        const float gradientMag = gradient.magnitude;

        //Determine what transfer function this voxel belongs to.
        const int tfNum = determine2DTF(intensityVal, gradientMag);
        //If it belongs to a transfer funtion, determine the color and opacity.
        //Otherwise we can just go to the next voxel on this ray.
        if (tfNum >= 0) {

            //Get the opacity and color of this voxel
            const float opacity = getTF2DOpacity(intensityVal, gradientMag, tfNum);
            const glm::vec4 color = m_config.TFunctions2D->at(tfNum).color;

            //If the opacity of this voxel is 0, just go to the next voxel in along this ray.
            if (opacity != 0.0f) {
                glm::vec3 ci;
                if (m_config.volumeShading) {
                    ci = computePhongShading(
                        glm::vec3(color.r, color.g, color.b), // Color
                        gradient, // gradient
                        //glm::normalize(m_pCamera->position() - samplePos),  // Light source
                        glm::normalize(-ray.direction), // is the same as above, but quicker
                        glm::normalize(-ray.direction)); // camera direction
                } else {
                    ci = glm::vec3(color.r, color.g, color.b);
                }

                ci = ci * opacity;

                //Composite the color value and opacity of previous voxels with this voxel
                composite_color = composite_color + (1 - opacity) * ci;
                composite_opacity = composite_opacity + (1 - composite_opacity) * opacity;
            }
        }

    }

    //If the total opacity for this pixel is 0 we want to display a empty (black) pixel.
    if (composite_opacity == 0.0f) {
        composite_color = glm::vec3(0.0f, 0.0f, 0.0f);
    }

    glm::vec4 returnVal(composite_color, composite_opacity);

    return returnVal;
}

// ======= TODO: IMPLEMENT ========
// Compute Phong Shading given the voxel color (material color), the gradient, the light vector and view vector.
// You can find out more about the Phong shading model at:
// https://en.wikipedia.org/wiki/Phong_reflection_model
//
// Use the given color for the ambient/specular/diffuse (you are allowed to scale these constants by a scalar value).
// You are free to choose any specular power that you'd like.
glm::vec3 Renderer::computePhongShading(const glm::vec3& color, const volume::GradientVoxel& gradient, const glm::vec3& L, const glm::vec3& V)
{

    // define the constants
    const float ka(0.1f), kd(0.7f), ks(0.2f);
    const int alpha(10);

    // if the normal vector is clos to 0, the normal vector can be disregarded, in this case only ambient reflection is used.
    if (gradient.magnitude < 0.001) {
        return ka * color;
    }
    const glm::vec3 N = glm::dot(gradient.dir, L) < 0 ? glm::normalize(-gradient.dir) : glm::normalize(gradient.dir);
    // Get the normalized normal vector (gradient) from voxel
    //const glm::vec3 N = glm::normalize(gradient.dir);

    // compute reflection ray, and determine collinearity between reflection ray and viewer ray for specular reflection
    const glm::vec3 R = 2.0f * glm::dot(N, L) * N - L;
    const float VR = pow(glm::dot(V, R), alpha);
    
    glm::vec3 Itotal = 
        ka * color +                    // Ambient reflection
        kd * glm::dot(L, N) * color +   // Diffuse reflection
        ks * VR * color;                // specular reflection

    //return the sum of the ambient, diffuse and specular reflection.
    return Itotal;
}

// ======= DO NOT MODIFY THIS FUNCTION ========
// Looks up the color+opacity corresponding to the given volume value from the 1D tranfer function LUT (m_config.tfColorMap).
// The value will initially range from (m_config.tfColorMapIndexStart) to (m_config.tfColorMapIndexStart + m_config.tfColorMapIndexRange) .
glm::vec4 Renderer::getTFValue(float val) const
{
    // Map value from [m_config.tfColorMapIndexStart, m_config.tfColorMapIndexStart + m_config.tfColorMapIndexRange) to [0, 1) .
    const float range01 = (val - m_config.tfColorMapIndexStart) / m_config.tfColorMapIndexRange;
    const size_t i = std::min(static_cast<size_t>(range01 * static_cast<float>(m_config.tfColorMap.size())), m_config.tfColorMap.size()-1);
    return m_config.tfColorMap[i];
}

// ======= SELF ADDED =============
//Determine which 2D transfer functions that the point defined by the intensity and gradientMagnitude lies in.
//Return the index of that transfer function
//If it does not lie in any transfer function, return -1
int Renderer::determine2DTF(float intensity, float gradientMagnitude) const
{
    glm::vec2 target_point(intensity, gradientMagnitude);

    for (int i = 0; i < m_config.TFunctions2D->size(); ++i) {
        TFunction tf = m_config.TFunctions2D->at(i);

        //If this point clearly does not lie inside the transfer function, skip this transfer function.
        if (intensity < (tf.t_intensity - tf.t_radius)) {
            continue;
        }
        if (intensity > (tf.t_intensity + tf.t_radius)) {
            continue;
        }

        glm::vec2 p1(tf.t_intensity, tf.t_minGradient * m_pGradientVolume->maxMagnitude());
        glm::vec2 p2(tf.t_intensity - tf.t_radius, tf.t_maxGradient * m_pGradientVolume->maxMagnitude());
        glm::vec2 p3(tf.t_intensity + tf.t_radius, tf.t_maxGradient * m_pGradientVolume->maxMagnitude());

        //Define the matrix's that are created by the different points
        glm::mat2x2 mt2(target_point, p2);
        glm::mat2x2 mt3(target_point, p3);

        glm::mat2x2 m12(p1, p2);
        glm::mat2x2 m13(p1, p3);
        glm::mat2x2 m23(p2, p3);

        const float dt2 = glm::determinant(mt2);
        const float dt3 = glm::determinant(mt3);

        const float d12 = glm::determinant(m12);
        const float d13 = glm::determinant(m13);
        const float d23 = glm::determinant(m23);

        //To determine if the point lies inside the triangle we check if the convex hull contains 3 or 4 points.
        //If the convex hull contains 4 points the point lies outside the triangle.
        //If the convex hull contains 3 points the point lies inside the triangle.
        //Here is a definition of how this works:
        //https://mathworld.wolfram.com/TriangleInterior.html
        const float a = ((dt3 - d13) / d23);
        const float b = -1 * ((dt2 - d12) / d23);

        if ((a > 0) && (b > 0)) {
            if ((a + b) > 0) {
                return i;
            }
        }
    }
     return -1;
}

// ======= TODO: IMPLEMENT ========
// This function should return an opacity value for the given intensity and gradient according to the 2D transfer function.
// Calculate whether the values are within the radius/intensity triangle defined in the 2D transfer function widget.
// If so: return a tent weighting as described in the assignment
// Otherwise: return 0.0f
//
// The 2D transfer function settings can be accessed through m_config.TF2DIntensity and m_config.TF2DRadius.
float  Renderer::getTF2DOpacity(float intensity, float gradientMagnitude, int tfNum) const
{
    TFunction tf = m_config.TFunctions2D->at(tfNum);

    float intDifferance = abs(tf.t_intensity - intensity);
    float normalizedDifferance = intDifferance / tf.t_radius;

    float intValue = 1 - normalizedDifferance;

    return intValue;
}

// This function computes if a ray intersects with the axis-aligned bounding box around the volume.
// If the ray intersects then tmin/tmax are set to the distance at which the ray hits/exists the
// volume and true is returned. If the ray misses the volume the the function returns false.
//
// If you are interested you can learn about it at.
// https://www.scratchapixel.com/lessons/3d-basic-rendering/minimal-ray-tracer-rendering-simple-shapes/ray-box-intersection
bool Renderer::instersectRayVolumeBounds(Ray& ray, const Bounds& bounds) const
{
    const glm::vec3 invDir = 1.0f / ray.direction;
    const glm::bvec3 sign = glm::lessThan(invDir, glm::vec3(0.0f));

    float tmin = (bounds.lowerUpper[sign[0]].x - ray.origin.x) * invDir.x;
    float tmax = (bounds.lowerUpper[!sign[0]].x - ray.origin.x) * invDir.x;
    const float tymin = (bounds.lowerUpper[sign[1]].y - ray.origin.y) * invDir.y;
    const float tymax = (bounds.lowerUpper[!sign[1]].y - ray.origin.y) * invDir.y;

    if ((tmin > tymax) || (tymin > tmax))
        return false;
    tmin = std::max(tmin, tymin);
    tmax = std::min(tmax, tymax);

    const float tzmin = (bounds.lowerUpper[sign[2]].z - ray.origin.z) * invDir.z;
    const float tzmax = (bounds.lowerUpper[!sign[2]].z - ray.origin.z) * invDir.z;

    if ((tmin > tzmax) || (tzmin > tmax))
        return false;

    ray.tmin = std::max(tmin, tzmin);
    ray.tmax = std::min(tmax, tzmax);
    return true;
}

// This function inserts a color into the framebuffer at position x,y
void Renderer::fillColor(int x, int y, const glm::vec4& color)
{
    const size_t index = static_cast<size_t>(m_config.renderResolution.x * y + x);
    m_frameBuffer[index] = color;
}
}