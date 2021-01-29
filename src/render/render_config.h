#pragma once
#include <array>
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>
#include <cstring> // memcmp  // macOS change TH
#include <vector>

namespace render {

enum class RenderMode {
    RenderSlicer,
    RenderMIP,
    RenderIso,
    RenderComposite,
    RenderTF2D
};

struct TFunction {
    //std::array<float, 3> intensities { 0.0f, 0.0f, 0.0f };
    //std::array<float, 3> gradMag { 0.0f, 0.0f, 0.0f };

    float t_intensity;
    float t_radius;
    float t_minGradient, t_maxGradient;

    glm::vec4 color;
};

struct RenderConfig {
    RenderMode renderMode { RenderMode::RenderSlicer };
    glm::ivec2 renderResolution;

    float stepSize = 1.0f;

    bool volumeShading { false };
    float isoValue { 95.0f };

    // 1D transfer function.
    std::array<glm::vec4, 256> tfColorMap;
    // Used to convert from a value to an index in the color map.
    // index = (value - start) / range * tfColorMap.size();
    float tfColorMapIndexStart;
    float tfColorMapIndexRange;

    // 2D transfer function.
    float TF2DIntensity;
    float TF2DRadius;
    glm::vec4 TF2DColor;

    std::vector<TFunction>* TFunctions2D = new std::vector<TFunction>();
};

// NOTE(Mathijs): should be replaced by C++20 three-way operator (aka spaceship operator) if we require C++ 20 support from Linux users (GCC10 / Clang10).
inline bool operator==(const RenderConfig& lhs, const RenderConfig& rhs)
{
    return std::memcmp(&lhs, &rhs, sizeof(RenderConfig)) == 0;
}
inline bool operator!=(const RenderConfig& lhs, const RenderConfig& rhs)
{
    return !(lhs == rhs);
}

}