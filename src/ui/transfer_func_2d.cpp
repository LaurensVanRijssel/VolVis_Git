#include "transfer_func_2d.h"
#include <algorithm>
#include <array>
#include <filesystem>
#include <glm/common.hpp>
#include <glm/gtc/type_ptr.hpp>
#include <glm/gtx/string_cast.hpp>
#include <imgui.h>
#include <iostream>
#include <vector>

static ImVec2 glmToIm(const glm::vec2& v);
static glm::vec2 ImToGlm(const ImVec2& v);
static std::vector<glm::vec4> createHistogramImage(
    const volume::Volume& volume, const volume::GradientVolume& gradient, const glm::ivec2& res);

namespace ui {

// Radius of the three points in the histogram image.
static constexpr float pointRadius = 8.0f;
static constexpr glm::ivec2 widgetSize { 475, 300 };

TransferFunction2DWidget::TransferFunction2DWidget(const volume::Volume& volume, const volume::GradientVolume& gradient)
    : m_intensity(68.0f)
    , m_maxIntensity(volume.maximum())
    , m_radius(38.0f)
    , m_color(0.0f, 0.8f, 0.6f, 0.3f)
    , m_interactingPoint(-1)
    , m_histogramImg(0)
{
    const glm::ivec2 res = glm::ivec2(volume.maximum(), gradient.maxMagnitude() + 1);
    const auto imgData = createHistogramImage(volume, gradient, res);

    //Initialize the first Transfer Function
    TFunction initial_TF;
    /*initial_TF.intensities[0] = 68.0f;
    initial_TF.intensities[1] = 68.0f - 38.0f;
    initial_TF.intensities[2] = 68.0f + 38.0f;

    initial_TF.gradMag[0] = 0.0f;
    initial_TF.gradMag[1] = 1.0f;
    initial_TF.gradMag[2] = 1.0f;*/
    initial_TF.t_intensity = 68.0f;
    initial_TF.t_radius = 38.0f;
    initial_TF.t_minGradient = 0.0f;
    initial_TF.t_maxGradient = 1.0f;
    initial_TF.color = glm::vec4(0.0f, 0.8f, 0.6f, 0.3f);
    m_tfunctions.push_back(initial_TF);
    currentlySelectedTF = 0;
    
    glGenTextures(1, &m_histogramImg);
    glBindTexture(GL_TEXTURE_2D, m_histogramImg);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glBindTexture(GL_TEXTURE_2D, m_histogramImg);
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, res.x, res.y, 0, GL_RGBA, GL_FLOAT, imgData.data());

}

// Draw the widget and handle interactions
void TransferFunction2DWidget::draw()
{
    const ImGuiIO& io = ImGui::GetIO();

    ImGui::Text("2D Transfer Function");
    ImGui::Text("Click and drag points to alter the m_radius or m_intensity");

    // Histogram image is positioned to the right of the content region.
    const glm::vec2 canvasSize { widgetSize.x, widgetSize.y - 20 };
    glm::vec2 canvasPos = ImToGlm(ImGui::GetCursorScreenPos()); // this is the imgui draw cursor, not mouse cursor
    const float xOffset = (ImToGlm(ImGui::GetContentRegionAvail()).x - canvasSize.x);
    canvasPos.x += xOffset; // center widget

    // Draw side text (imgui cannot center-align text so we have to do it ourselves).
    ImVec2 cursorPos = ImGui::GetCursorPos();
    ImGui::SetCursorPos(ImVec2(cursorPos.x + 25, cursorPos.y + canvasSize.y / 2 - 20 - ImGui::GetFontSize()));
    //ImGui::Text("");
    ImGui::Separator();
    ImGui::SetCursorPos(ImVec2(cursorPos.x + 15, cursorPos.y + canvasSize.y / 2 - 20));
    ImGui::Text("Gradient");
    ImGui::SetCursorPos(ImVec2(cursorPos.x + 5, cursorPos.y + canvasSize.y / 2 - 20 + ImGui::GetFontSize()));
    ImGui::Text("Magnitude");
    ImGui::SetCursorPos(cursorPos);

    // Draw box and histogram image.
    ImDrawList* drawList = ImGui::GetWindowDrawList();
    drawList->PushClipRect(glmToIm(canvasPos), glmToIm(canvasPos + canvasSize));
    drawList->AddRect(glmToIm(canvasPos), glmToIm(canvasPos + canvasSize), ImColor(180, 180, 180, 255));

    cursorPos = ImVec2(ImGui::GetCursorPosX() + xOffset, ImGui::GetCursorPosY());

    // Draw histogram image that we uploaded to the GPU using OpenGL.
    ImGui::SetCursorPos(cursorPos);

    // NOTE(Mathijs): reinterpret casting the pointer is undefined behavior according to the standard.
    // Use memcpy for now and move to std::bit_cast (C++20) when more compilers support it.
    //
    // From the standard:
    // When it is needed to interpret the bytes of an object as a value of a different type, std::memcpy or std::bit_cast (since C++20)can be used:
    // https://en.cppreference.com/w/cpp/language/reinterpret_cast
    ImTextureID imguiTexture;
    std::memcpy(&imguiTexture, &m_histogramImg, sizeof(m_histogramImg));
    ImGui::Image(imguiTexture, glmToIm(canvasSize - glm::vec2(1)));

    // Detect and handle mouse interaction.
    if (!io.MouseDown[0] && !io.MouseDown[1]) {
        m_interactingPoint = -1;
    }

    // Place an invisible button on top of the histogram. IsItemHovering returns whether the cursor is
    // hovering over the last added item which in this case is the invisble button. This way we can
    // easily detect whether the cursor is inside the histogram image.
    ImGui::SetCursorPos(cursorPos);
    ImGui::InvisibleButton("tfn_canvas", glmToIm(canvasSize));

    // Mouse position within the histogram image.
    const ImVec2 bbMin = ImGui::GetItemRectMin();
    const ImVec2 bbMax = ImGui::GetItemRectMax();
    const glm::vec2 clippedMousePos {
        std::min(std::max(io.MousePos.x, bbMin.x), bbMax.x),
        std::min(std::max(io.MousePos.y, bbMin.y), bbMax.y)
    };

    const float normalizedIntensity = m_intensity / m_maxIntensity;
    const float normalizedRadius = m_radius / m_maxIntensity;

    const glm::vec2 viewScale { canvasSize.x, -canvasSize.y };
    const glm::vec2 viewOffset { canvasPos.x, canvasPos.y + canvasSize.y };

    std::vector<glm::vec2> points;

    for (int i = 0; i < m_tfunctions.size(); ++i) {
        TFunction transferfunction = m_tfunctions[i];

        /*const glm::vec2 p1(transferfunction.intensities[0] / m_maxIntensity, transferfunction.gradMag[0]);
        const glm::vec2 p2(transferfunction.intensities[1] / m_maxIntensity, transferfunction.gradMag[1]);
        const glm::vec2 p3(transferfunction.intensities[2] / m_maxIntensity, transferfunction.gradMag[2]);*/

        float t_normalizedIntensity = transferfunction.t_intensity / m_maxIntensity;
        float t_normalizedRadius = transferfunction.t_radius / m_maxIntensity;
        const glm::vec2 p1(t_normalizedIntensity, transferfunction.t_minGradient);
        const glm::vec2 p2(t_normalizedIntensity - t_normalizedRadius, transferfunction.t_maxGradient);
        const glm::vec2 p3(t_normalizedIntensity + t_normalizedRadius, transferfunction.t_maxGradient);

        points.push_back(p1);
        points.push_back(p2);
        points.push_back(p3);

        const ImVec2 point1Pos = glmToIm(p1 * viewScale + viewOffset);
        const ImVec2 point2Pos = glmToIm(p2 * viewScale + viewOffset);
        const ImVec2 point3Pos = glmToIm(p3 * viewScale + viewOffset);

        //glm::

        //ImU32 test = transferfunction.color;
        //float r = transferfunction.color[0];
        ImVec4 imv4(transferfunction.color.r, transferfunction.color.g, transferfunction.color.b, 1.0f);
        ImU32 test = ImGui::ColorConvertFloat4ToU32(imv4);
        ImU32 colorTest = 0xFFFFFFFF;

        drawList->AddLine(point1Pos, point2Pos, test);
        drawList->AddLine(point1Pos, point3Pos, test);
        drawList->AddLine(point2Pos, point3Pos, test);

        drawList->AddCircleFilled(point1Pos, pointRadius, test);
        drawList->AddCircleFilled(point2Pos, pointRadius, test);
        drawList->AddCircleFilled(point3Pos, pointRadius, test);

    }

    /*const std::array points = {
        glm::vec2(normalizedIntensity, 0.f),
        glm::vec2(normalizedIntensity - normalizedRadius, 1.f),
        glm::vec2(normalizedIntensity + normalizedRadius, 1.f)
    };*/

    if (ImGui::IsItemHovered() && (io.MouseDown[0] || io.MouseDown[1])) {
        const glm::vec2 mousePos = glm::clamp((clippedMousePos - viewOffset) / viewScale, 0.0f, 1.0f);

        if (io.MouseDown[0]) {
            // Left mouse button is down.
            if (m_interactingPoint == -1) {
                // No point is currently selected; check if the user clicked any of the points.
                for (size_t i = 0; i < points.size(); i++) {
                    const glm::vec2 pointPosition = points[i] * viewScale + viewOffset;
                    const glm::vec2 diff = pointPosition - clippedMousePos;
                    float distanceSquared = glm::dot(diff, diff);
                    if (distanceSquared < pointRadius * pointRadius) {
                        m_interactingPoint = int(i);
                        break;
                    }
                }
            }

            if (m_interactingPoint != -1) {

                float new_intensity = mousePos.x * m_maxIntensity;
                float new_gradient = mousePos.y;

                //Determine which transferfunction this point belongs too.
                //Since there are 3 points per transfer function we divide this number by 3.
                int tf_num = glm::floor(m_interactingPoint / 3);
                int point_num = m_interactingPoint % 3;

                /*m_tfunctions[tf_num].intensities[point_num] = new_intensity;
                m_tfunctions[tf_num].gradMag[point_num] = new_gradient;*/

                // A point was already selected.
                switch (point_num) {
                case 0: {
                    // m_intensity point
                    //m_intensity = mousePos.x * m_maxIntensity;
                    m_tfunctions[tf_num].t_intensity = new_intensity;
                    m_tfunctions[tf_num].t_minGradient = new_gradient;
                } break;
                case 1: {
                    // left m_radius point
                    //m_radius = std::max(m_intensity - mousePos.x * m_maxIntensity, 1.0f);
                    m_tfunctions[tf_num].t_radius = std::max(m_tfunctions[tf_num].t_intensity - new_intensity, 1.0f);
                    m_tfunctions[tf_num].t_maxGradient = std::max(new_gradient, m_tfunctions[tf_num].t_minGradient + 0.05f);
                } break;
                case 2: {
                    // right m_radius point
                    //m_radius = std::max(mousePos.x * m_maxIntensity - m_intensity, 1.0f);
                    m_tfunctions[tf_num].t_radius = std::max(new_intensity - m_tfunctions[tf_num].t_intensity, 1.0f);
                    m_tfunctions[tf_num].t_maxGradient = std::max(new_gradient, m_tfunctions[tf_num].t_minGradient + 0.05f);
                } break;
                };

                currentlySelectedTF = tf_num;
            }
        }
    }

    /*for (int i = 0; i < m_tfunctions.size(); ++i) {
        TFunction transferfunction = m_tfunctions[i];
        const ImVec2 point1Pos = glmToIm(glm::vec2(transferfunction.intensity1 / m_maxIntensity, transferfunction.gradMag1) * viewScale + viewOffset);
        const ImVec2 point2Pos = glmToIm(glm::vec2(transferfunction.intensity2 / m_maxIntensity, transferfunction.gradMag2) * viewScale + viewOffset);
        const ImVec2 point3Pos = glmToIm(glm::vec2(transferfunction.intensity3 / m_maxIntensity, transferfunction.gradMag3) * viewScale + viewOffset);

        drawList->AddLine(point1Pos, point2Pos, 0xFFFFFFFF);
        drawList->AddLine(point1Pos, point3Pos, 0xFFFFFFFF);

        drawList->AddCircleFilled(point1Pos, pointRadius, 0xFFFFFFFF);
        drawList->AddCircleFilled(point2Pos, pointRadius, 0xFFFFFFFF);
        drawList->AddCircleFilled(point3Pos, pointRadius, 0xFFFFFFFF);
    }*/

    // Draw m_intensity and m_radius points.
    /*const ImVec2 leftRadiusPointPos = glmToIm(glm::vec2(normalizedIntensity - normalizedRadius, 1.f) * viewScale + viewOffset);
    const ImVec2 rightRadiusPointPos = glmToIm(glm::vec2(normalizedIntensity + normalizedRadius, 1.f) * viewScale + viewOffset);
    const ImVec2 intensityPointPos = glmToIm(glm::vec2(normalizedIntensity, 0.f) * viewScale + viewOffset);

    drawList->AddLine(leftRadiusPointPos, intensityPointPos, 0xFFFFFFFF);
    drawList->AddLine(intensityPointPos, rightRadiusPointPos, 0xFFFFFFFF);

    drawList->AddCircleFilled(leftRadiusPointPos, pointRadius, 0xFFFFFFFF);
    drawList->AddCircleFilled(intensityPointPos, pointRadius, 0xFFFFFFFF);
    drawList->AddCircleFilled(rightRadiusPointPos, pointRadius, 0xFFFFFFFF);*/

    //const ImVec2 testPoint = glmToIm(glm::vec2(0.5f, 0.5f) * viewScale + viewOffset);
    //drawList->AddCircleFilled(testPoint, pointRadius, 0xFFFFFFFF);

    drawList->PopClipRect();

    // Bottom text.
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + xOffset + canvasSize.x / 2 - 40);
    ImGui::Text("Voxel Value");

    if (ImGui::Button("Add Function")) {
        const ImVec2 testPoint = glmToIm(glm::vec2(0.8f, 0.8f) * viewScale + viewOffset);
        drawList->AddCircleFilled(testPoint, pointRadius, 0xFFFFFFFF);

        TFunction newTF;
        newTF.t_intensity = 68.0f;
        newTF.t_radius = 38.0f;
        newTF.t_minGradient = 0.0f;
        newTF.t_maxGradient = 1.0f;
        newTF.color = glm::vec4(0.0f, 0.8f, 0.6f, 0.3f);
        m_tfunctions.push_back(newTF);

        currentlySelectedTF = m_tfunctions.size() - 1;
    }
    if (ImGui::Button("Remove Function")) {
        m_tfunctions.erase(m_tfunctions.begin() + currentlySelectedTF);

        currentlySelectedTF = 0;
    }

    // Draw text boxes and color picker.
    ImGui::NewLine();
    int inputFill = 0;
    ImGui::PushItemWidth(0.1f); // These 3 coming lines show nothing but vertically align the text of the intensity and radius boxes
    ImGui::InputScalar("", ImGuiDataType_U32, &inputFill);
    ImGui::SameLine();
    ImGui::Text("Intensity: ");
    ImGui::SameLine();
    ImGui::PushItemWidth(50.f);
    ImGui::InputScalar("", ImGuiDataType_Float, &m_tfunctions[currentlySelectedTF].t_intensity, NULL, NULL, "%.2f", ImGuiInputTextFlags_ReadOnly);
    ImGui::SameLine();
    ImGui::Text("Radius: ");
    ImGui::SameLine();
    ImGui::PushItemWidth(50.f);
    ImGui::InputScalar("", ImGuiDataType_Float, &m_tfunctions[currentlySelectedTF].t_radius, NULL, NULL, "%.2f", ImGuiInputTextFlags_ReadOnly);
    ImGui::NewLine();
    ImGui::Text("Min Gradient: ");
    ImGui::SameLine();
    ImGui::PushItemWidth(50.f);
    ImGui::InputScalar("", ImGuiDataType_Float, &m_tfunctions[currentlySelectedTF].t_minGradient, NULL, NULL, "%.2f", ImGuiInputTextFlags_ReadOnly);
    ImGui::SameLine();
    ImGui::Text("Max Gradient: ");
    ImGui::SameLine();
    ImGui::PushItemWidth(50.f);
    ImGui::InputScalar("", ImGuiDataType_Float, &m_tfunctions[currentlySelectedTF].t_maxGradient, NULL, NULL, "%.2f", ImGuiInputTextFlags_ReadOnly);

    ImGui::NewLine();
    ImGui::SetCursorPosX(ImGui::GetCursorPosX() + xOffset / 2);
    ImGui::PushItemWidth(ImGui::GetContentRegionAvailWidth() * 0.4f);
    ImGui::ColorPicker4("Color", glm::value_ptr(m_tfunctions[currentlySelectedTF].color));
}

void TransferFunction2DWidget::updateRenderConfig(render::RenderConfig& renderConfig)
{
    
    /*renderConfig.TF2DIntensity = m_intensity;
    renderConfig.TF2DRadius = m_radius;*/
    //renderConfig.TFunctions2D = m_tfunctions;
    
    renderConfig.TFunctions2D->clear();

    for (int i = 0; i < m_tfunctions.size(); ++i) {
        TFunction tf = m_tfunctions[i];
        render::TFunction newTF;
        /*newTF.intensities[0] = tf.intensities[0];
        newTF.intensities[1] = tf.intensities[1];
        newTF.intensities[2] = tf.intensities[2];

        newTF.gradMag[0] = tf.gradMag[0];
        newTF.gradMag[1] = tf.gradMag[1];
        newTF.gradMag[2] = tf.gradMag[2];*/
        newTF.t_intensity = tf.t_intensity;
        newTF.t_radius = tf.t_radius;
        newTF.t_minGradient = tf.t_minGradient;
        newTF.t_maxGradient = tf.t_maxGradient;

        newTF.color = tf.color;

        renderConfig.TFunctions2D->push_back(newTF);
    }

    renderConfig.TF2DColor = m_color;
}
}

static ImVec2 glmToIm(const glm::vec2& v)
{
    return ImVec2(v.x, v.y);
}

static glm::vec2 ImToGlm(const ImVec2& v)
{
    return glm::vec2(v.x, v.y);
}

// Compute a histogram texture from the volume and gradient data
static std::vector<glm::vec4> createHistogramImage(
    const volume::Volume& volume, const volume::GradientVolume& gradient, const glm::ivec2& res)
{
    const size_t numPixels = static_cast<size_t>(res.x * res.y);
    std::vector<int> bins(numPixels, 0);
    for (int z = 0; z < volume.dims().z; z++) {
        for (int y = 0; y < volume.dims().y; y++) {
            for (int x = 0; x < volume.dims().x; x++) {
                const size_t imgX = static_cast<size_t>(volume.getVoxel(x, y, z));
                const size_t imgY = static_cast<size_t>(res.y - 1) - static_cast<size_t>(gradient.getGradientVoxel(x, y, z).magnitude);
                const size_t index = imgX + imgY * static_cast<size_t>(res.x);
                bins[index]++;
            }
        }
    }

    const int maxCount = *std::max_element(std::begin(bins), std::end(bins));
    const float factor = 1.0f / std::log(float(maxCount));

    std::vector<glm::vec4> imageData(numPixels, glm::vec4(0.0f));
    std::transform(std::begin(bins), std::end(bins), std::begin(imageData),
        [&](int count) {
            return glm::vec4(1.0f, 1.0f, 1.0f, std::log(float(count)) * factor);
        });
    return imageData;
}
