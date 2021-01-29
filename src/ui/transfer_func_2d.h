#pragma once
#include "render/render_config.h"
#include "volume/gradient_volume.h"
#include "volume/volume.h"
#include <GL/glew.h> // Include before glfw3
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

namespace ui {

struct TFunction {
    //std::array<float, 3> intensities {0.0f, 0.0f, 0.0f};
    //std::array<float, 3> gradMag { 0.0f, 0.0f, 0.0f };

    float t_intensity;
    float t_radius;
    float t_minGradient, t_maxGradient;

    glm::vec4 color;
};

class TransferFunction2DWidget {
public:
    TransferFunction2DWidget(const volume::Volume& volume, const volume::GradientVolume& gradient);

    void draw();
    void updateRenderConfig(render::RenderConfig& renderConfig);

private:
    int currentlySelectedTF;
    std::vector<TFunction> m_tfunctions;

    float m_intensity, m_maxIntensity;
    float m_radius;
    glm::vec4 m_color;

    int m_interactingPoint;
    GLuint m_histogramImg;

    bool updated;
};
}
