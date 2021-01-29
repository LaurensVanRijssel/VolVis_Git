#pragma once
#include "render/render_config.h"
#include "volume/gradient_volume.h"
#include "volume/volume.h"
#include <GL/glew.h> // Include before glfw3
#include <glm/vec2.hpp>
#include <glm/vec4.hpp>

namespace ui {

//Define a structure for the transfer function
//This allows for multiple transfer functions to be defined.
//This struct is the exact same as the one defined in render_config.h to allow them to share data easily.
struct TFunction {
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
