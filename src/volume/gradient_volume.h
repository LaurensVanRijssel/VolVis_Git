#pragma once
#include "volume.h"
#include <glm/vec3.hpp>
#include <glm/vec4.hpp>
#include <string>
#include <vector>

namespace volume {
struct GradientVoxel {
    glm::vec3 dir;
    float magnitude;
};

class GradientVolume {
public:
    // DO NOT REMOVE
    InterpolationMode interpolationMode { InterpolationMode::NearestNeighbour };

public:
    GradientVolume(const Volume& volume);

    GradientVoxel getGradientVoxel(const glm::vec3& coord) const;
    GradientVoxel getGradientVoxel(int x, int y, int z) const;

    float minMagnitude() const;
    float maxMagnitude() const;
    glm::ivec3 dims() const;

protected:
    GradientVoxel getGradientVoxelNN(const glm::vec3& coord) const;
    GradientVoxel getGradientVoxelLinearInterpolate(const glm::vec3& coord) const;
    static GradientVoxel linearInterpolate(const GradientVoxel& g0, const GradientVoxel& g1, float factor);

protected:
    const glm::ivec3 m_dim;
    const std::vector<GradientVoxel> m_data;
    const float m_minMagnitude, m_maxMagnitude;
};
}
