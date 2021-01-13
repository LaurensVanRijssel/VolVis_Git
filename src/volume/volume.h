#pragma once
#include <filesystem>
#include <glm/vec2.hpp>
#include <glm/vec3.hpp>
#include <string>
#include <vector>

namespace volume {

enum class InterpolationMode {
    NearestNeighbour = 0,
    Linear,
    Cubic
};

class Volume {
public:
    // DO NOT REMOVE
    InterpolationMode interpolationMode { InterpolationMode::NearestNeighbour };

public:
    Volume(const std::filesystem::path& file);
    Volume(std::vector<uint16_t> data, const glm::ivec3& dim);

    float minimum() const;
    float maximum() const;
    std::vector<int> histogram() const;
    glm::ivec3 dims() const;
    std::string_view fileName() const;

    float getVoxelInterpolate(const glm::vec3& coord) const;
    float getVoxel(int x, int y, int z) const;

protected:
    float getVoxelNN(const glm::vec3& coord) const;
    static float linearInterpolate(float g0, float g1, float factor);
    float getVoxelLinearInterpolate(const glm::vec3& coord) const;
    static float weight(float x);
    static float cubicInterpolate(float g0, float g1, float g2, float g3, float factor);
    float bicubicInterpolateXY(const glm::vec2& xyCoord, int z) const;
    float getVoxelTriCubicInterpolate(const glm::vec3& coord) const;

private:
    void loadFile(const std::filesystem::path& file);

protected:
    const std::string m_fileName;
    size_t m_elementSize;
    glm::ivec3 m_dim;

    std::vector<uint16_t> m_data;

    float m_minimum, m_maximum;
    std::vector<int> m_histogram;
};
}
