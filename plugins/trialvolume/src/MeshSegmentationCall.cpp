#include "trialvolume/MeshSegmentationCall.h"

namespace megamol::trialvolume {

std::shared_ptr<std::vector<MeshSegmentation::Segment>> MeshSegmentationCall::GetSegments() const {
    return segments_;
}

void MeshSegmentationCall::SetSegments(const std::shared_ptr<std::vector<MeshSegmentation::Segment>>& segmentation) {
    segments_ = segmentation;
}

std::shared_ptr<std::vector<float>> MeshSegmentationCall::GetBaseVertices() const {
    return base_vertices_;
}

void MeshSegmentationCall::SetBaseVertices(const std::shared_ptr<std::vector<float>>& base_vertices) {
    base_vertices_ = base_vertices;
}

std::shared_ptr<std::vector<float>> MeshSegmentationCall::GetBaseNormals() const {
    return base_normals_;
}

void MeshSegmentationCall::SetBaseNormals(const std::shared_ptr<std::vector<float>>& base_normals) {
    base_normals_ = base_normals;
}

std::shared_ptr<std::vector<unsigned int>> MeshSegmentationCall::GetBaseIndices() const {
    return base_indices_;
}

void MeshSegmentationCall::SetBaseIndices(const std::shared_ptr<std::vector<unsigned int>>& base_indices) {
    base_indices_ = base_indices;
}
} // namespace megamol::trialvolume
