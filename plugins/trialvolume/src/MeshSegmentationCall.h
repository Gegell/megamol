#pragma once

#include "MeshSegmentation.h"

#include "mmcore/AbstractGetDataCall.h"

namespace megamol::trialvolume {

class MeshSegmentationCall : public core::AbstractGetDataCall {
public:
    typedef core::factories::CallAutoDescription<MeshSegmentationCall> segmentation_description;

    /** The name of the objects of this description. */
    static const char* ClassName() {
        return "MeshSegmentationCall";
    }

    /** The description of the objects of this type. */
    static const char* Description() {
        return "Call transporting a segmented mesh.";
    }

    /** The number of functions used for this call. */
    static const size_t FunctionCount() {
        return 1;
    }

    /** The names of the functions used for this call. */
    static const char* FunctionName(size_t index) {
        switch (index) {
        case 0:
            return "getSegmentation";
        default:
            return nullptr;
        }
    }

    /** Get the segmentation data. */
    std::shared_ptr<std::vector<MeshSegmentation::Segment>> GetSegments() const;

    /** Set the segmentation data. */
    void SetSegments(const std::shared_ptr<std::vector<MeshSegmentation::Segment>> &segments);

    /** Get the base vertex data. */
    std::shared_ptr<std::vector<float>> GetBaseVertices() const;

    /** Set the base vertex data. */
    void SetBaseVertices(const std::shared_ptr<std::vector<float>> &base_vertices);

    /** Get the base normal data. */
    std::shared_ptr<std::vector<float>> GetBaseNormals() const;

    /** Set the base normal data. */
    void SetBaseNormals(const std::shared_ptr<std::vector<float>> &base_normals);

    /** Get the base index data. */
    std::shared_ptr<std::vector<unsigned int>> GetBaseIndices() const;

    /** Set the base index data. */
    void SetBaseIndices(const std::shared_ptr<std::vector<unsigned int>> &base_indices);


private:
    /** The segments found in the mesh. */
    std::shared_ptr<std::vector<MeshSegmentation::Segment>> segments_;

    std::shared_ptr<std::vector<float>> base_vertices_;
    std::shared_ptr<std::vector<float>> base_normals_;
    std::shared_ptr<std::vector<unsigned int>> base_indices_;
};

}