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
    std::shared_ptr<std::vector<MeshSegmentation::Segment>> GetSegments();

    /** Set the segmentation data. */
    void SetSegments(std::shared_ptr<std::vector<MeshSegmentation::Segment>> &segments);

private:
    /** The segments found in the mesh. */
    std::shared_ptr<std::vector<MeshSegmentation::Segment>> segments_;
};

}