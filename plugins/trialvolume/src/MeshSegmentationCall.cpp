#include "MeshSegmentationCall.h"

namespace megamol::trialvolume {

std::shared_ptr<std::vector<MeshSegmentation::Segment>> MeshSegmentationCall::GetSegments() {
    return segments_;
}

void MeshSegmentationCall::SetSegments(std::shared_ptr<std::vector<MeshSegmentation::Segment>> &segmentation) {
    segments_ = segmentation;
}

}