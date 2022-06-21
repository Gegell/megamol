#include "MeshSegmentationCall.h"

namespace megamol::trialvolume {

std::shared_ptr<std::vector<MeshSegmentation::Segment>> MeshSegmentationCall::GetSegmentation() {
    return segments_;
}

void MeshSegmentationCall::SetSegmentation(std::shared_ptr<std::vector<MeshSegmentation::Segment>> &segmentation) {
    segments_ = segmentation;
}

}