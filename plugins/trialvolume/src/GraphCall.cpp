#include "trialvolume/GraphCall.h"

using namespace megamol::trialvolume;

std::shared_ptr<std::vector<AbstractNodeData>> GraphCall::GetNodes() const {
    return nodes_;
}

void GraphCall::SetNodes(const std::shared_ptr<std::vector<AbstractNodeData>>& nodes) {
    this->nodes_ = nodes;
}

std::shared_ptr<std::vector<AbstractEdgeData>> GraphCall::GetEdges() const {
    return edges_;
}

void GraphCall::SetEdges(const std::shared_ptr<std::vector<AbstractEdgeData>>& edges) {
    this->edges_ = edges;
}
