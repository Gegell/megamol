#include "trialvolume/GraphCall.h"

using namespace megamol::trialvolume;

std::shared_ptr<GraphData> GraphCall::GetGraph() const {
    return graph_;
}

void GraphCall::SetGraph(const std::shared_ptr<GraphData>& graph) {
    this->graph_ = graph;
}
