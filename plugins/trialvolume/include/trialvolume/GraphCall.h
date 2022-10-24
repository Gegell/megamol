#pragma once

#include "mmstd/data/AbstractGetDataCall.h"
#include "trialvolume/GraphData.h"

namespace megamol::trialvolume {

class GraphCall : public core::AbstractGetDataCall {
public:
    typedef core::factories::CallAutoDescription<GraphCall> graph_description;

    /** The name of the objects of this description. */
    static const char* ClassName() {
        return "GraphCall";
    }

    /** The description of the objects of this type. */
    static const char* Description() {
        return "Call transporting a graph data structure.";
    }

    /** The number of functions used for this call. */
    static const size_t FunctionCount() {
        return 1;
    }

    /** The names of the functions used for this call. */
    static const char* FunctionName(size_t index) {
        switch (index) {
        case 0:
            return "getData";
        default:
            return nullptr;
        }
    }

    /** Get graph */
    std::shared_ptr<GraphData> GetGraph() const;

    /** Set graph */
    void SetGraph(const std::shared_ptr<GraphData>& graph);

private:
    std::shared_ptr<GraphData> graph_;
};


/** Description class typedef */
typedef core::factories::CallAutoDescription<GraphCall> GraphCallDescription;

} // namespace megamol::trialvolume
