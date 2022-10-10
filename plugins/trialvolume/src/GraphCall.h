#pragma once

#include "mmstd/data/AbstractGetDataCall.h"

namespace megamol::trialvolume {

class AbstractNodeData {
public:
    size_t id;

    virtual ~AbstractNodeData() = default;

    virtual std::string ToDot() const {
        return std::to_string(id);
    };
};

class AbstractEdgeData {
public:
    std::shared_ptr<AbstractNodeData> source;
    std::shared_ptr<AbstractNodeData> target;

    virtual ~AbstractEdgeData() = default;

    virtual std::string ToDot() const {
        return std::to_string(source->id) + " -> " + std::to_string(target->id);
    };
};

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

    /** Get nodes */
    std::shared_ptr<std::vector<AbstractNodeData>> GetNodes() const;

    /** Set nodes */
    void SetNodes(const std::shared_ptr<std::vector<AbstractNodeData>>& nodes);

    /** Get edges */
    std::shared_ptr<std::vector<AbstractEdgeData>> GetEdges() const;

    /** Set edges */
    void SetEdges(const std::shared_ptr<std::vector<AbstractEdgeData>>& edges);

private:
    std::shared_ptr<std::vector<AbstractNodeData>> nodes_;
    std::shared_ptr<std::vector<AbstractEdgeData>> edges_;
};

} // namespace megamol::trialvolume
