#pragma once

#include <iostream>
#include <vector>
#include <string>

namespace megamol::trialvolume {

class DotGenerator {
public:
    DotGenerator() = default;
    virtual ~DotGenerator() = default;

    virtual void writeDot(std::ostream& os) const = 0;
    virtual void writeDotAttributes(std::ostream& os) const {};
};

class TsvGenerator {
public:
    TsvGenerator() = default;
    virtual ~TsvGenerator() = default;

    virtual void writeTsv(std::ostream& os) const = 0;
    virtual void writeTsvHeaders(std::ostream& os) const = 0;
};

class AbstractNodeData : public DotGenerator, public TsvGenerator {
public:
    uint64_t id;

    virtual ~AbstractNodeData() = default;

    virtual std::string getId() const {
        return std::to_string(id);
    }

    virtual void writeDot(std::ostream& os) const {
        os << getId() << " [";
        writeDotAttributes(os);
        os << "];";
    };

    virtual void writeTsv(std::ostream& os) const {
        os << getId();
    };

    virtual void writeTsvHeaders(std::ostream& os) const {
        os << "id";
    };
};

class AbstractEdgeData : public DotGenerator, public TsvGenerator {
public:
    std::shared_ptr<AbstractNodeData> source;
    std::shared_ptr<AbstractNodeData> target;

    virtual ~AbstractEdgeData() = default;

    virtual void writeDot(std::ostream& os) const {
        os << source->getId() << " -> " << target->getId() << " [";
        writeDotAttributes(os);
        os << "];";
    };

    virtual void writeTsv(std::ostream& os) const {
        os << source->getId() << "\t" << target->getId();
    };

    virtual void writeTsvHeaders(std::ostream& os) const {
        os << "from\tto";
    };
};

class GraphData : public DotGenerator {
public:
    std::vector<std::shared_ptr<AbstractNodeData>> nodes_;
    std::vector<std::shared_ptr<AbstractEdgeData>> edges_;

    virtual ~GraphData() = default;

    virtual void writeDot(std::ostream& os) const {
        os << "digraph G {" << std::endl;
        for (auto& node : nodes_) {
            node->writeDot(os);
            os << std::endl;
        }
        for (auto& edge : edges_) {
            edge->writeDot(os);
            os << std::endl;
        }
        os << "}" << std::endl;
    };
};

} // namespace megamol::trialvolume
