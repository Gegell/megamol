#pragma once

#include "trialvolume/GraphData.h"
#include "vislib/math/Cuboid.h"
#include "vislib/math/Vector.h"

#include <string>
#include <vector>

namespace megamol::trialvolume {

class ClusterInfo : public AbstractNodeData {
public:
    void writeDotAttributes(std::ostream& os) const override {
        os << "label=\"" << frame << "_" << frame_local_id << " (" << total_mass << ")\", __bounds=\"["
           << bounding_box.Left() << "," << bounding_box.Bottom() << "," << bounding_box.Back() << ","
           << bounding_box.Right() << "," << bounding_box.Top() << "," << bounding_box.Front() << "]\""
           << ", __center_of_mass=\"[" << center_of_mass.X() << "," << center_of_mass.Y() << "," << center_of_mass.Z()
           << "]\""
           << ", __velocity=\"[" << velocity.X() << "," << velocity.Y() << "," << velocity.Z() << "]\""
           << ", __frame=" << frame << ", __local_id=" << frame_local_id << ", __total_mass=" << total_mass;
    }

    void writeTsvHeaders(std::ostream& os) const override {
        // TODO actually use the base class
        os << "graph_id\tframe_id\tlocal_id\ttotal_mass\tcenter_of_mass\tvelocity\tbbox";
    }

    void writeTsv(std::ostream& os) const override {
        os << frame << "_" << frame_local_id << "\t" << frame << "\t" << frame_local_id << "\t" << total_mass << "\t"
           << center_of_mass.X() << "," << center_of_mass.Y() << "," << center_of_mass.Z() << "\t" << velocity.X()
           << "," << velocity.Y() << "," << velocity.Z() << "\t" << bounding_box.Left() << "," << bounding_box.Bottom()
           << "," << bounding_box.Back() << "," << bounding_box.Right() << "," << bounding_box.Top() << ","
           << bounding_box.Front();
    }

    std::string getId() const override {
        return std::to_string(frame) + "_" + std::to_string(frame_local_id);
    }

    float total_mass;

    vislib::math::Cuboid<float> bounding_box;
    vislib::math::Vector<float, 3> center_of_mass;
    vislib::math::Vector<float, 3> velocity;

    uint32_t frame;
    uint32_t frame_local_id;
};

class ClusterConnection : public AbstractEdgeData {
public:
    void writeDotAttributes(std::ostream& os) const override {
        os << "label=\"" << kept << "\", __kept=" << kept;
    }

    void writeTsvHeaders(std::ostream& os) const override {
        AbstractEdgeData::writeTsvHeaders(os);
        os << "\tkept";
    }

    void writeTsv(std::ostream& os) const override {
        AbstractEdgeData::writeTsv(os);
        os << "\t" << kept;
    }

    float kept;
};

class ClusterGraph : public GraphData, public TsvGenerator {
public:
    ClusterGraph() = default;
    ~ClusterGraph() = default;

    void writeTsvHeaders(std::ostream& os) const override {
        os << "frame_id\ttimestep";
    }

    void writeTsv(std::ostream& os) const override {
        for (auto const& info : frames_) {
            os << info.frame_id << "\t" << info.timestep << std::endl;
        }
    }

    struct frame_info_t {
        uint32_t frame_id;
        float timestep;
    };

    std::vector<frame_info_t> frames_;
};

} // namespace megamol::trialvolume
