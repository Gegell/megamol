#pragma once

#include "MeshSegmentation.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "vislib/math/Cuboid.h"

#include <Eigen/Dense>

namespace megamol::trialvolume {

class SegmentationAnalysis : public core::Module {
public:
    struct SegmentMetadata {
        int id;

        Eigen::Vector3f centroid;
        vislib::math::Cuboid<float> bounds;

        float volume;
        float surface_area;

        Eigen::Vector3f singular_vals;
    };

    /** Return the name of this module. */
    static const char* ClassName() {
        return "SegmentationAnalysis";
    }

    /** Return the description of this module. */
    static const char* Description() {
        return "Compute metrics on a separated mesh.";
    }

    /** Answer whether this module is available on the current system. */
    static bool IsAvailable() {
        return true;
    }

    /** Constructor. */
    SegmentationAnalysis();

    /** Destructor. */
    ~SegmentationAnalysis() override;

private:
    /** Gets called every time the module is created. */
    bool create() override;

    /** Gets called every time the module is released. */
    void release() override;

    /**
     * Computes the metrics for the connected segmented mesh call.
     */
    bool computeMetrics(const MeshSegmentation::Segment& segment, SegmentMetadata& metadata);

    /** Callback for the manual triggered analysis. */
    bool buttonPressedCallback(core::param::ParamSlot& slot);

    float signedVolumeOfTriangle(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3);
    float surfaceAreaOfTriangle(Eigen::Vector3f p1, Eigen::Vector3f p2, Eigen::Vector3f p3);

    /** The slot for the mesh data. */
    core::CallerSlot mesh_slot_;

    /** Button to trigger the analysis. */
    core::param::ParamSlot button_slot_;
};

} // namespace megamol::trialvolume