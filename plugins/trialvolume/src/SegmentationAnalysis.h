#pragma once

#include "MeshSegmentation.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "vislib/math/Cuboid.h"

namespace megamol::trialvolume {

class SegmentationAnalysis : public core::Module {
public:
    struct SegmentMetadata {
        int id;

        float centroid[3];
        vislib::math::Cuboid<float> bounds;

        float volume;
        float surface_area;

        float singular_vals[3];
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
    SegmentMetadata computeMetrics(const MeshSegmentation::Segment& segment, const int id);

    /** Callback for the manual triggered analysis. */
    bool buttonPressedCallback(core::param::ParamSlot& slot);

    /** The slot for the mesh data. */
    core::CallerSlot mesh_slot_;

    /** Button to trigger the analysis. */
    core::param::ParamSlot button_slot_;
};

} // namespace megamol::trialvolume