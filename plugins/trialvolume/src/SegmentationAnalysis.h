#pragma once

#include "MeshSegmentation.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"

#include "datatools/table/TableDataCall.h"

namespace megamol::trialvolume {

class SegmentationAnalysis : public core::Module {
public:
    struct SegmentMetadata {
        int id;
        int num_vertices;
        int num_triangles;

        float centroid[3];
        float extents[3];

        float volume;
        float surface_area;
        float sphericity;

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

    /** Computes the new metrics if new data is available. */
    bool recalculateMetrics();

    /** Callback for the manual triggered analysis. */
    bool buttonPressedCallback(core::param::ParamSlot& slot);

    /** Callback for the tabular data call. */
    bool analyzeSegmentsCallback(core::Call& call);

    /** Callback for getting the hash. */
    bool getHashCallback(core::Call& call);

    /** The slot for the mesh data. */
    core::CallerSlot mesh_slot_;

    /** Tabular data output. */
    core::CalleeSlot tabular_output_slot_;

    /** Button to trigger the analysis. */
    core::param::ParamSlot button_slot_;

    /** Store the segmentation metrics. */
    std::vector<SegmentMetadata> segment_metrics_;

    /** The hash of the last computation. */
    size_t hash_;

    /** The table columns */
    std::vector<datatools::table::TableDataCall::ColumnInfo> table_columns_;

    /** The table content */
    std::vector<float> table_content_;
};

} // namespace megamol::trialvolume