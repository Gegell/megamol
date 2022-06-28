#pragma once

#include "MeshSegmentation.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"

#include "datatools/table/TableDataCall.h"
#include "mesh/MeshDataCall.h"

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
    SegmentMetadata computeMetrics(
        const MeshSegmentation::Segment& segment,
        const int id,
        const std::vector<float>& vertices,
        const std::vector<unsigned int>& indices);

    /** Computes the new metrics if new data is available. */
    bool recalculateMetrics();

    /** Refresh the input, if necessary. */
    bool refreshInput();

    /** Check input and recalculate metrics, if necessary. */
    bool ensureFreshMetrics();

    /** Callback for the tabular data call. */
    bool analyzeSegmentsCallback(core::Call& call);

    /** Callback for the mesh data call. */
    bool meshDataCallCallback(core::Call& call);

    /** Callback for the mesh metadata call. */
    bool meshMetadataCallCallback(core::Call& call);

    /** Callback for updating the transfer function. */
    bool transferFunctionCallback(core::param::ParamSlot& param);

    /** Callback for getting the hash. */
    bool getHashCallback(core::Call& call);

    /** The slot for the mesh data. */
    core::CallerSlot mesh_slot_;

    struct input_t {
        std::shared_ptr<std::vector<megamol::trialvolume::MeshSegmentation::Segment>> segments;
        std::shared_ptr<std::vector<float>> vertices;
        std::shared_ptr<std::vector<unsigned int>> indices;

        size_t hash;
    } input_data_;

    /** The slot for the data asociated with the mesh. */
    core::CalleeSlot mesh_data_slot_;

    /** Tabular data output. */
    core::CalleeSlot tabular_output_slot_;

    /** Transfer function for the data sets. */
    core::param::ParamSlot transfer_function_slot_;

    /** Store the segmentation metrics. */
    std::vector<SegmentMetadata> segment_metrics_;

    /** The hash of the last computation. */
    size_t hash_;

    /** The table columns */
    std::vector<datatools::table::TableDataCall::ColumnInfo> table_columns_;

    /** The table content */
    std::vector<float> table_content_;

    /** Output data sets. */
    std::array<std::shared_ptr<mesh::MeshDataCall::data_set>, 5> output_data_sets_;
};

} // namespace megamol::trialvolume