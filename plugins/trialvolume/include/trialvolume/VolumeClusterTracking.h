#pragma once

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

namespace megamol::trialvolume {

class VolumeClusterTracking : public core::Module {

public:
    /** Store general information of the clusters */
    struct ClusterMetadata_t {
        unsigned int local_time_cluster_id;
        unsigned int frame_id;
        std::map<size_t, size_t> parents;
        float total_mass;

        vislib::math::Cuboid<float> bounding_box;
        vislib::math::Vector<float, 3> center_of_mass;
        vislib::math::Vector<float, 3> velocity;
    };

    /** Answer the name of the objects of this description. */
    static const char* ClassName() {
        return "VolumeClusterTracking";
    }

    /** Gets a human readable description of the module. */
    static const char* Description() {
        return "Tracks clusters over time with volumetric data.";
    }

    /** Answers whether this module is available on the current system. */
    static bool IsAvailable() {
        return true;
    }

    /** Constructor */
    VolumeClusterTracking();

    /** Destructor */
    ~VolumeClusterTracking() override;

private:
    /** Implementation of 'Create'. */
    bool create() override;

    /** Implementation of 'Release'. */
    void release() override;

    /** Check if any parameter is dirty */
    inline bool anythingDirty() const {
        return false;
    }

    /** Resets the dirty flags */
    inline void resetDirtyFlags() {}

    /** Retrieves the cluster track data */
    bool getDataCallback(core::Call& call);

    /** Retrieves the cluster track extent */
    bool getExtentCallback(core::Call& call);

    /** The callback for the manual start button */
    bool buttonCallback(core::param::ParamSlot& slot);

    /** Computes the newest tracks, if any are available */
    void computeTracks();

    /** Generate a corresponding .dot file for the current tracks */
    bool generateDotFile(bool silent=false);

    /** Generate the corresponding .tsv files for the current tracks */
    bool generateTsvFiles(bool silent=false);

    /** The slot for the cluster track call */
    core::CalleeSlot out_cluster_track_slot_;

    /** The slot for the volume cluster id call */
    core::CallerSlot in_cluster_id_slot_;

    /** The slot for the volume velocity call */
    core::CallerSlot in_velocity_slot_;

    /** The slot for the time step size */
    core::param::ParamSlot time_step_size_;

    /** The slot for the manual start button */
    core::param::ParamSlot start_button_;

    /** The minimal amount of connections needed to connect two clusters */
    core::param::ParamSlot min_connection_count_;

    /** The frame range and step size for the cluster track call */
    core::param::ParamSlot frame_start_param_;
    core::param::ParamSlot frame_end_param_;
    core::param::ParamSlot frame_step_param_;
    core::param::ParamSlot frame_range_limit_param_;

    /** The file name for the .dot file */
    core::param::ParamSlot dot_file_name_;

    /** The name for the .tsv directory */
    core::param::ParamSlot tsv_directory_;

    /** Store the cluster metadata for each frame */
    std::vector<std::vector<ClusterMetadata_t>> cluster_metadata_;

    /** The hash of the last cluster call */
    size_t hash_;
};

} // namespace megamol::trialvolume
