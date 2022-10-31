#pragma once

#include "trialvolume/TrackingData.h"

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

namespace megamol::trialvolume {

class VolumeClusterTracking : public core::Module {

public:
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

    /**
     * Global unique ID that can e.g. be used for hash calculation.
     *
     * @return Unique ID
     */
    static inline size_t GUID() {
        return 0x72bfcc51ab0eaaf9ull;
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
    // bool generateDotFile(bool silent=false);

    /** Generate the corresponding .tsv files for the current tracks */
    // bool generateTsvFiles(bool silent=false);

    /** The slot for the cluster track call, receiving request from upper node */
    core::CalleeSlot out_poll_cluster_track_slot_;

    /** The slot for the cluster track call, pushing updates to lower nodes */
    core::CallerSlot out_push_cluster_track_slot_;

    /** The slot for the volume cluster id call */
    core::CallerSlot in_cluster_id_slot_;

    /** The slot for the volume velocity call */
    core::CallerSlot in_velocity_slot_;

    /** The slot for the timestamp information */
    core::CallerSlot in_timestamp_slot_;

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

    /** Store the cluster metadata for each frame */
    ClusterGraph graph_data_;

    /** Store the timestep of every frame */
    std::vector<float> frame_timesteps_;

    /** The hash of the last cluster call */
    size_t hash_;

    /** Store whether we already have one computed result */
    bool has_result_;
};

} // namespace megamol::trialvolume
