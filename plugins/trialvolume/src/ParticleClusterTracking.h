#pragma once

#include <map>
#include <vector>

#include "geometry_calls/MultiParticleDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "vislib/math/Cuboid.h"

namespace megamol::trialvolume {

class ParticleClusterTracking : public core::Module {
public:
    typedef size_t track_unique_id_t;

    /** Store general information of the clusters */
    struct ClusterMetadata_t {
        unsigned int local_time_cluster_id;
        unsigned int frame_id;
        std::map<size_t, size_t> parents;
        unsigned int num_particles;

        vislib::math::Cuboid<float> bounding_box;
    };

    /** Answer the name of the objects of this description. */
    static const char* ClassName() {
        return "ParticleClusterTracking";
    }

    /** Gets a human readable description of the module. */
    static const char* Description() {
        return "Tracks the particle IColor clusters over time.";
    }

    /** Answers whether this module is available on the current system. */
    static bool IsAvailable() {
        return true;
    }

    /** Constructor */
    ParticleClusterTracking();

    /** Destructor */
    ~ParticleClusterTracking() override;

private:
    /** Implementation of 'Create'. */
    bool create(void) override;

    /** Implementation of 'Release'. */
    void release(void) override;

    /** Gets the extent of the tracked clusters */
    bool getExtentCallback(core::Call& call);

    /** Gets the data of the tracked clusters */
    bool getDataCallback(core::Call& call);

    /** Start via button press callback */
    bool buttonCallback(core::param::ParamSlot& param);

    /** Computes the newest tracks, if any are available */
    void computeTracks(void);

    /** Generate a corresponding .dot file for the current tracks */
    bool generateDotFile(bool silent=false);

    /** The slot for the cluster track call */
    megamol::core::CalleeSlot out_cluster_track_slot_;

    /** The slot for the particle cluster call, for one time step */
    megamol::core::CallerSlot in_cluster_slot_;

    /** The param slot for the manual start button */
    megamol::core::param::ParamSlot start_button_;

    /** The minimal amount of particles needed to connect two clusters */
    megamol::core::param::ParamSlot min_connection_count_;

    /** The frame range and step size for the cluster track call */
    megamol::core::param::ParamSlot frame_start_param_;
    megamol::core::param::ParamSlot frame_end_param_;
    megamol::core::param::ParamSlot frame_step_param_;
    megamol::core::param::ParamSlot frame_range_limit_param_;

    /** The file name for the .dot file */
    megamol::core::param::ParamSlot dot_file_name_;

    /** Store a the cluster metadata for each time slice */
    std::vector<std::vector<ClusterMetadata_t>> cluster_metadata_;

    /** The hash of the last cluster call */
    size_t hash_;
};
} // namespace megamol::trialvolume
