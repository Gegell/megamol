#pragma once

#include <vector>
#include <set>

#include "geometry_calls/MultiParticleDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "vislib/math/Cuboid.h"

namespace megamol::trialvolume {

class ParticleClusterTracking : public core::Module {
public:
    typedef size_t track_unique_id_t;

    /** Store general information of the clusters */
    struct ClusterMetadata_t {
        unsigned int local_time_cluster_id;
        unsigned int frame_id;
        std::set<track_unique_id_t> parents;

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

    /** Computes the newest tracks, if any are available */
    void computeTracks(void);

    /** The slot for the cluster track call */
    megamol::core::CalleeSlot out_cluster_track_slot_;

    /** The slot for the particle cluster call, for one time step */
    megamol::core::CallerSlot in_cluster_slot_;

    /** Store a mapping from the track unique id to the metadata of the time slices */
    std::vector<ClusterMetadata_t> cluster_metadata_;
};
} // namespace megamol::trialvolume
