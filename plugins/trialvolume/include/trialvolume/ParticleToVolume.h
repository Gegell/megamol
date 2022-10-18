/**
 * MegaMol
 * Copyright (c) 2016, MegaMol Dev Team
 * All rights reserved.
 */

#pragma once

#include <vector>

#include "trialvolume/BaseParticleToVolume.h"
#include "geometry_calls/MultiParticleDataCall.h"
#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "vislib/String.h"

namespace megamol::trialvolume {

/**
 * Loader module for splatting particle data into a volume
 */
class ParticleToVolume : public BaseParticleToVolume {
public:
    enum SplattingMethod : int {
        SPLAT_METHOD_KERNEL = 0,
        SPLAT_METHOD_NATURAL_NEIGHBOR = 1
    };

    /**
     * Answer the name of the objects of this description.
     *
     * TUTORIAL: Mandatory method for every module or call that states the name of the class.
     * This name should be unique.
     *
     * @return The name of the objects of this description.
     */
    static const char* ClassName() {
        return "ParticleToVolume";
    }

    /**
     * Answers whether this module is available on the current system.
     *
     * TUTORIAL: Mandatory method for every module.
     *
     * @return 'true' if the module is available, 'false' otherwise.
     */
    static bool IsAvailable() {
        return true;
    }

    /** Constructor */
    ParticleToVolume();

    /** Destructor */
    ~ParticleToVolume() override;

private:
    /**
     * Implementation of 'Create'.
     *
     * TUTORIAL: The overwritten version of this method gets called right after an object of this class has been
     *instantiated.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() override;

    bool computeVolume(geocalls::MultiParticleDataCall* caller);

    bool computeNaturalNeighborhood(geocalls::MultiParticleDataCall* caller);

    bool computeKernel(geocalls::MultiParticleDataCall* caller);

    /**
     * Implementation of 'Release'.
     */
    void release() override;

    /**
     * Check if any of the parameters has changed.
     */
    inline bool anythingDirty() const {
        // TODO expand if more parameters are added
        return BaseParticleToVolume::anythingDirty() || splatting_method_slot_.IsDirty();
    }

    /**
     * Reset the dirty flags.
     */
    inline void resetDirtyFlags() {
        // TODO expand if more parameters are added
        BaseParticleToVolume::resetDirtyFlags();
        splatting_method_slot_.ResetDirty();
    }

    /** The slot specifying the splatting method */
    core::param::ParamSlot splatting_method_slot_;
};

} // namespace megamol::trialvolume
