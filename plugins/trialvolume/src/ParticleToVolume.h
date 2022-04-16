/**
 * MegaMol
 * Copyright (c) 2016, MegaMol Dev Team
 * All rights reserved.
 */

#ifndef MEGAMOL_TRIALVOLUME_PARTICLETOVOLUME_H
#define MEGAMOL_TRIALVOLUME_PARTICLETOVOLUME_H

#include <vector>

#include "geometry_calls/MultiParticleDataCall.h"
#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "vislib/String.h"
#include "vislib/math/Cuboid.h"

namespace megamol::trialvolume {

/**
 * Loader module for splatting particle data into a volume
 */
class ParticleToVolume : public core::Module {
public:
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
     * Gets a human readable description of the module.
     *
     * TUTORIAL: Mandatory method for every module or call that returns a description.
     *
     * @return A human readable description of the module.
     */
    static const char* Description() {
        return "Splatts particle data into a volume.";
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

    /**
     * Gets the data from the source.
     *
     * TUTORIAL: This method computes the final data and writes them to the calling call.
     * Data that was written by getExtentCallback should be written again, just in case...
     *
     * @param caller The calling call.
     * @return 'true' on success, 'false' on failure.
     */
    bool getDataCallback(core::Call& caller);

    /**
     * Gets the data extents from the source.
     *
     * TUTORIAL: This method computes the extents of the data set, namely the bounding box and other relevant values,
     * and writes them into the calling call.
     *
     * @param caller The calling call.
     * @return 'true' on success, 'false' on failure.
     */
    bool getExtentCallback(core::Call& caller);

    bool dummyCallback(core::Call& caller);

    bool createVolume(geocalls::MultiParticleDataCall* caller);

    /**
     * Implementation of 'Release'.
     */
    void release() override;

    /**
     * Check if any of the parameters has changed.
     */
    inline bool anythingDirty() const {
        // TODO expand if more parameters are added
        return this->voxelSizeSlot.IsDirty();
    }

    /**
     * Reset the dirty flags.
     */
    inline void resetDirtyFlags() {
        // TODO expand if more parameters are added
        this->voxelSizeSlot.ResetDirty();
    }

    /** The slot for a single volume voxel sidelength */
    core::param::ParamSlot voxelSizeSlot;

    /** The slot for requesting data */
    core::CalleeSlot outDataSlot;

    /** The slot accessing the original particle data */
    core::CallerSlot inParticleDataSlot;

    /** The data update hash */
    std::size_t dataHash = 0;

    /** The last hash of the particle data */
    std::size_t inDataHash;

    /** The bounding box */
    vislib::math::Cuboid<float> bbox;

    /** Pointer to the sphere parameters. */
    std::vector<float> spheres;

    /** The volume data */
    std::vector<float> volume;

    /** The volume metadata */
    megamol::geocalls::VolumetricDataCall::Metadata metadata;

    /** Last time the data was updated */
    unsigned int time = 0;
};

} // namespace megamol::trialvolume

#endif // MEGAMOL_TRIALVOLUME_PARTICLETOVOLUME_H
