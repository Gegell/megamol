/**
 * MegaMol
 * Copyright (c) 2016, MegaMol Dev Team
 * All rights reserved.
 */

#pragma once

#include <vector>

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
class ParticleToVolume : public core::Module {
public:
    enum SplattingMethod : int {
        SPLAT_METHOD_KERNEL = 0,
        SPLAT_METHOD_NATURAL_NEIGHBOR = 1
    };

    enum KernelType : int {
        KERNEL_TYPE_NEAREST = 0,
        KERNEL_TYPE_BUMP = 1
    };

    enum KernelMetric : int {
        KERNEL_METRIC_EUCLIDEAN = 0,
        KERNEL_METRIC_MANHATTAN = 1,
        KERNEL_METRIC_CHEBYSHEV = 2
    };

    enum KernelBoundary : int {
        KERNEL_BOUNDARY_CLIP = 0,
        KERNEL_BOUNDARY_CLAMP = 1,
        KERNEL_BOUNDARY_WRAP = 2
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
     * Gets the velocity of the source.
     *
     * @param caller The calling call.
     * @return 'true' on success, 'false' on failure.
     */
    bool getVelocityCallback(core::Call& caller);

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

    bool assertData(geocalls::VolumetricDataCall &caller);

    bool createVolume(geocalls::MultiParticleDataCall* caller);

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
        return splatting_method_slot_.IsDirty()
            || voxel_size_slot_.IsDirty()
            || kernel_type_slot_.IsDirty()
            || kernel_metric_slot_.IsDirty()
            || kernel_radius_slot_.IsDirty()
            || kernel_boundary_slot_.IsDirty();
    }

    /**
     * Reset the dirty flags.
     */
    inline void resetDirtyFlags() {
        // TODO expand if more parameters are added
        splatting_method_slot_.ResetDirty();
        voxel_size_slot_.ResetDirty();
        kernel_type_slot_.ResetDirty();
        kernel_metric_slot_.ResetDirty();
        kernel_radius_slot_.ResetDirty();
        kernel_boundary_slot_.ResetDirty();
    }

    /** The slot specifying the splatting method */
    core::param::ParamSlot splatting_method_slot_;

    /** The slot specifying the kernel type */
    core::param::ParamSlot kernel_type_slot_;

    /** The slot specifying the kernel metric */
    core::param::ParamSlot kernel_metric_slot_;

    /** The slot specifying the kernel radius */
    core::param::ParamSlot kernel_radius_slot_;

    /** The slot specifying the kernal boundary handling */
    core::param::ParamSlot kernel_boundary_slot_;

    /** The slot for a single volume voxel sidelength */
    core::param::ParamSlot voxel_size_slot_;

    /** The slot for requesting data */
    core::CalleeSlot out_density_slot_;

    /** The slot for the velocity volume data output */
    core::CalleeSlot out_velocity_slot_;

    /** The slot accessing the original particle data */
    core::CallerSlot in_particle_data_slot_;

    /** The data update hash */
    std::size_t data_hash_ = 0;

    /** The last hash of the particle data */
    std::size_t in_data_hash_;

    /** The volume density data */
    std::vector<float> density_;

    /** The volume velocity data */
    std::vector<float> velocity_;

    /** The minimum value of the volume */
    float min_value_ = 0.0f;

    /** The maximum value of the volume */
    float max_value_ = 0.0f;

    /** The mimimum velocities of the volume */
    std::array<float, 3> min_velocity_ = {0.0f, 0.0f, 0.0f};

    /** The maximum velocities of the volume */
    std::array<float, 3> max_velocity_ = {0.0f, 0.0f, 0.0f};

    /** The number of voxels in the x direction */
    size_t x_cells_;

    /** The number of voxels in the y direction */
    size_t y_cells_;

    /** The number of voxels in the z direction */
    size_t z_cells_;

    /** The volume density metadata */
    megamol::geocalls::VolumetricDataCall::Metadata metadata_density_;

    /** The volume velocity metadata */
    megamol::geocalls::VolumetricDataCall::Metadata metadata_velocity_;

    /** Last time the data was updated */
    unsigned int time_ = 0;
};

} // namespace megamol::trialvolume
