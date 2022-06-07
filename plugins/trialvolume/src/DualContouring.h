#pragma once

#include <vector>

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"
#include "vislib/math/Cuboid.h"

namespace megamol::trialvolume {

/**
 * Loader module for splatting particle data into a volume
 */
class DualContouring : public core::Module {
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
        return "DualContouring";
    }

    /**
     * Gets a human readable description of the module.
     *
     * TUTORIAL: Mandatory method for every module or call that returns a description.
     *
     * @return A human readable description of the module.
     */
    static const char* Description() {
        return "Performs dual contouring on a volume to generate a surface mesh.";
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
    DualContouring();

    /** Destructor */
    ~DualContouring() override;

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
     * Computes the dual contouring of the volume.
     *
     * TUTORIAL: This method computes the final data and writes them to the calling call.
     * Data that was written by getExtentCallback should be written again, just in case...
     *
     * @param caller The calling call.
     * @return 'true' on success, 'false' on failure.
     */
    bool getTriangleSurfaceCallback(core::Call& caller);

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

    /**
     * Implementation of 'Release'.
     */
    void release() override;

    /**
     * Check if any of the parameters has changed.
     */
    inline bool anythingDirty() const {
        // TODO expand if more parameters are added
        return isoLevelSlot.IsDirty();
    }

    /**
     * Reset the dirty flags.
     */
    inline void resetDirtyFlags() {
        // TODO expand if more parameters are added
        isoLevelSlot.ResetDirty();
    }

    /**
     * Re-computes the surface along the iso line.
     */
    bool computeSurface(geocalls::VolumetricDataCall& volumeDataCall);

    // TODO move to seperate util file
    inline bool sameSign(float a, float b) {
        return (a >= 0.0f) == (b >= 0.0f);
    }

    // TODO move to seperate util file
    inline size_t toFlatIndex(size_t x, size_t y, size_t z, const geocalls::VolumetricDataCall::Metadata* metadata) {
        return x + (y + z * (metadata->Resolution[1] - 1)) * (metadata->Resolution[0] - 1);
    }


    /** The slot specifying the iso surface level */
    core::param::ParamSlot isoLevelSlot;

    /** The slot for requesting data */
    core::CalleeSlot outTriangleSurfaceSlot;

    /** The slot accessing the original volume data */
    core::CallerSlot inVolumeDataSlot;

    /** The data update hash */
    std::size_t dataHash = 0;

    /** The hash of when the input data was last read */
    std::size_t inDataHash;

    /** The bounding box */
    vislib::math::Cuboid<float> bbox;

    /** Last time the data was updated */
    unsigned int time = 0;

    /** The mesh data */
    std::shared_ptr<std::vector<float>> vertex_buffer, normal_buffer;
    std::shared_ptr<std::vector<unsigned int>> index_buffer;
};

} // namespace megamol::trialvolume
