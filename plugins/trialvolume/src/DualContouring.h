#pragma once

#include <vector>

#include <Eigen/Dense>

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
    struct QEF {
        Eigen::MatrixX3f edge_surface_points;
        Eigen::MatrixX3f normals;
        Eigen::VectorXf b;
    };
    
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
        return iso_level_slot_.IsDirty();
    }

    /**
     * Reset the dirty flags.
     */
    inline void resetDirtyFlags() {
        // TODO expand if more parameters are added
        iso_level_slot_.ResetDirty();
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

    /**
     * Compute the intersection point on the edge for 2 given relative heights.
     */
    float zeroCrossingLocation(float level1, float level2);

    /**
     * Compute the approximate normal at a point in the volume utilizing the forward difference.
     */
    Eigen::Vector3f forwardDiffNormalAtVolumePoint(const size_t x, const size_t y, const size_t z, const geocalls::VolumetricDataCall &volumeDataCall); 

    /**
     * Find the best position in the cell to place the vertex at.
     */
    Eigen::Vector3f findBestVertex(const size_t x, const size_t y, const size_t z, const geocalls::VolumetricDataCall &volumeDataCall, QEF& qef);

    /** The slot specifying the iso surface level */
    core::param::ParamSlot iso_level_slot_;

    /** The absolute iso surface level */
    float absolute_iso_level_;

    /** The slot for requesting data */
    core::CalleeSlot out_triangle_surface_slot_;

    /** The slot accessing the original volume data */
    core::CallerSlot in_volume_data_slot_;

    /** The data update hash */
    std::size_t data_hash_ = 0;

    /** The hash of when the input data was last read */
    std::size_t in_data_hash_;

    /** The bounding box */
    vislib::math::Cuboid<float> bbox_;

    /** Last time the data was updated */
    unsigned int time_ = 0;

    /** The mesh data */
    std::shared_ptr<std::vector<float>> vertex_buffer_, normal_buffer_;
    std::shared_ptr<std::vector<unsigned int>> index_buffer_;
};

} // namespace megamol::trialvolume
