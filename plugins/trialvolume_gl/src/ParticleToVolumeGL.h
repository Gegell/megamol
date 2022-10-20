#pragma once

#include "mmcore_gl/utility/ShaderFactory.h"

#include "trialvolume/BaseParticleToVolume.h"

#include <glowl/glowl.h>

#include <vector>

namespace megamol::trialvolume_gl {

class ParticleToVolumeGL : public trialvolume::BaseParticleToVolume {

public:
    static const char* ClassName() {
        return "ParticleToVolumeGL";
    }

    static bool IsAvailable() {
        return true;
    }

    ParticleToVolumeGL();
    ~ParticleToVolumeGL() override;

    /**
     * Request resources to ask for OpenGL state
     */
    virtual std::vector<std::string> requested_lifetime_resources() override;

private:
    struct VoxelData {
        glm::vec3 velocity;
        GLfloat density;
    };

    /**
     * Implementation of 'Create'.
     *
     * TUTORIAL: The overwritten version of this method gets called right after an object of this class has been
     *instantiated.
     *
     * @return 'true' on success, 'false' otherwise.
     */
    bool create() final;

    /**
     * Implementation of 'Release'.
     */
    void release() final;

    bool computeVolume(geocalls::MultiParticleDataCall* caller);

    void bindOutputBuffers(std::vector<VoxelData>& voxels);

    void bindInputBuffer(geocalls::MultiParticleDataCall* caller);

    std::unique_ptr<glowl::GLSLProgram> calc_volume_program_;

    GLint splat_workgroup_size_[3];
    GLint max_workgroup_count_[3];

    GLuint volume_buffer_;

    GLuint particle_buffer_;

    glm::uvec3 buffer_dimensions_;

    size_t num_particles_;
};

} // namespace megamol::trialvolume_gl
