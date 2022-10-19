#include "ParticleToVolumeGL.h"

#include "OpenGL_Context.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore_gl/utility/ShaderFactory.h"

#include <GL/glu.h>
#include <glowl/glowl.h>

using namespace megamol::trialvolume_gl;
using megamol::core::utility::log::Log;

#define checkGLError                                                                                   \
    {                                                                                                  \
        GLenum errCode = glGetError();                                                                 \
        if (errCode != GL_NO_ERROR)                                                                    \
            std::cout << "Error in line " << __LINE__ << ": " << gluErrorString(errCode) << std::endl; \
    }

bool ParticleToVolumeGL::create(void) {
    auto const& ogl_ctx = frontend_resources.get<frontend_resources::OpenGL_Context>();
    if (!ogl_ctx.isVersionGEQ(4, 3)) {
        Log::DefaultLog.WriteError("OpenGL version 4.3 or higher required.");
        return false;
    }
    if (!ogl_ctx.isExtAvailable("GL_NV_shader_atomic_float")) {
        Log::DefaultLog.WriteError("GL_NV_shader_atomic_float extension required.");
        return false;
    }

    auto shader_options = msf::ShaderFactoryOptionsOpenGL(GetCoreInstance()->GetShaderPaths());
    // HACK, should not be necessary, e.g. should be set automatically, but is not
    shader_options.addDefinition("GL_NV_shader_atomic_float");

    try {
        calc_volume_program_ = core::utility::make_glowl_shader(
            "trialvol_splat_volume", shader_options, "trialvolume_gl/splat_volume.comp.glsl");
    } catch (std::exception const& e) {
        Log::DefaultLog.WriteError("Unable to compile compute shader for ParticleToVolumeGL: %s", e.what());
        return false;
    }

    glGetProgramiv(calc_volume_program_->getHandle(), GL_COMPUTE_WORK_GROUP_SIZE, splat_workgroup_size_);

    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 0, &max_workgroup_count_[0]);
    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 1, &max_workgroup_count_[1]);
    glGetIntegeri_v(GL_MAX_COMPUTE_WORK_GROUP_COUNT, 2, &max_workgroup_count_[2]);

    glGenBuffers(1, &volume_density_buffer_);
    glGenBuffers(1, &volume_velocity_buffer_);
    glGenBuffers(1, &particle_position_buffer_);
    glGenBuffers(1, &particle_velocity_buffer_);

    return true;
}

void ParticleToVolumeGL::release(void) {
    // TODO release any data here
    glDeleteBuffers(1, &volume_density_buffer_);
    glDeleteBuffers(1, &volume_velocity_buffer_);
    glDeleteBuffers(1, &particle_position_buffer_);
    glDeleteBuffers(1, &particle_velocity_buffer_);
}

std::vector<std::string> ParticleToVolumeGL::requested_lifetime_resources() {
    std::vector<std::string> resources = BaseParticleToVolume::requested_lifetime_resources();
    resources.emplace_back("OpenGL_Context");
    return resources;
}

ParticleToVolumeGL::ParticleToVolumeGL() : BaseParticleToVolume() {}

ParticleToVolumeGL::~ParticleToVolumeGL() {
    Release();
}

void ParticleToVolumeGL::bindOutputBuffers() {
    if (buffer_dimensions_ == glm::uvec3(x_cells_, y_cells_, z_cells_))
        return;

    buffer_dimensions_ = glm::uvec3(x_cells_, y_cells_, z_cells_);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_density_buffer_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * density_.size(), density_.data(), GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, volume_density_buffer_);
    checkGLError;

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_velocity_buffer_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * velocity_.size(), velocity_.data(), GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 3, volume_velocity_buffer_);
    checkGLError;

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, 0);
    checkGLError;
}

void ParticleToVolumeGL::bindInputBuffer(geocalls::MultiParticleDataCall* caller) {
    num_particles_ = 0;
    for (size_t i = 0; i < caller->GetParticleListCount(); ++i) {
        auto const& pl = caller->AccessParticles(i);
        num_particles_ += pl.GetCount();
    }

    if (num_particles_ == 0)
        return;

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, particle_position_buffer_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * num_particles_ * 3, nullptr, GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particle_position_buffer_);
    checkGLError;

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, particle_velocity_buffer_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(float) * num_particles_ * 3, nullptr, GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 1, particle_velocity_buffer_);
    checkGLError;


    size_t offset = 0;
    std::vector<float> buffer;
    for (size_t i = 0; i < caller->GetParticleListCount(); ++i) {
        auto const& particles = caller->AccessParticles(i);

        const float* pos_ptr;
        if (particles.GetVertexDataType() == geocalls::MultiParticleDataCall::Particles::VERTDATA_FLOAT_XYZ) {
            pos_ptr = static_cast<const float*>(particles.GetVertexData());
        } else {
            buffer.resize(particles.GetCount() * 3);
            pos_ptr = buffer.data();
            auto const& ps = particles.GetParticleStore();
            auto const& xAcc = ps.GetXAcc();
            auto const& yAcc = ps.GetYAcc();
            auto const& zAcc = ps.GetZAcc();
            for (size_t j = 0; j < particles.GetCount(); ++j) {
                buffer[j * 3 + 0] = xAcc->Get_f(j);
                buffer[j * 3 + 1] = yAcc->Get_f(j);
                buffer[j * 3 + 2] = zAcc->Get_f(j);
            }
        }
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, particle_position_buffer_);
        checkGLError;
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset, sizeof(float) * particles.GetCount() * 3, pos_ptr);
        checkGLError;

        const float* vel_ptr;
        if (particles.GetDirDataType() == geocalls::MultiParticleDataCall::Particles::DIRDATA_FLOAT_XYZ) {
            vel_ptr = static_cast<const float*>(particles.GetDirData());
        } else {
            buffer.resize(particles.GetCount() * 3);
            vel_ptr = buffer.data();
            auto const& ps = particles.GetParticleStore();
            auto const& dxAcc = ps.GetDXAcc();
            auto const& dyAcc = ps.GetDYAcc();
            auto const& dzAcc = ps.GetDZAcc();
            for (size_t j = 0; j < particles.GetCount(); ++j) {
                buffer[j * 3 + 0] = dxAcc->Get_f(j);
                buffer[j * 3 + 1] = dyAcc->Get_f(j);
                buffer[j * 3 + 2] = dzAcc->Get_f(j);
            }
        }
        glBindBuffer(GL_SHADER_STORAGE_BUFFER, particle_velocity_buffer_);
        checkGLError;
        glBufferSubData(GL_SHADER_STORAGE_BUFFER, offset, sizeof(float) * particles.GetCount() * 3, vel_ptr);
        checkGLError;

        offset += sizeof(float) * particles.GetCount() * 3;
    }
}

bool ParticleToVolumeGL::computeVolume(geocalls::MultiParticleDataCall* caller) {
    calc_volume_program_->use();

    bindInputBuffer(caller);
    bindOutputBuffers();

    auto const& bbox = caller->AccessBoundingBoxes().ObjectSpaceBBox();

    calc_volume_program_->setUniform("numCells", buffer_dimensions_);
    calc_volume_program_->setUniform("numParticles", static_cast<GLuint>(num_particles_));
    calc_volume_program_->setUniform("bboxMin", glm::vec3(bbox.Left(), bbox.Bottom(), bbox.Back()));
    calc_volume_program_->setUniform("bboxMax", glm::vec3(bbox.Right(), bbox.Top(), bbox.Front()));
    checkGLError;

    calc_volume_program_->setUniform("kernel.type",     static_cast<GLuint>(kernel_type_slot_.Param<core::param::EnumParam>()->Value()));
    calc_volume_program_->setUniform("kernel.boundary", static_cast<GLuint>(kernel_boundary_slot_.Param<core::param::EnumParam>()->Value()));
    calc_volume_program_->setUniform("kernel.metric",   static_cast<GLuint>(kernel_metric_slot_.Param<core::param::EnumParam>()->Value()));
    calc_volume_program_->setUniform("kernel.radius", kernel_radius_slot_.Param<core::param::FloatParam>()->Value());
    checkGLError;

    // Actually start the computation
    // TODO adjust the workgroup size, to not be just 1D
    glm::uvec3 workgroup_size(splat_workgroup_size_[0], splat_workgroup_size_[1], splat_workgroup_size_[2]);
    glm::uvec3 num_workgroups((num_particles_ + workgroup_size.x - 1) / workgroup_size.x, 1, 1);
    glDispatchCompute(num_workgroups.x, num_workgroups.y, num_workgroups.z);
    checkGLError;

    glMemoryBarrier(GL_SHADER_STORAGE_BARRIER_BIT);
    checkGLError;
    glFinish();

    // Load the data back from the GPU into RAM
    // glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_density_buffer_);
    // checkGLError;
    // auto* data = static_cast<float*>(glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY));
    // checkGLError;
    // std::memcpy(density_.data(), data, sizeof(float) * density_.size());
    // glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    // checkGLError;

    // glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_velocity_buffer_);
    // checkGLError;
    // data = static_cast<float*>(glMapBuffer(GL_SHADER_STORAGE_BUFFER, GL_READ_ONLY));
    // checkGLError;
    // std::memcpy(velocity_.data(), data, sizeof(float) * velocity_.size());
    // glUnmapBuffer(GL_SHADER_STORAGE_BUFFER);
    // checkGLError;
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_density_buffer_);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(float) * density_.size(), density_.data());
    checkGLError;
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_velocity_buffer_);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(float) * velocity_.size(), velocity_.data());
    checkGLError;


    return true;
}
