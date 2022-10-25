#include "ParticleToVolumeGL.h"

#include "OpenGL_Context.h"
#include "mmcore/CoreInstance.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore_gl/utility/ShaderFactory.h"

#include <GL/glu.h>
#include <glowl/glowl.h>
#include <omp.h>

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

    glGenBuffers(1, &volume_buffer_);
    glGenBuffers(1, &particle_buffer_);

    return true;
}

void ParticleToVolumeGL::release(void) {
    // TODO release any data here
    glDeleteBuffers(1, &volume_buffer_);
    glDeleteBuffers(1, &particle_buffer_);
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

void ParticleToVolumeGL::bindOutputBuffers(std::vector<VoxelData>& output_buffer) {
    if (buffer_dimensions_ == glm::uvec3(x_cells_, y_cells_, z_cells_))
        return;

    buffer_dimensions_ = glm::uvec3(x_cells_, y_cells_, z_cells_);

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_buffer_);
    glBufferData(
        GL_SHADER_STORAGE_BUFFER, sizeof(VoxelData) * output_buffer.size(), output_buffer.data(), GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 2, volume_buffer_);
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

    struct ParticleData {
        glm::vec3 position;
        GLfloat padding3;
        glm::vec3 velocity;
        GLfloat padding7;
    };

    std::vector<ParticleData> buffer;
    for (size_t i = 0; i < caller->GetParticleListCount(); ++i) {
        auto const& particles = caller->AccessParticles(i);
        auto const& ps = particles.GetParticleStore();
        auto const& xAcc = ps.GetXAcc();
        auto const& yAcc = ps.GetYAcc();
        auto const& zAcc = ps.GetZAcc();
        auto const& dxAcc = ps.GetDXAcc();
        auto const& dyAcc = ps.GetDYAcc();
        auto const& dzAcc = ps.GetDZAcc();

        for (size_t j = 0; j < particles.GetCount(); ++j) {
            buffer.push_back({glm::vec3(xAcc->Get_f(j), yAcc->Get_f(j), zAcc->Get_f(j)), 0.0f,
                glm::vec3(dxAcc->Get_f(j), dyAcc->Get_f(j), dzAcc->Get_f(j)), 0.0f});
        }
    }

    glBindBuffer(GL_SHADER_STORAGE_BUFFER, particle_buffer_);
    glBufferData(GL_SHADER_STORAGE_BUFFER, sizeof(ParticleData) * buffer.size(), buffer.data(), GL_DYNAMIC_COPY);
    glBindBufferBase(GL_SHADER_STORAGE_BUFFER, 0, particle_buffer_);
    checkGLError;
}

bool ParticleToVolumeGL::computeVolume(geocalls::MultiParticleDataCall* caller) {
    calc_volume_program_->use();

    std::vector<VoxelData> data_buffer(density_.size());

    auto const start = std::chrono::high_resolution_clock::now();

    bindInputBuffer(caller);
    bindOutputBuffers(data_buffer);

    auto const bindTime = std::chrono::high_resolution_clock::now();

    // Clear the output buffers
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_buffer_);
    glClearBufferData(GL_SHADER_STORAGE_BUFFER, GL_RGBA32F, GL_RGBA, GL_FLOAT, nullptr);
    checkGLError;

    auto const clearTime = std::chrono::high_resolution_clock::now();

    auto const& bbox = caller->AccessBoundingBoxes().ObjectSpaceBBox();

    calc_volume_program_->setUniform("numCells", buffer_dimensions_);
    calc_volume_program_->setUniform("numParticles", static_cast<GLuint>(num_particles_));
    calc_volume_program_->setUniform("bboxMin", glm::vec3(bbox.Left(), bbox.Bottom(), bbox.Back()));
    calc_volume_program_->setUniform("bboxMax", glm::vec3(bbox.Right(), bbox.Top(), bbox.Front()));
    checkGLError;

    calc_volume_program_->setUniform(
        "kernel.type", static_cast<GLuint>(kernel_type_slot_.Param<core::param::EnumParam>()->Value()));
    calc_volume_program_->setUniform(
        "kernel.boundary", static_cast<GLuint>(kernel_boundary_slot_.Param<core::param::EnumParam>()->Value()));
    calc_volume_program_->setUniform(
        "kernel.metric", static_cast<GLuint>(kernel_metric_slot_.Param<core::param::EnumParam>()->Value()));
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

    auto const endTime = std::chrono::high_resolution_clock::now();

    // Load the data back from the GPU into RAM
    glBindBuffer(GL_SHADER_STORAGE_BUFFER, volume_buffer_);
    glGetBufferSubData(GL_SHADER_STORAGE_BUFFER, 0, sizeof(VoxelData) * data_buffer.size(), data_buffer.data());
    checkGLError;
    auto const readTime = std::chrono::high_resolution_clock::now();

#pragma omp parallel for
    for (auto i = 0; i < data_buffer.size(); ++i) {
        auto const& d = data_buffer[i];
        auto const avg_vel = d.density > 0.0f ? d.velocity / d.density : glm::vec3(0.0f);
        density_[i] = d.density;
        velocity_[i] = avg_vel.x;
        velocity_[i + 1] = avg_vel.y;
        velocity_[i + 2] = avg_vel.z;
    }

    auto const normalizeTime = std::chrono::high_resolution_clock::now();

    // Print some timing information
#if TRIALVOLUME_VERBOSE
    Log::DefaultLog.WriteInfo("ParticleToVolumeGL: Bind time: %f ms",
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(bindTime - start).count());
    Log::DefaultLog.WriteInfo("ParticleToVolumeGL: Clear time: %f ms",
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(clearTime - bindTime).count());
    Log::DefaultLog.WriteInfo("ParticleToVolumeGL: Compute time: %f ms",
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(endTime - clearTime).count());
    Log::DefaultLog.WriteInfo("ParticleToVolumeGL: Read time: %f ms",
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(readTime - endTime).count());
    Log::DefaultLog.WriteInfo("ParticleToVolumeGL: Normalize time: %f ms",
        std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(normalizeTime - readTime).count());
#endif

    return true;
}
