/**
 * MegaMol
 * Copyright (c) 2018, MegaMol Dev Team
 * All rights reserved.
 */

#include "TrackingGraphRenderer.h"

#include "trialvolume/GraphCall.h"
#include "trialvolume/TrackingData.h"

#include "mmcore/CoreInstance.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore_gl/utility/ShaderFactory.h"
#include "mmstd_gl/renderer/CallRender3DGL.h"
#include "vislib/math/Matrix.h"
#include "vislib/math/ShallowMatrix.h"

#include <GL/glu.h>

#include <map>
#include <vector>

using namespace megamol;
using namespace megamol::trialvolume_gl;
using megamol::core::utility::log::Log;

/*
 * TrackingGraphRenderer::TrackingGraphRenderer
 */
TrackingGraphRenderer::TrackingGraphRenderer()
        : mmstd_gl::Renderer3DModuleGL()
        , in_graph_data_slot_("inData", "The input data slot for the graph data.")
        , line_width_slot_("line width", "Width of the connecting lines") {
    this->in_graph_data_slot_.SetCompatibleCall<trialvolume::GraphCallDescription>();
    this->MakeSlotAvailable(&this->in_graph_data_slot_);

    this->line_width_slot_.SetParameter(new core::param::FloatParam(1.0f, 0.01f, 1000.0f));
    this->MakeSlotAvailable(&this->line_width_slot_);

    last_data_hash_ = 0;
    vbo = 0;
    va = 0;
    ibo = 0;
}

/*
 * TrackingGraphRenderer::~TrackingGraphRenderer
 */
TrackingGraphRenderer::~TrackingGraphRenderer() {
    this->Release();
}

/*
 * TrackingGraphRenderer::create
 */
bool TrackingGraphRenderer::create() {

    // TUTORIAL Shader creation should always happen in the create method of a renderer.

    using namespace megamol::core::utility::log;

    auto const shader_options = ::msf::ShaderFactoryOptionsOpenGL(this->GetCoreInstance()->GetShaderPaths());

    try {
        line_shader_ = core::utility::make_glowl_shader(
            "trialvol_lines", shader_options, "trialvolume_gl/lines.vert.glsl", "trialvolume_gl/lines.frag.glsl");

    } catch (std::exception& e) {
        Log::DefaultLog.WriteError(("TrackingGraphRenderer: " + std::string(e.what())).c_str());
        return false;
    }

    return true;
}

/*
 * TrackingGraphRenderer::GetExtents
 */
bool TrackingGraphRenderer::GetExtents(mmstd_gl::CallRender3DGL& call) {
    mmstd_gl::CallRender3DGL* cr3d = dynamic_cast<mmstd_gl::CallRender3DGL*>(&call);
    if (cr3d == nullptr)
        return false;

    auto* gc = this->in_graph_data_slot_.CallAs<trialvolume::GraphCall>();
    if (gc == nullptr)
        return false;
    // TODO add the extent call to the graph graph data.
    // if (!(*gc)(GraphCall::CallForGetExtent))

    // cr3d->AccessBoundingBoxes() = gc->AccessBoundingBoxes();
    // cr3d->SetTimeFramesCount(gc->FrameCount());

    return true;
}

/*
 * TrackingGraphRenderer::release
 */
void TrackingGraphRenderer::release() {
    if (va != 0) {
        glDeleteVertexArrays(1, &va);
    }
    if (vbo != 0) {
        glDeleteBuffers(1, &vbo);
    }
    if (ibo != 0) {
        glDeleteBuffers(1, &ibo);
    }
}

/*
 * TrackingGraphRenderer::Render
 */
bool TrackingGraphRenderer::Render(mmstd_gl::CallRender3DGL& call) {
    mmstd_gl::CallRender3DGL* cr3d = dynamic_cast<mmstd_gl::CallRender3DGL*>(&call);
    if (cr3d == nullptr)
        return false;

    // before rendering, call all necessary data
    auto* gc = this->in_graph_data_slot_.CallAs<trialvolume::GraphCall>();
    if (gc == nullptr)
        return false;
    // if (!(*gc)(GraphCall::CallForGetExtent))
    //     return false;
    if (!(*gc)(0)) // GraphCall::CallForGetData
        return false;

    auto const graph = gc->GetGraph();
    auto const cluster_graph = dynamic_cast<trialvolume::ClusterGraph*>(graph.get());
    if (cluster_graph == nullptr) {
        return false;
    }

    // only reload the vertex array if the data has changed
    if (gc->DataHash() != last_data_hash_) {
        last_data_hash_ = gc->DataHash();

        if (va == 0 || vbo == 0 || ibo == 0) { // generate new buffers only if they do not exist
            glGenVertexArrays(1, &va);
            glGenBuffers(1, &vbo);
            glGenBuffers(1, &ibo);
        }

        // get the data
        const auto& clusters = cluster_graph->nodes_;
        const auto& cluster_connections = cluster_graph->edges_;

        if (clusters.size() == 0) {
            return false;
        }
        if (cluster_connections.size() == 0) {
            return false;
        }

        // Copy the data into a contiguous array
        struct cluster_data_t {
            glm::vec3 bbox_min;
            GLuint frame;
            glm::vec3 bbox_max;
            GLuint frame_local_id;
            glm::vec3 center_of_mass;
            GLfloat total_mass;
            glm::vec3 velocity;
            GLfloat padding_15;
        };

        std::vector<cluster_data_t> cluster_data;
        std::map<uint64_t, size_t> id_to_index;
        cluster_data.reserve(clusters.size());

        glm::vec3 min_pos = glm::vec3(std::numeric_limits<float>::max());
        glm::vec3 max_pos = glm::vec3(std::numeric_limits<float>::min());

        for (auto const& cluster : clusters) {
            auto const& cluster_info = dynamic_cast<trialvolume::ClusterInfo*>(cluster.get());
            if (cluster_info == nullptr) {
                continue;
            }

            id_to_index[cluster_info->id] = cluster_data.size();

            cluster_data_t data;
            data.bbox_min.x = cluster_info->bounding_box.Left();
            data.bbox_min.y = cluster_info->bounding_box.Bottom();
            data.bbox_min.z = cluster_info->bounding_box.Back();
            data.frame = cluster_info->frame;
            data.bbox_max.x = cluster_info->bounding_box.Right();
            data.bbox_max.y = cluster_info->bounding_box.Top();
            data.bbox_max.z = cluster_info->bounding_box.Front();
            data.frame_local_id = cluster_info->frame_local_id;
            data.center_of_mass.x = cluster_info->center_of_mass.X();
            data.center_of_mass.y = cluster_info->center_of_mass.Y();
            data.center_of_mass.z = cluster_info->center_of_mass.Z();
            data.total_mass = cluster_info->total_mass;
            data.velocity.x = cluster_info->velocity.X();
            data.velocity.y = cluster_info->velocity.Y();
            data.velocity.z = cluster_info->velocity.Z();
            data.padding_15 = 0.0f;
            cluster_data.push_back(data);

            min_pos = glm::min(min_pos, data.bbox_min);
            max_pos = glm::max(max_pos, data.bbox_max);
        }

        bbox_ = vislib::math::Cuboid<float>(min_pos.x, min_pos.y, min_pos.z, max_pos.x, max_pos.y, max_pos.z);

        // load the data into the vertex buffer
        glBindVertexArray(va);
        glBindBuffer(GL_ARRAY_BUFFER, vbo);
        glEnableVertexAttribArray(0); // bbox_min
        glEnableVertexAttribArray(1); // bbox_max
        glEnableVertexAttribArray(2); // center_of_mass
        glEnableVertexAttribArray(3); // velocity
        glEnableVertexAttribArray(4); // total_mass
        glEnableVertexAttribArray(5); // frame
        glEnableVertexAttribArray(6); // frame_local_id

        glBufferData(GL_ARRAY_BUFFER, sizeof(cluster_data_t) * cluster_data.size(), cluster_data.data(), GL_STATIC_DRAW);

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), (void*)offsetof(cluster_data_t, bbox_min));
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), (void*)offsetof(cluster_data_t, bbox_max));
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), (void*)offsetof(cluster_data_t, center_of_mass));
        glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), (void*)offsetof(cluster_data_t, velocity));
        glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), (void*)offsetof(cluster_data_t, total_mass));
        glVertexAttribIPointer(5, 1, GL_UNSIGNED_INT, sizeof(cluster_data_t), (void*)offsetof(cluster_data_t, frame));
        glVertexAttribIPointer(6, 1, GL_UNSIGNED_INT, sizeof(cluster_data_t), (void*)offsetof(cluster_data_t, frame_local_id));

        // Generate the index buffer
        std::vector<GLuint> indices;
        indices.reserve(cluster_connections.size() * 2);
        for (auto const& connection : cluster_connections) {
            auto const& cluster_info_1 = dynamic_cast<trialvolume::ClusterInfo*>(connection->source.get());
            auto const& cluster_info_2 = dynamic_cast<trialvolume::ClusterInfo*>(connection->target.get());
            if (cluster_info_1 == nullptr || cluster_info_2 == nullptr) {
                continue;
            }

            indices.push_back(id_to_index[cluster_info_1->id]);
            indices.push_back(id_to_index[cluster_info_2->id]);
        }
        num_indices_ = indices.size();

        // load the data into the index buffer
        glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, ibo);
        glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(GLuint) * indices.size(), indices.data(), GL_STATIC_DRAW);

        // unbind the buffers
        glBindVertexArray(0);
    }

    cr3d->AccessBoundingBoxes().SetBoundingBox(bbox_);

    core::view::Camera cam = cr3d->GetCamera();

    auto view = cam.getViewMatrix();
    auto proj = cam.getProjectionMatrix();
    auto mvp = proj * view;
    auto cam_pose = cam.get<core::view::Camera::Pose>();

    // start the rendering

    // Scale the point size with the parameter
    glPointSize(this->line_width_slot_.Param<core::param::FloatParam>()->Value());

    // Use the line shader
    this->line_shader_->use();

    // set all uniforms for the shaders
    this->line_shader_->setUniform("mvp", mvp);
    this->line_shader_->setUniform("view", view);
    this->line_shader_->setUniform("proj", proj);
    this->line_shader_->setUniform("camRight", cam_pose.right.x, cam_pose.right.y, cam_pose.right.z);
    this->line_shader_->setUniform("camUp", cam_pose.up.x, cam_pose.up.y, cam_pose.up.z);
    this->line_shader_->setUniform("camPos", cam_pose.position.x, cam_pose.position.y, cam_pose.position.z);
    this->line_shader_->setUniform("camDir", cam_pose.direction.x, cam_pose.direction.y, cam_pose.direction.z);
    this->line_shader_->setUniform("scalingFactor", this->line_width_slot_.Param<core::param::FloatParam>()->Value());

    // draw the lines
    glBindVertexArray(va);
    glEnable(GL_DEPTH_TEST);

    glDrawElements(GL_LINES, num_indices_, GL_UNSIGNED_INT, nullptr);

    glBindVertexArray(0);
    glDisable(GL_DEPTH_TEST);

    glUseProgram(0);

    return true;
}
