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

using namespace megamol;
using namespace megamol::trialvolume_gl;

/*
 * TrackingGraphRenderer::TrackingGraphRenderer
 */
TrackingGraphRenderer::TrackingGraphRenderer()
        : mmstd_gl::Renderer3DModuleGL()
        , in_graph_data_slot_("inData", "The input data slot for the graph data.")
        , line_width_slot_("line width", "Width of the connecting lines") {
    this->in_graph_data_slot_.SetCompatibleCall<GraphCall>();
    this->MakeSlotAvailable(&this->in_graph_data_slot_);

    this->line_width_slot_.SetParameter(new core::param::FloatParam(1.0f, 0.01f, 1000.0f));
    this->MakeSlotAvailable(&this->line_width_slot_);

    last_data_hash_ = 0;
    vbo = 0;
    va = 0;
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

    GraphCall* gc = this->in_graph_data_slot_.CallAs<GraphCall>();
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
}

/*
 * TrackingGraphRenderer::Render
 */
bool TrackingGraphRenderer::Render(mmstd_gl::CallRender3DGL& call) {
    mmstd_gl::CallRender3DGL* cr3d = dynamic_cast<mmstd_gl::CallRender3DGL*>(&call);
    if (cr3d == nullptr)
        return false;

    // before rendering, call all necessary data
    GraphCall* gc = this->in_graph_data_slot_.CallAs<GraphCall>();
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

        if (va == 0 || vbo == 0) { // generate new buffers only if they do not exist
            glGenVertexArrays(1, &va);
            glGenBuffers(1, &vbo);
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
        cluster_data.reserve(clusters.size());
        for (auto const& cluster : clusters) {
            cluster_data_t data;
            data.bbox_min.x = cluster->bounding_box.getLeft();
            data.bbox_min.y = cluster->bounding_box.getBottom();
            data.bbox_min.z = cluster->bounding_box.getBack();
            data.frame = cluster->frame;
            data.bbox_max.x = cluster->bounding_box.getRight();
            data.bbox_max.y = cluster->bounding_box.getTop();
            data.bbox_max.z = cluster->bounding_box.getFront();
            data.frame_local_id = cluster->frame_local_id;
            data.center_of_mass.x = cluster->center_of_mass.X();
            data.center_of_mass.x = cluster->center_of_mass.Y();
            data.center_of_mass.x = cluster->center_of_mass.Z();
            data.total_mass = cluster->total_mass;
            data.velocity.x = cluster->velocity.X();
            data.velocity.y = cluster->velocity.Y();
            data.velocity.z = cluster->velocity.Z();
            data.padding_15 = 0.0f;
            cluster_data.push_back(data);
        }

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

        glBufferData(GL_ARRAY_BUFFER, sizeof(cluster_data_t) * clusters.size(), cluster_data.data());

        glVertexAttribPointer(0, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), 0);
        glVertexAttribPointer(1, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), 4);
        glVertexAttribPointer(2, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), 8);
        glVertexAttribPointer(3, 3, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), 12);
        glVertexAttribPointer(4, 1, GL_FLOAT, GL_FALSE, sizeof(cluster_data_t), 11);
        glVertexAttribPointer(5, 1, GL_UINT, GL_FALSE, sizeof(cluster_data_t), 3);
        glVertexAttribPointer(6, 1, GL_UINT, GL_FALSE, sizeof(cluster_data_t), 7);

        glBindBuffer(GL_ARRAY_BUFFER, 0);
        glBindVertexArray(0);
    }

    // TODO add bbox to graph data?
    // cr3d->AccessBoundingBoxes() = gc->AccessBoundingBoxes();

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

    glBindVertexArray(va);
    glEnableVertexAttribArray(0);
    glEnableVertexAttribArray(1);
    glEnableVertexAttribArray(2);
    glEnableVertexAttribArray(3);
    glEnableVertexAttribArray(4);
    glEnableVertexAttribArray(5);
    glEnableVertexAttribArray(6);

    // set all uniforms for the shaders
    this->line_shader_->setUniform("mvp", mvp);
    this->line_shader_->setUniform("view", view);
    this->line_shader_->setUniform("proj", proj);
    this->line_shader_->setUniform("camRight", cam_pose.right.x, cam_pose.right.y, cam_pose.right.z);
    this->line_shader_->setUniform("camUp", cam_pose.up.x, cam_pose.up.y, cam_pose.up.z);
    this->line_shader_->setUniform("camPos", cam_pose.position.x, cam_pose.position.y, cam_pose.position.z);
    this->line_shader_->setUniform("camDir", cam_pose.direction.x, cam_pose.direction.y, cam_pose.direction.z);
    this->line_shader_->setUniform("scalingFactor", this->line_width_slot_.Param<core::param::FloatParam>()->Value());

    glEnable(GL_DEPTH_TEST);

    // draw one point for each sphere
    glDrawArrays(GL_LINES, 0, static_cast<GLsizei>(cluster_connections.size()));
    glBindVertexArray(0);

    glDisable(GL_DEPTH_TEST);

    glUseProgram(0);

    return true;
}
