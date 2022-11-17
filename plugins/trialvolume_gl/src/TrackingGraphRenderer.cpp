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
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/ColorParam.h"
#include "mmcore_gl/utility/ShaderFactory.h"
#include "mmstd_gl/renderer/CallRender3DGL.h"
#include "mmstd_gl/renderer/CallGetTransferFunctionGL.h"
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
        , in_tf_slot_("inTransferfunction", "The slot for the transfer function module")
        , line_width_slot_("line width", "Width of the connecting lines")
        , draw_connections_slot_("drawConnections", "Draw the connections between the nodes")
        , draw_bboxes_slot_("drawBBoxes", "Draw the bounding boxes of the nodes")
        , bbox_color_slot_("bboxColor", "Color of the bounding boxes")
        , filter_min_mass_slot_("filter::min_mass", "Filter the nodes by mass")
        , color_mode_slot_("colorMode", "The color mode for the nodes") {
    this->in_graph_data_slot_.SetCompatibleCall<trialvolume::GraphCallDescription>();
    this->MakeSlotAvailable(&this->in_graph_data_slot_);

    this->in_tf_slot_.SetCompatibleCall<mmstd_gl::CallGetTransferFunctionGLDescription>();
    this->MakeSlotAvailable(&this->in_tf_slot_);

    this->line_width_slot_.SetParameter(new core::param::FloatParam(1.0f, 0.01f, 1000.0f));
    this->MakeSlotAvailable(&this->line_width_slot_);

    this->draw_connections_slot_.SetParameter(new core::param::BoolParam(true));
    this->MakeSlotAvailable(&this->draw_connections_slot_);

    this->draw_bboxes_slot_.SetParameter(new core::param::BoolParam(true));
    this->MakeSlotAvailable(&this->draw_bboxes_slot_);

    this->bbox_color_slot_.SetParameter(new core::param::ColorParam(0.0f, 0.0f, 0.0f, 1.0f));
    this->MakeSlotAvailable(&this->bbox_color_slot_);

    this->filter_min_mass_slot_.SetParameter(new core::param::FloatParam(0.0f, 0.0f));
    this->MakeSlotAvailable(&this->filter_min_mass_slot_);

    auto color_mode_param = new core::param::EnumParam(static_cast<int>(ColorMode::VELOCITY_DIR));
    color_mode_param->SetTypePair(static_cast<int>(ColorMode::VELOCITY_DIR), "Velocity direction");
    color_mode_param->SetTypePair(static_cast<int>(ColorMode::VELOCITY_MAG), "Velocity magnitude");
    color_mode_param->SetTypePair(static_cast<int>(ColorMode::TOTAL_MASS), "Total Mass");
    color_mode_param->SetTypePair(static_cast<int>(ColorMode::LOCAL_ID), "Local ID");
    color_mode_param->SetTypePair(static_cast<int>(ColorMode::FRAME), "Frame");
    this->color_mode_slot_.SetParameter(color_mode_param);
    this->MakeSlotAvailable(&this->color_mode_slot_);

    last_data_hash_ = 0;
    vbo = 0;
    va = 0;
    ibo = 0;
    grey_tf_ = 0;
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
            "trialvol_lines", shader_options, "trialvolume_gl/lines.vert.glsl", "trialvolume_gl/lines.geom.glsl", "trialvolume_gl/lines.frag.glsl");

        bbox_shader_ = core::utility::make_glowl_shader(
            "trialvol_bbox", shader_options, "trialvolume_gl/bbox.geom.glsl", "trialvolume_gl/bbox.vert.glsl", "trialvolume_gl/bbox.frag.glsl");
    } catch (std::exception& e) {
        Log::DefaultLog.WriteError(("TrackingGraphRenderer: " + std::string(e.what())).c_str());
        return false;
    }

    // Fallback transfer function texture
    glGenTextures(1, &this->grey_tf_);
    unsigned char tex[6] = {0, 0, 0, 255, 255, 255};
    glBindTexture(GL_TEXTURE_1D, this->grey_tf_);
    glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, 2, 0, GL_RGB, GL_UNSIGNED_BYTE, tex);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_NEAREST);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP);
    glBindTexture(GL_TEXTURE_1D, 0);

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

    cr3d->AccessBoundingBoxes().SetBoundingBox(bbox_);
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
    if (grey_tf_ != 0) {
        glDeleteTextures(1, &grey_tf_);
    }
}

// Joinked from moldyn::SphereRenderer
bool TrackingGraphRenderer::enableTransferFunctionTexture(glowl::GLSLProgram& prgm) {
    mmstd_gl::CallGetTransferFunctionGL* cgtf = in_tf_slot_.CallAs<mmstd_gl::CallGetTransferFunctionGL>();
    if ((cgtf != nullptr) && (*cgtf)(0)) {
        cgtf->BindConvenience(prgm, GL_TEXTURE0, 0);
    } else {
        glEnable(GL_TEXTURE_1D);
        glActiveTexture(GL_TEXTURE0);
        glBindTexture(GL_TEXTURE_1D, grey_tf_);
        glUniform1i(prgm.getUniformLocation("tfTexture"), 0);
        glUniform2fv(prgm.getUniformLocation("tfRange"), 1, static_cast<GLfloat*>(tf_range_.data()));
    }
    return true;
}

// Joinked from moldyn::SphereRenderer
bool TrackingGraphRenderer::disableTransferFunctionTexture(void) {
    mmstd_gl::CallGetTransferFunctionGL* cgtf = in_tf_slot_.CallAs<mmstd_gl::CallGetTransferFunctionGL>();
    if (cgtf != nullptr) {
        cgtf->UnbindConvenience();
    } else {
        glBindTexture(GL_TEXTURE_1D, 0);
        glDisable(GL_TEXTURE_1D);
    }
    return true;
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
        last_frame_ = 0;
        max_local_id_ = 0;
        max_total_mass_ = 0.0f;
        max_velocity_ = 0.0f;

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
            last_frame_ = std::max(last_frame_, static_cast<int>(data.frame));
            max_local_id_ = std::max(max_local_id_, static_cast<int>(data.frame_local_id));
            max_total_mass_ = std::max(max_total_mass_, data.total_mass);
            max_velocity_ = std::max(max_velocity_, glm::length(data.velocity));
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
    cr3d->SetTimeFramesCount(last_frame_ + 1);

    core::view::Camera cam = cr3d->GetCamera();

    auto view = cam.getViewMatrix();
    auto proj = cam.getProjectionMatrix();
    auto mvp = proj * view;

    bool draw_bboxes = draw_bboxes_slot_.Param<core::param::BoolParam>()->Value();
    bool draw_connections = draw_connections_slot_.Param<core::param::BoolParam>()->Value();

    bool draw_anything = draw_bboxes || draw_connections;

    if (!draw_anything) {
        return true;
    }

    float const min_mass = filter_min_mass_slot_.Param<core::param::FloatParam>()->Value();
    auto const color_mode = color_mode_slot_.Param<core::param::EnumParam>()->Value();
    auto const bbox_color = bbox_color_slot_.Param<core::param::ColorParam>()->Value();

    // start the rendering
    glBindVertexArray(va);
    glEnable(GL_DEPTH_TEST);
    glEnable(GL_LINE_SMOOTH);
    glEnable(GL_BLEND);
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);

    // Scale the point size with the parameter
    glLineWidth(this->line_width_slot_.Param<core::param::FloatParam>()->Value());

    if (draw_bboxes) {
        // use the bbox shader
        bbox_shader_->use();
        bbox_shader_->setUniform("mvp", mvp);
        bbox_shader_->setUniform("time", cr3d->Time());
        bbox_shader_->setUniform("min_mass", min_mass);
        bbox_shader_->setUniform("line_color", glm::vec4(bbox_color[0], bbox_color[1], bbox_color[2], bbox_color[3]));

        glDrawElements(GL_LINES, num_indices_, GL_UNSIGNED_INT, nullptr);
    }

    if (draw_connections) {
        // Use the line shader
        line_shader_->use();
        // set all uniforms for the shaders
        line_shader_->setUniform("mvp", mvp);
        line_shader_->setUniform("min_mass", min_mass);
        line_shader_->setUniform("color_mode", color_mode);
        line_shader_->setUniform("max_mass", max_total_mass_);
        line_shader_->setUniform("max_frame", static_cast<float>(last_frame_));
        line_shader_->setUniform("max_frame_local_id", static_cast<float>(max_local_id_));
        line_shader_->setUniform("max_velocity", max_velocity_);
        line_shader_->setUniform("time", cr3d->Time());

        enableTransferFunctionTexture(*line_shader_);
        glDrawElements(GL_LINES, num_indices_, GL_UNSIGNED_INT, nullptr);
        disableTransferFunctionTexture();
    }

    glBindVertexArray(0);
    glDisable(GL_DEPTH_TEST);
    glDisable(GL_LINE_SMOOTH);
    glDisable(GL_BLEND);

    glUseProgram(0);

    return true;
}
