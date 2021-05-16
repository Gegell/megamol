#include "stdafx.h"
#include "TriangleMeshRenderer3D.h"

#include "mesh/MeshDataCall.h"
#include "mesh/TriangleMeshCall.h"

#include "mesh/GPUMeshCollection.h"
#include "mesh/MeshCalls.h"

#include "mmcore/BoundingBoxes_2.h"
#include "mmcore/Call.h"
#include "mmcore/param/ColorParam.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/FlexEnumParam.h"
#include "mmcore/param/TransferFunctionParam.h"
#include "mmcore/utility/DataHash.h"
#include "mmcore/view/CallClipPlane.h"

#include "glad/glad.h"

#include "glowl/VertexLayout.hpp"

#include <memory>
#include <sstream>
#include <string>
#include <utility>
#include <vector>

namespace megamol {
namespace mesh {

    TriangleMeshRenderer3D::TriangleMeshRenderer3D()
            : triangle_mesh_slot("get_triangle_mesh", "Triangle mesh input")
            , mesh_data_slot("get_mesh_data", "Mesh data input")
            , clip_plane_slot("clip_plane", "Clip plane for clipping the rendered triangle mesh")
            , data_set("data_set", "Data set used for coloring the triangles")
            , default_color("default_color", "Default color if no dataset is specified")
            , btf_filename_slot("BTF filename", "Shader path")
            , triangle_mesh_hash(-681861)
            , triangle_mesh_changed(false)
            , mesh_data_hash(-3618424)
            , mesh_data_changed(false)
            , rhs_gpu_tasks_version(0)
            , material_collection(nullptr)
            , btf_file_changed(false) {

        // Connect input slots
        this->triangle_mesh_slot.SetCompatibleCall<TriangleMeshCall::triangle_mesh_description>();
        this->MakeSlotAvailable(&this->triangle_mesh_slot);

        this->mesh_data_slot.SetCompatibleCall<MeshDataCall::mesh_data_description>();
        this->MakeSlotAvailable(&this->mesh_data_slot);

        this->clip_plane_slot.SetCompatibleCall<core::view::CallClipPlaneDescription>();
        this->MakeSlotAvailable(&this->clip_plane_slot);

        // Connect parameter slots
        this->data_set << new core::param::FlexEnumParam("");
        this->MakeSlotAvailable(&this->data_set);

        this->default_color << new core::param::ColorParam(0.7f, 0.7f, 0.7f, 1.0f);
        this->MakeSlotAvailable(&this->default_color);

        this->btf_filename_slot << new core::param::FilePathParam("triangle_mesh");
        this->MakeSlotAvailable(&this->btf_filename_slot);

        // Disconnect inherited slots
        this->SetSlotUnavailable(&this->m_mesh_slot);
    }

    TriangleMeshRenderer3D::~TriangleMeshRenderer3D() {
        this->Release();
    }

    bool TriangleMeshRenderer3D::create() {
        mesh::AbstractGPURenderTaskDataSource::create();

        this->material_collection = std::make_shared<GPUMaterialCollection>();
        this->material_collection->addMaterial(this->instance(), "triangle_mesh", "triangle_mesh");

        return true;
    }

    void TriangleMeshRenderer3D::release() {}

    bool TriangleMeshRenderer3D::get_input_data() {
        auto tmc_ptr = this->triangle_mesh_slot.CallAs<TriangleMeshCall>();
        auto mdc_ptr = this->mesh_data_slot.CallAs<MeshDataCall>();

        if (tmc_ptr == nullptr) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "Triangle mesh input is not connected. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);

            return false;
        }

        auto& tmc = *tmc_ptr;

        if (!tmc(0)) {
            if (tmc.DataHash() != this->triangle_mesh_hash) {
                megamol::core::utility::log::Log::DefaultLog.WriteError(
                    "Error getting triangle mesh. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);

                this->triangle_mesh_hash = tmc.DataHash();
            }

            return false;
        }

        if (tmc.DataHash() != this->triangle_mesh_hash) {
            this->render_data.vertices = tmc.get_vertices();
            this->render_data.normals = tmc.get_normals();
            this->render_data.indices = tmc.get_indices();

            this->triangle_mesh_hash = tmc.DataHash();
            this->triangle_mesh_changed = true;
        }

        if (this->btf_filename_slot.IsDirty()) {
            this->btf_filename_slot.ResetDirty();

            auto vislib_filename = this->btf_filename_slot.Param<core::param::FilePathParam>()->Value();
            std::string filename(vislib_filename.PeekBuffer());

            this->material_collection->clear();
            this->material_collection->addMaterial(this->instance(), filename, filename);

            this->btf_file_changed = true;
        }

        if (mdc_ptr != nullptr && (*mdc_ptr)(0) && !mdc_ptr->get_data_sets().empty() &&
            mdc_ptr->DataHash() != this->mesh_data_hash) {

            this->data_set.Param<core::param::FlexEnumParam>()->ClearValues();

            for (const auto& data_set_name : mdc_ptr->get_data_sets()) {
                this->data_set.Param<core::param::FlexEnumParam>()->AddValue(data_set_name);
            }

            this->mesh_data_hash = mdc_ptr->DataHash();
        }

        if (mdc_ptr != nullptr) {
            auto data = mdc_ptr->get_data(this->data_set.Param<core::param::FlexEnumParam>()->Value());

            if (data != nullptr && this->data_set.IsDirty()) {
                this->render_data.values = data;

                this->data_set.ResetDirty();
                this->mesh_data_changed = true;
            }
        }

        if (this->render_data.values == nullptr || this->default_color.IsDirty()) {
            this->render_data.values = std::make_shared<MeshDataCall::data_set>();

            this->render_data.values->min_value = 0.0f;
            this->render_data.values->max_value = 1.0f;

            const auto color = this->default_color.Param<core::param::ColorParam>()->Value();
            this->default_color.ResetDirty();

            std::stringstream ss;
            ss << "{\"Interpolation\":\"LINEAR\",\"Nodes\":["
               << "[" << color[0] << "," << color[1] << "," << color[2] << "," << color[3]
               << ",0.0,0.05000000074505806],"
               << "[" << color[0] << "," << color[1] << "," << color[2] << "," << color[3]
               << ",1.0,0.05000000074505806]]"
               << ",\"TextureSize\":2,\"ValueRange\":[0.0,1.0]}";

            this->render_data.values->transfer_function = ss.str();
            this->render_data.values->transfer_function_dirty = true;

            this->render_data.values->data =
                std::make_shared<std::vector<GLfloat>>(this->render_data.vertices->size() / 3, 1.0f);

            this->mesh_data_changed = true;
        }

        return true;
    }

    bool TriangleMeshRenderer3D::get_input_extent() {
        auto tmc_ptr = this->triangle_mesh_slot.CallAs<TriangleMeshCall>();
        auto mdc_ptr = this->triangle_mesh_slot.CallAs<MeshDataCall>();

        if (tmc_ptr == nullptr) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "Triangle mesh input is not connected. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);

            return false;
        }

        auto& tmc = *tmc_ptr;

        if (!tmc(1)) {
            if (tmc.DataHash() != this->triangle_mesh_hash) {
                megamol::core::utility::log::Log::DefaultLog.WriteError(
                    "Error getting extents for the triangle mesh. [%s, %s, line %d]\n", __FILE__, __FUNCTION__,
                    __LINE__);

                this->triangle_mesh_hash = tmc.DataHash();
            }

            return false;
        }

        if (tmc.get_dimension() != TriangleMeshCall::dimension_t::THREE) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "Input triangle mesh must be three-dimensional. [%s, %s, line %d]\n", __FILE__, __FUNCTION__, __LINE__);

            return false;
        }

        this->bounding_box = tmc.get_bounding_box();

        return true;
    }

    bool TriangleMeshRenderer3D::getDataCallback(core::Call& call) {
        auto& grtc = static_cast<mesh::CallGPURenderTaskData&>(call);

        if (!get_input_data()) {
            return false;
        }

        // Get render task from rhs
        CallGPURenderTaskData* rhs_rtc = this->m_renderTask_rhs_slot.CallAs<CallGPURenderTaskData>();

        std::vector<std::shared_ptr<GPURenderTaskCollection>> rt_collections;

        if (rhs_rtc != nullptr) {
            if (!(*rhs_rtc)(0)) {
                return false;
            }

            if (rhs_rtc->hasUpdate()) {
                ++this->rhs_gpu_tasks_version;
            }

            rt_collections = rhs_rtc->getData();
        }

        rt_collections.push_back(this->m_rendertask_collection.first);

        const std::string identifier("triangle_mesh");

        if (this->triangle_mesh_changed || this->mesh_data_changed || this->btf_file_changed) {
            clearRenderTaskCollection();

            // Create mesh
            this->render_data.mesh = std::make_shared<mesh::GPUMeshCollection>();

            using vbi_t = typename std::vector<GLfloat>::iterator;
            using ibi_t = typename std::vector<GLuint>::iterator;

            std::vector<glowl::VertexLayout> vertex_descriptors{
                glowl::VertexLayout(3 * sizeof(float), {glowl::VertexLayout::Attribute(3, GL_FLOAT, GL_FALSE, 0)}),
                glowl::VertexLayout(1 * sizeof(float), {glowl::VertexLayout::Attribute(1, GL_FLOAT, GL_FALSE, 0)})};

            if (this->render_data.normals != nullptr) {
                vertex_descriptors.push_back(
                    glowl::VertexLayout(3 * sizeof(float), {glowl::VertexLayout::Attribute(3, GL_FLOAT, GL_TRUE, 0)}));
            }

            std::vector<std::pair<vbi_t, vbi_t>> vertex_buffer{
                {this->render_data.vertices->begin(), this->render_data.vertices->end()},
                {this->render_data.values->data->begin(), this->render_data.values->data->end()}};

            if (this->render_data.normals != nullptr) {
                vertex_buffer.push_back({this->render_data.normals->begin(), this->render_data.normals->end()});
            }

            std::pair<ibi_t, ibi_t> index_buffer{this->render_data.indices->begin(), this->render_data.indices->end()};

            this->render_data.mesh->template addMesh<vbi_t, ibi_t>(identifier, vertex_descriptors, vertex_buffer,
                index_buffer, GL_UNSIGNED_INT, GL_STATIC_DRAW, GL_TRIANGLES);

            // Create render task
            const auto& mesh_data = this->render_data.mesh->getSubMeshData().at(identifier);

            std::memcpy(&this->render_data.per_draw_data[per_draw_data_t::offset_min_value],
                &this->render_data.values->min_value, per_draw_data_t::size_min_value);
            std::memcpy(&this->render_data.per_draw_data[per_draw_data_t::offset_max_value],
                &this->render_data.values->max_value, per_draw_data_t::size_max_value);

            this->m_rendertask_collection.first->clear();
            this->m_rendertask_collection.first->addRenderTask(identifier,
                this->material_collection->getMaterials().begin()->second.shader_program, mesh_data.mesh->mesh,
                mesh_data.sub_mesh_draw_command, this->render_data.per_draw_data);
        }

        if (this->render_data.values->transfer_function_dirty || this->triangle_mesh_changed ||
            this->mesh_data_changed || this->btf_file_changed) {
            // Create texture for transfer function
            std::vector<GLfloat> texture_data;
            int transfer_function_size, _unused__height;

            const auto valid_tf = core::param::TransferFunctionParam::GetTextureData(
                this->render_data.values->transfer_function, texture_data, transfer_function_size, _unused__height);

            if (!valid_tf) {
                return false;
            }

            if (this->render_data.transfer_function != 0) {
                glDeleteTextures(1, &this->render_data.transfer_function);
            }

            glGenTextures(1, &this->render_data.transfer_function);

            glActiveTexture(GL_TEXTURE0);
            glBindTexture(GL_TEXTURE_1D, this->render_data.transfer_function);

            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
            glTexParameteri(GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);

            glTexImage1D(GL_TEXTURE_1D, 0, GL_RGBA, static_cast<GLsizei>(transfer_function_size), 0, GL_RGBA, GL_FLOAT,
                static_cast<GLvoid*>(texture_data.data()));

            glBindTexture(GL_TEXTURE_1D, 0);

            const auto transfer_function_handle = glGetTextureHandleARB(this->render_data.transfer_function);
            glMakeTextureHandleResidentARB(transfer_function_handle);

            // Update per draw data
            std::memcpy(&this->render_data.per_draw_data[per_draw_data_t::offset_tf], &transfer_function_handle,
                per_draw_data_t::size_tf);

            this->m_rendertask_collection.first->updatePerDrawData(identifier, this->render_data.per_draw_data);
        }

        {
            auto cp = this->clip_plane_slot.CallAs<core::view::CallClipPlane>();

            if (cp != nullptr && (*cp)(0)) {
                // Set clip plane flag to enabled
                const int use_plane = 1;

                std::memcpy(&this->render_data.per_draw_data[per_draw_data_t::offset_plane_bool], &use_plane,
                    per_draw_data_t::size_plane_bool);

                // Get clip plane
                const auto plane = cp->GetPlane();
                const std::array<float, 4> abcd_plane{plane.A(), plane.B(), plane.C(), plane.D()};

                std::memcpy(&this->render_data.per_draw_data[per_draw_data_t::offset_plane], abcd_plane.data(),
                    per_draw_data_t::size_plane);
            } else {
                // Set clip plane flag to disabled
                const int use_plane = 0;

                std::memcpy(&this->render_data.per_draw_data[per_draw_data_t::offset_plane_bool], &use_plane,
                    per_draw_data_t::size_plane_bool);
            }

            this->m_rendertask_collection.first->updatePerDrawData(identifier, this->render_data.per_draw_data);
        }

        this->triangle_mesh_changed = false;
        this->mesh_data_changed = false;
        this->render_data.values->transfer_function_dirty = false;
        this->btf_file_changed = false;

        grtc.setData(rt_collections,
            core::utility::DataHash(this->triangle_mesh_hash, this->mesh_data_hash, this->rhs_gpu_tasks_version));

        return true;
    }

    bool TriangleMeshRenderer3D::getMetaDataCallback(core::Call& call) {
        auto& grtc = static_cast<mesh::CallGPURenderTaskData&>(call);

        if (!get_input_extent()) {
            return false;
        }

        core::BoundingBoxes_2 bbox;
        bbox.SetBoundingBox(this->bounding_box);
        bbox.SetClipBox(this->bounding_box);

        grtc.setMetaData(core::Spatial3DMetaData{1, 0, bbox});

        return true;
    }

} // namespace mesh
} // namespace megamol
