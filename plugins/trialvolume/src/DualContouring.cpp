#include "DualContouring.h"

#include <vector>

#include "geometry_calls/VolumetricDataCall.h"
#include "mesh/TriangleMeshCall.h"
#include "mmcore/param/FloatParam.h"

using namespace megamol;

bool trialvolume::DualContouring::create(void) {
    return true;
}

void trialvolume::DualContouring::release(void) {
    // TODO: release any data here
}

trialvolume::DualContouring::DualContouring()
        : iso_level_slot_("isoLevel", "The isovalue to use for the dual contouring")
        , out_triangle_surface_slot_("outTriangleSurface", "The triangle surface data")
        , in_volume_data_slot_("inVolumeData", "The volume data") {
    // Setup the input volume data slot
    in_volume_data_slot_.SetCompatibleCall<geocalls::VolumetricDataCallDescription>();
    MakeSlotAvailable(&in_volume_data_slot_);

    // Setup the output triangle surface data slot
    out_triangle_surface_slot_.SetCallback(mesh::TriangleMeshCall::ClassName(),
        mesh::TriangleMeshCall::FunctionName(0u), &DualContouring::getTriangleSurfaceCallback);
    out_triangle_surface_slot_.SetCallback(mesh::TriangleMeshCall::ClassName(),
        mesh::TriangleMeshCall::FunctionName(1u), &DualContouring::getExtentCallback);
    MakeSlotAvailable(&out_triangle_surface_slot_);

    // Setup the iso level slot
    iso_level_slot_ << new core::param::FloatParam(0.5f, 0.0f, 1.0f, 0.05f);
    MakeSlotAvailable(&iso_level_slot_);

    // Setup the generated data containers
    index_buffer_ = std::make_shared<std::vector<unsigned int>>();
    vertex_buffer_ = std::make_shared<std::vector<float>>();
    normal_buffer_ = std::make_shared<std::vector<float>>();
}

trialvolume::DualContouring::~DualContouring() {
    release();
}

bool trialvolume::DualContouring::getTriangleSurfaceCallback(core::Call& caller) {
    auto* volumeDataCall = in_volume_data_slot_.CallAs<geocalls::VolumetricDataCall>();
    if (volumeDataCall == nullptr) {
        return false;
    }

    // Check if the data has changed
    if (volumeDataCall->DataHash() != in_data_hash_ || anythingDirty()) {
        in_data_hash_ = volumeDataCall->DataHash();
        resetDirtyFlags();
        computeSurface(*volumeDataCall);
        data_hash_++;

        core::utility::log::Log::DefaultLog.WriteInfo("[DualContouring] Triangles: %d", index_buffer_->size() / 3);
    }

    auto* triangleMeshCall = dynamic_cast<mesh::TriangleMeshCall*>(&caller);

    triangleMeshCall->set_vertices(vertex_buffer_);
    triangleMeshCall->set_indices(index_buffer_);
    triangleMeshCall->set_normals(normal_buffer_);
    triangleMeshCall->set_bounding_box(bbox_);
    triangleMeshCall->set_dimension(mesh::TriangleMeshCall::dimension_t::THREE);
    triangleMeshCall->SetDataHash(data_hash_);

    return true;
}

bool trialvolume::DualContouring::computeSurface(geocalls::VolumetricDataCall& volumeDataCall) {

    volumeDataCall(geocalls::VolumetricDataCall::IDX_GET_DATA);

    auto const metadata = volumeDataCall.GetMetadata();
    if (metadata == nullptr) {
        core::utility::log::Log::DefaultLog.WriteError("[DualContouring] No metadata available");
        return false;
    }

    // Clear previous data
    index_buffer_->clear();
    vertex_buffer_->clear();
    normal_buffer_->clear();

    // Compute the vertices
    for (auto z = 0u; z < metadata->Resolution[2] - 1; ++z) {
        for (auto y = 0u; y < metadata->Resolution[1] - 1; ++y) {
            for (auto x = 0u; x < metadata->Resolution[0] - 1; ++x) {
                // FIXME Assume the coordinates are uniformly spaced
                vertex_buffer_->push_back(
                    (x + 0.5f) * metadata->Extents[0] / metadata->Resolution[0] + metadata->Origin[0]);
                vertex_buffer_->push_back(
                    (y + 0.5f) * metadata->Extents[1] / metadata->Resolution[1] + metadata->Origin[1]);
                vertex_buffer_->push_back(
                    (z + 0.5f) * metadata->Extents[2] / metadata->Resolution[2] + metadata->Origin[2]);

                // FIXME TEMPORARY: Compute the normals via forward differencing
                auto const& voxel = volumeDataCall.GetAbsoluteVoxelValue(x, y, z);
                auto const dx = volumeDataCall.GetAbsoluteVoxelValue(x + 1, y, z) - voxel;
                auto const dy = volumeDataCall.GetAbsoluteVoxelValue(x, y + 1, z) - voxel;
                auto const dz = volumeDataCall.GetAbsoluteVoxelValue(x, y, z + 1) - voxel;
                auto const length = std::sqrt(dx * dx + dy * dy + dz * dz);
                normal_buffer_->push_back(dx / length);
                normal_buffer_->push_back(dy / length);
                normal_buffer_->push_back(dz / length);
            }
        }
    }

    // Set the bounding box
    // bbox = volumeDataCall->AccessBoundingBoxes().ObjectSpaceBBox();
    bbox_ = vislib::math::Cuboid<float>(metadata->Origin[0], metadata->Origin[1], metadata->Origin[2],
        metadata->Origin[0] + metadata->Extents[0], metadata->Origin[1] + metadata->Extents[1],
        metadata->Origin[2] + metadata->Extents[2]);

    // Compute the triangles
    auto isoLevel = iso_level_slot_.Param<core::param::FloatParam>()->Value();
    auto quad_buffer = std::vector<unsigned int>();
    for (auto z = 0u; z < metadata->Resolution[2] - 1; ++z) {
        for (auto y = 0u; y < metadata->Resolution[1] - 1; ++y) {
            for (auto x = 0u; x < metadata->Resolution[0] - 1; ++x) {
                if (x > 0 && y > 0) {
                    auto level1 = volumeDataCall.GetRelativeVoxelValue(x, y, z) - isoLevel;
                    auto level2 = volumeDataCall.GetRelativeVoxelValue(x, y, z + 1) - isoLevel;
                    if (!sameSign(level1, level2)) {
                        // Push a quad as 2 triangles
                        quad_buffer.clear();
                        quad_buffer.push_back(toFlatIndex(x - 1, y - 0, z, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 1, y - 1, z, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 0, y - 0, z, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 0, y - 0, z, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 1, y - 1, z, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 0, y - 1, z, metadata));
                        // Reverse if other is inside
                        if (level1 < 0.0f) {
                            index_buffer_->insert(index_buffer_->end(), quad_buffer.rbegin(), quad_buffer.rend());
                        } else {
                            index_buffer_->insert(index_buffer_->end(), quad_buffer.begin(), quad_buffer.end());
                        }
                    }
                }
                if (x > 0 && z > 0) {
                    auto level1 = volumeDataCall.GetRelativeVoxelValue(x, y, z) - isoLevel;
                    auto level2 = volumeDataCall.GetRelativeVoxelValue(x, y + 1, z) - isoLevel;
                    if (!sameSign(level1, level2)) {
                        // Push a quad as 2 triangles
                        quad_buffer.clear();
                        quad_buffer.push_back(toFlatIndex(x - 1, y, z - 0, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 0, y, z - 0, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 1, y, z - 1, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 0, y, z - 0, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 0, y, z - 1, metadata));
                        quad_buffer.push_back(toFlatIndex(x - 1, y, z - 1, metadata));
                        // Reverse if other is inside
                        if (level1 < 0.0f) {
                            index_buffer_->insert(index_buffer_->end(), quad_buffer.rbegin(), quad_buffer.rend());
                        } else {
                            index_buffer_->insert(index_buffer_->end(), quad_buffer.begin(), quad_buffer.end());
                        }
                    }
                }
                if (y > 0 && z > 0) {
                    auto level1 = volumeDataCall.GetRelativeVoxelValue(x, y, z) - isoLevel;
                    auto level2 = volumeDataCall.GetRelativeVoxelValue(x + 1, y, z) - isoLevel;
                    if (!sameSign(level1, level2)) {
                        // Push a quad as 2 triangles
                        quad_buffer.clear();
                        quad_buffer.push_back(toFlatIndex(x, y - 1, z - 0, metadata));
                        quad_buffer.push_back(toFlatIndex(x, y - 1, z - 1, metadata));
                        quad_buffer.push_back(toFlatIndex(x, y - 0, z - 0, metadata));
                        quad_buffer.push_back(toFlatIndex(x, y - 0, z - 0, metadata));
                        quad_buffer.push_back(toFlatIndex(x, y - 1, z - 1, metadata));
                        quad_buffer.push_back(toFlatIndex(x, y - 0, z - 1, metadata));
                        // Reverse if other is inside
                        if (level1 < 0.0f) {
                            index_buffer_->insert(index_buffer_->end(), quad_buffer.rbegin(), quad_buffer.rend());
                        } else {
                            index_buffer_->insert(index_buffer_->end(), quad_buffer.begin(), quad_buffer.end());
                        }
                    }
                }
            }
        }
    }
    // TODO Compute normals

    // HACK add a single triangle to the mesh to avoid empty meshes, as the
    // mesh library is VERY unhappy with empty meshes
    if (index_buffer_->size() == 0) {
        index_buffer_->push_back(0);
        index_buffer_->push_back(1);
        index_buffer_->push_back(2);
    }
}

bool trialvolume::DualContouring::getExtentCallback(core::Call& caller) {
    auto triangleMeshCall = dynamic_cast<mesh::TriangleMeshCall*>(&caller);
    if (triangleMeshCall != nullptr) {
        triangleMeshCall->set_bounding_box(bbox_);
        triangleMeshCall->set_dimension(mesh::TriangleMeshCall::dimension_t::THREE);
    }
    return true;
}
