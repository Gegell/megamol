#include "DualContouring.h"

#include <vector>

#include <Eigen/Dense>

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

Eigen::Vector3f trialvolume::DualContouring::forwardDiffNormalAtVolumePoint(
    size_t x, size_t y, size_t z, const geocalls::VolumetricDataCall& volumeDataCall) {
    // Clamp the indices to the volume dimensions
    x = std::clamp(x, 0ull, volumeDataCall.GetMetadata()->Resolution[0] - 2u);
    y = std::clamp(y, 0ull, volumeDataCall.GetMetadata()->Resolution[1] - 2u);
    z = std::clamp(z, 0ull, volumeDataCall.GetMetadata()->Resolution[2] - 2u);

    auto const vol = volumeDataCall.GetAbsoluteVoxelValue(x, y, z);

    Eigen::Vector3f normal;
    normal << volumeDataCall.GetAbsoluteVoxelValue(x + 1, y, z) - vol,
        volumeDataCall.GetAbsoluteVoxelValue(x, y + 1, z) - vol,
        volumeDataCall.GetAbsoluteVoxelValue(x, y, z + 1) - vol;
    normal.normalize();
    return normal;
}

Eigen::Vector3f trialvolume::DualContouring::findBestVertex(
    const size_t x, const size_t y, const size_t z, const geocalls::VolumetricDataCall& volumeDataCall) {
    // HACK just return the center of the cell for now. :(
    // return Eigen::Vector3f(x + 0.5f, y + 0.5f, z + 0.5f);

    auto isoLevel = iso_level_slot_.Param<core::param::FloatParam>()->Value();

    Eigen::MatrixX3f edge_surface_points(12, 3);
    Eigen::MatrixX3f normals(12, 3);

    auto total_edges = 0;
    // Iterate over all edges of the voxel
    for (auto dy = 0; dy <= 1; ++dy) {
        for (auto dx = 0; dx <= 1; ++dx) {
            auto level1 = volumeDataCall.GetRelativeVoxelValue(x + dx, y + dy, z) - isoLevel;
            auto level2 = volumeDataCall.GetRelativeVoxelValue(x + dx, y + dy, z + 1) - isoLevel;
            if (!sameSign(level1, level2)) {
                // Calculate the intersection point
                auto isoCrossing = zeroCrossingLocation(level1, level2);

                // Approximate the normal on said point
                Eigen::Vector3f normal1 = forwardDiffNormalAtVolumePoint(x + dx, y + dy, z, volumeDataCall);
                Eigen::Vector3f normal2 = forwardDiffNormalAtVolumePoint(x + dx, y + dy, z + 1, volumeDataCall);

                Eigen::Vector3f pos = Eigen::Vector3f(x + dx, y + dy, z) + isoCrossing * Eigen::Vector3f(0, 0, 1);
                // Eigen::Vector3f normal = isoCrossing * normal1 + (1.0f - isoCrossing) * normal2;
                // normal.normalize();

                edge_surface_points.row(total_edges) = pos;
                normals.row(total_edges) = normal1;
                total_edges++;
            }
        }
    }
    for (auto dz = 0; dz <= 1; ++dz) {
        for (auto dx = 0; dx <= 1; ++dx) {
            auto level1 = volumeDataCall.GetRelativeVoxelValue(x + dx, y, z + dz) - isoLevel;
            auto level2 = volumeDataCall.GetRelativeVoxelValue(x + dx, y + 1, z + dz) - isoLevel;
            if (!sameSign(level1, level2)) {
                // Calculate the intersection point
                auto isoCrossing = zeroCrossingLocation(level1, level2);

                // Approximate the normal on said point
                auto normal1 = forwardDiffNormalAtVolumePoint(x + dx, y, z + dz, volumeDataCall);
                auto normal2 = forwardDiffNormalAtVolumePoint(x + dx, y + 1, z + dz, volumeDataCall);

                auto pos = Eigen::Vector3f(x + dx, y, z + dz) + isoCrossing * Eigen::Vector3f(0, 1, 0);
                // auto normal = isoCrossing * normal1 + (1.0f - isoCrossing) * normal2;
                // normal.normalize();

                edge_surface_points.row(total_edges) = pos;
                normals.row(total_edges) = normal1;
                total_edges++;
            }
        }
    }
    for (auto dz = 0; dz <= 1; ++dz) {
        for (auto dy = 0; dy <= 1; ++dy) {
            auto level1 = volumeDataCall.GetRelativeVoxelValue(x, y + dy, z + dz) - isoLevel;
            auto level2 = volumeDataCall.GetRelativeVoxelValue(x + 1, y + dy, z + dz) - isoLevel;
            if (!sameSign(level1, level2)) {
                // Calculate the intersection point
                auto isoCrossing = zeroCrossingLocation(level1, level2);

                // Approximate the normal on said point
                auto normal1 = forwardDiffNormalAtVolumePoint(x, y + dy, z + dz, volumeDataCall);
                auto normal2 = forwardDiffNormalAtVolumePoint(x + 1, y + dy, z + dz, volumeDataCall);

                auto pos = Eigen::Vector3f(x, y + dy, z + dz) + isoCrossing * Eigen::Vector3f(1, 0, 0);
                // auto normal = isoCrossing * normal1 + (1.0f - isoCrossing) * normal2;
                // normal.normalize();

                edge_surface_points.row(total_edges) = pos;
                normals.row(total_edges) = normal1;
                total_edges++;
            }
        }
    }

    if (total_edges == 0) {
        return Eigen::Vector3f(x + 0.5f, y + 0.5f, z + 0.5f);
    }
    edge_surface_points.conservativeResize(total_edges, Eigen::NoChange);
    normals.conservativeResize(total_edges, Eigen::NoChange);

    // Solve the least squares problem
    Eigen::VectorXf b(edge_surface_points.rows());
    for (auto i = 0; i < edge_surface_points.rows(); ++i) {
        b(i) = edge_surface_points.row(i).dot(normals.row(i));
    }

    // auto pos = normals.fullPivHouseholderQr().solve(b);
    // BUG currently outputting a nan/inf position from the solver.
    // I don't quite get why this is happening. Seperate CAS showed it should be resolvable / not explode.

    auto pos = (normals.transpose() * normals).ldlt().solve(normals.transpose() * b);

    auto constrained = pos.cwiseMax(Eigen::Vector3i(x, y, z).cast<float>()).cwiseMin(Eigen::Vector3i(x + 1, y + 1, z + 1).cast<float>());

    // std::cout << "pos: " << pos.transpose() << "\tb: " << b.transpose() << std::endl;
    // std::cout << "normals: " << normals.transpose() << std::endl;
    return constrained;
}

float trialvolume::DualContouring::zeroCrossingLocation(const float level1, const float level2) {
    return (0.0f - level1) / (level2 - level1);
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

                auto pos = findBestVertex(x, y, z, volumeDataCall);
                vertex_buffer_->push_back(
                    pos(0) * metadata->Extents[0] / metadata->Resolution[0] + metadata->Origin[0]);
                vertex_buffer_->push_back(
                    pos(1) * metadata->Extents[1] / metadata->Resolution[1] + metadata->Origin[1]);
                vertex_buffer_->push_back(
                    pos(2) * metadata->Extents[2] / metadata->Resolution[2] + metadata->Origin[2]);

                // FIXME TEMPORARY: Compute the normals via forward differencing
                auto normal = forwardDiffNormalAtVolumePoint(x, y, z, volumeDataCall);
                normal_buffer_->push_back(normal(0));
                normal_buffer_->push_back(normal(1));
                normal_buffer_->push_back(normal(2));
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
