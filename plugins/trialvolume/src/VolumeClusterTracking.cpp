#include "VolumeClusterTracking.h"

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/param/BoolParam.h"
#include "mmcore/param/ButtonParam.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/IntParam.h"

#include <algorithm>
#include <fstream>

using namespace megamol::trialvolume;

VolumeClusterTracking::VolumeClusterTracking()
        : in_cluster_id_slot_("in_cluster_id", "The cluster id data")
        , in_velocity_slot_("in_velocity", "The velocity data")
        , out_cluster_track_slot_("out_cluster_id", "The cluster tracking data")
        , start_button_("start", "Start tracking")
        , time_step_size_("time_step_size", "The time step size")
        , min_connection_count_("min_connection_count", "The minimum amount of mass needed to form a connection.")
        , frame_range_limit_param_("frame::limit_range", "Limit the frame range in which to track the clusters")
        , frame_start_param_("frame::start", "The first frame to track")
        , frame_end_param_("frame::end", "The last frame to track")
        , frame_step_param_("frame::step", "The step size for the frame range")
        , dot_file_name_("dot_file_name", "The file name for the .dot file") {
    // Setup the input slots
    in_cluster_id_slot_.SetCompatibleCall<geocalls::VolumetricDataCallDescription>();
    MakeSlotAvailable(&in_cluster_id_slot_);

    in_velocity_slot_.SetCompatibleCall<geocalls::VolumetricDataCallDescription>();
    MakeSlotAvailable(&in_velocity_slot_);

    // Setup the output slots
    // TODO: This needs to be its own call transporting the graph data

    // Setup the minimum connection mass
    min_connection_count_.SetParameter(new core::param::FloatParam(1, 0));
    MakeSlotAvailable(&min_connection_count_);

    // Setup the time step size
    time_step_size_.SetParameter(new core::param::FloatParam(1, 0));
    MakeSlotAvailable(&time_step_size_);

    // Setup the frame range
    frame_range_limit_param_.SetParameter(new core::param::BoolParam(false));
    frame_start_param_.SetParameter(new core::param::IntParam(0, 0));
    frame_end_param_.SetParameter(new core::param::IntParam(0, 0));
    frame_step_param_.SetParameter(new core::param::IntParam(1, 1));
    MakeSlotAvailable(&frame_range_limit_param_);
    MakeSlotAvailable(&frame_start_param_);
    MakeSlotAvailable(&frame_end_param_);
    MakeSlotAvailable(&frame_step_param_);

    // Setup the dot file name
    dot_file_name_.SetParameter(
        new core::param::FilePathParam("out.dot", core::param::FilePathParam::Flag_File_ToBeCreated));
    MakeSlotAvailable(&dot_file_name_);

    // Setup manual start button
    start_button_.SetParameter(new core::param::ButtonParam());
    start_button_.SetUpdateCallback(&VolumeClusterTracking::buttonCallback);
    MakeSlotAvailable(&start_button_);
}

VolumeClusterTracking::~VolumeClusterTracking() {}

bool VolumeClusterTracking::create() {
    return true;
}

void VolumeClusterTracking::release() {}

bool VolumeClusterTracking::getDataCallback(core::Call& call) {
    assert(false); // TODO: Implement this
    return true;
}

bool VolumeClusterTracking::getExtentCallback(core::Call& call) {
    assert(false); // TODO: Implement this
    return true;
}

bool VolumeClusterTracking::buttonCallback(core::param::ParamSlot& slot) {
    computeTracks();
    generateDotFile();
    return true;
}

void VolumeClusterTracking::computeTracks() {
    // Check if the necessary data is available
    auto const cluster_id_call = in_cluster_id_slot_.CallAs<geocalls::VolumetricDataCall>();
    if (cluster_id_call == nullptr) {
        core::utility::log::Log::DefaultLog.WriteError("[VolumeClusterTracking] No cluster id data available.");
        return;
    }

    auto const velocity_call = in_velocity_slot_.CallAs<geocalls::VolumetricDataCall>();
    if (velocity_call == nullptr) {
        core::utility::log::Log::DefaultLog.WriteError("[VolumeClusterTracking] No velocity data available.");
        return;
    }

    // Capture the frame range by calling  the extents call
    cluster_id_call->SetFrameID(0, true);
    if (!(*cluster_id_call)(1)) {
        core::utility::log::Log::DefaultLog.WriteError(
            "[VolumeClusterTracking] Unable to fetch cluster id extents for frame 0.");
        return;
    }
    velocity_call->SetFrameID(0, true);
    if (!(*velocity_call)(1)) {
        core::utility::log::Log::DefaultLog.WriteError(
            "[VolumeClusterTracking] Unable to fetch velocity extents for frame 0.");
        return;
    }

    // Check that both calls correspond to one another (at least in terms of the frame count)
    if (cluster_id_call->FrameCount() != velocity_call->FrameCount()) {
        core::utility::log::Log::DefaultLog.WriteError(
            "[VolumeClusterTracking] The cluster id and velocity data have mismatching frame counts.");
        return;
    }
    // Check that both have the same dimensions
    auto const cluster_id_resolution = cluster_id_call->GetMetadata()->Resolution;
    auto const velocity_resolution = velocity_call->GetMetadata()->Resolution;
    for (auto i = 0; i < 3; ++i) {
        if (cluster_id_resolution[i] != velocity_resolution[i]) {
            core::utility::log::Log::DefaultLog.WriteError("[VolumeClusterTracking] The cluster id and velocity data "
                                                           "have mismatching resolutions in dimension %d. (%d vs. %d)",
                i, cluster_id_resolution[i], velocity_resolution[i]);
            return;
        }
    }
    // Check that the data types are correct
#if 0
    if (!(cluster_id_call->GetMetadata()->ScalarType == geocalls::VolumetricDataCall::ScalarType::UNSIGNED_INTEGER ||
            cluster_id_call->GetMetadata()->ScalarType == geocalls::VolumetricDataCall::ScalarType::SIGNED_INTEGER)) {
        core::utility::log::Log::DefaultLog.WriteError(
            "[VolumeClusterTracking] The cluster id data is not of type (unsigned) integer.");
        return;
    }
#endif
    if (velocity_call->GetMetadata()->ScalarType != geocalls::VolumetricDataCall::ScalarType::FLOATING_POINT) {
        core::utility::log::Log::DefaultLog.WriteError(
            "[VolumeClusterTracking] The velocity data is not of type float.");
        return;
    } else if (velocity_call->GetMetadata()->Components != 3) {
        core::utility::log::Log::DefaultLog.WriteError(
            "[VolumeClusterTracking] The velocity data does not have 3 components.");
        return;
    }


    // Possibly limit the frame range
    unsigned int frame_start, frame_end, frame_step;
    if (frame_range_limit_param_.Param<core::param::BoolParam>()->Value()) {
        frame_start = frame_start_param_.Param<core::param::IntParam>()->Value();
        frame_end = frame_end_param_.Param<core::param::IntParam>()->Value();
        frame_step = frame_step_param_.Param<core::param::IntParam>()->Value();
    } else {
        frame_start = 0;
        frame_end = cluster_id_call->FrameCount() - 1;
        frame_step = 1;
    }
    frame_end = std::min(frame_end, cluster_id_call->FrameCount() - 1);

    auto const dt = time_step_size_.Param<core::param::FloatParam>()->Value();

    // Precompute the grid offsets for equirectangular grids
    std::array<std::vector<float>, 3> grid_points;
    for (auto i = 0; i < 3; i++) {
        float total = 0.0f;
        for (auto j = 0; j < cluster_id_resolution[i]; j++) {
            grid_points[i].push_back(total);
            switch (cluster_id_call->GetMetadata()->GridType)
            {
            case geocalls::VolumetricDataCall::GridType::RECTILINEAR:
                total += cluster_id_call->GetMetadata()->SliceDists[i][j];
                break;
            case geocalls::VolumetricDataCall::GridType::CARTESIAN:
                total += cluster_id_call->GetMetadata()->SliceDists[i][0];
                break;
            default:
                core::utility::log::Log::DefaultLog.WriteError(
                    "[VolumeClusterTracking] Unsupported grid type.");
                return;
            }
        }
        grid_points[i].push_back(total);
    }

    // First reset the current tracks
    cluster_metadata_.clear();

    auto const index_extent = cluster_id_resolution[0] * cluster_id_resolution[1] * cluster_id_resolution[2];
    std::vector<size_t> prev_cluster_ids(index_extent, 0);
    std::vector<ClusterMetadata_t> prev_cluster_list;

    // Iterate over every time step t
    for (auto t = frame_start; t <= frame_end; t += frame_step) {
        // Fetch the current cluster id and velocity data
        bool found_frame_data = true;
        do {
            cluster_id_call->SetFrameID(t, true);
            if (!(*cluster_id_call)(1)) {
                core::utility::log::Log::DefaultLog.WriteError(
                    "[VolumeClusterTracking] Unable to fetch cluster id extents for frame %u.", t);
                found_frame_data = false;
                break;
            }
            if (!(*cluster_id_call)(0)) {
                core::utility::log::Log::DefaultLog.WriteError(
                    "[VolumeClusterTracking] Unable to fetch cluster id data for frame %u.", t);
                found_frame_data = false;
                break;
            }
        } while (cluster_id_call->FrameID() != t);
        do {
            velocity_call->SetFrameID(t, true);
            if (!(*velocity_call)(1)) {
                core::utility::log::Log::DefaultLog.WriteError(
                    "[VolumeClusterTracking] Unable to fetch velocity extents for frame %u.", t);
                found_frame_data = false;
                break;
            }
            if (!(*velocity_call)(0)) {
                core::utility::log::Log::DefaultLog.WriteError(
                    "[VolumeClusterTracking] Unable to fetch velocity data for frame %u.", t);
                found_frame_data = false;
                break;
            }
        } while (velocity_call->FrameID() != t);
        if (!found_frame_data) {
            continue;
        }

        // For each cell perform a single euler step to compute the new position
        // and compare the cluster id at the new position in the next time step
        // to the cluster id at the current position in the current time step.

        // 1. Find out how many clusters we have in the current time step
        auto const cluster_count = static_cast<size_t>(cluster_id_call->GetMetadata()->MaxValues[0]);

        // 2. Reserve the amount of new clusters we have in the current time step
        auto& current_cluster_list = cluster_metadata_.emplace_back(cluster_count);
        for (auto i = 0; i < cluster_count; i++) {
            current_cluster_list[i].local_time_cluster_id = i;
            current_cluster_list[i].frame_id = t;
        }

        auto const vel_ptr = static_cast<float*>(velocity_call->GetData());

        // 3. Estimate where each cell will be in the next time step, connect with future id.
        for (size_t z = 0; z < cluster_id_resolution[2]; ++z) {
            for (size_t y = 0; y < cluster_id_resolution[1]; ++y) {
                for (size_t x = 0; x < cluster_id_resolution[0]; ++x) {
                    auto const current_index =
                        x + y * cluster_id_resolution[0] + z * cluster_id_resolution[0] * cluster_id_resolution[1];

                    auto const x_vel = vel_ptr[current_index * 3 + 0];
                    auto const y_vel = vel_ptr[current_index * 3 + 1];
                    auto const z_vel = vel_ptr[current_index * 3 + 2];

                    // Compute the new position
                    auto const x_pos = grid_points[0][x] - x_vel * dt * frame_step;
                    auto const y_pos = grid_points[1][y] - y_vel * dt * frame_step;
                    auto const z_pos = grid_points[2][z] - z_vel * dt * frame_step;

                    // Find the new position in the grid (use std::lower_bound)
                    auto const x_it = std::lower_bound(grid_points[0].begin(), grid_points[0].end(), x_pos);
                    auto const y_it = std::lower_bound(grid_points[1].begin(), grid_points[1].end(), y_pos);
                    auto const z_it = std::lower_bound(grid_points[2].begin(), grid_points[2].end(), z_pos);

                    auto const new_x = std::distance(grid_points[0].begin(), x_it);
                    auto const new_y = std::distance(grid_points[1].begin(), y_it);
                    auto const new_z = std::distance(grid_points[2].begin(), z_it);


                    // Check if the new position is inside the volume
                    if (new_x < 0 || new_x >= cluster_id_resolution[0] || new_y < 0 ||
                        new_y >= cluster_id_resolution[1] || new_z < 0 || new_z >= cluster_id_resolution[2]) {
                        continue;
                    }

                    // Check if the cluster id at the new position is the same as the cluster id at the current position
                    auto const current_cluster_id = static_cast<size_t>(cluster_id_call->GetAbsoluteVoxelValue(new_x, new_y, new_z));
                    auto const prev_cluster_id = prev_cluster_ids[current_index];


                    if (current_cluster_id != 0) {
                        auto& current_cluster = current_cluster_list[current_cluster_id - 1];
                        auto const cell_pos = vislib::math::Vector<float, 3>(grid_points[0][x], grid_points[1][y], grid_points[2][z]);

                        if (current_cluster.total_mass == 0) {
                            current_cluster.bounding_box.Set(cell_pos.X(), cell_pos.Y(), cell_pos.Z(), cell_pos.X(), cell_pos.Y(), cell_pos.Z());
                        } else {
                            current_cluster.bounding_box.GrowToPoint(cell_pos.X(), cell_pos.Y(), cell_pos.Z());
                        }

                        current_cluster.center_of_mass += cell_pos;
                        current_cluster.velocity += vislib::math::Vector<float, 3>(x_vel, y_vel, z_vel);
                        current_cluster.total_mass += 1.0f; // TODO: Use density instead of 1.0f

                        if (t > frame_start && prev_cluster_id != 0) {
                            current_cluster.parents[prev_cluster_id - 1] += 1; // TODO: Use density instead of 1
                        }
                    }
                }
            }
        }

        // 4. Normalize the center of mass and velocity
        for (auto& cluster : current_cluster_list) {
            cluster.center_of_mass /= cluster.total_mass;
            auto const origin = cluster_id_call->GetMetadata()->Origin;
            cluster.center_of_mass += vislib::math::Vector<float, 3>(origin[0], origin[1], origin[2]);
            cluster.velocity /= cluster.total_mass;
        }

        // 5. Copy the next cluster ids to the current cluster ids
        // HACK this cast to a float pointer is still scuffed as heck, but works for my case now.
        auto const id_ptr = static_cast<float*>(cluster_id_call->GetData());
        prev_cluster_ids.clear();
        prev_cluster_ids.reserve(index_extent);
        for (size_t i = 0; i < index_extent; i++) {
            prev_cluster_ids.push_back(static_cast<size_t>(id_ptr[i]));
        }

        // 6. Report some statistics
        core::utility::log::Log::DefaultLog.WriteInfo(
            "[VolumeClusterTracking] Frame %u: %u clusters", t, cluster_count);

        // 7. Write the cluster data to the output file
        generateDotFile(true);
    }
}


bool VolumeClusterTracking::generateDotFile(bool silent) {
    auto filename = dot_file_name_.Param<core::param::FilePathParam>()->Value();
    if (filename.empty()) {
        if (!silent) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[ParticleClusterTracking] No filename specified for dot file.");
        }
        return false;
    }
    std::ofstream dot_file(filename.generic_u8string().c_str());
    dot_file << "digraph G {" << std::endl;

    // Iterate over all time steps
    for (size_t t = 0; t < cluster_metadata_.size(); t++) {
        auto& cluster_list = cluster_metadata_[t];
        // Iterate over all clusters in the current time step
        for (auto& cluster : cluster_list) {
            // Iterate over all parents of the current cluster and write the edges
            // connecting the current cluster to the previous time steps clusters
            for (auto& p : cluster.parents) {
                // Discard clusters with too few connections
                if (p.second < min_connection_count_.Param<core::param::FloatParam>()->Value()) {
                    continue;
                }
                auto parent_cluster = cluster_metadata_[t - 1][p.first];
                dot_file << "\"" << parent_cluster.frame_id << "_" << parent_cluster.local_time_cluster_id << "\"";
                dot_file << " -> ";
                dot_file << "\"" << cluster.frame_id << "_" << cluster.local_time_cluster_id << "\"";
                dot_file << " [label=\"" << p.second << "\", __kept=" << p.second << "];" << std::endl;
            }
            // Output the current cluster as a node
            dot_file << "\"" << cluster.frame_id << "_" << cluster.local_time_cluster_id << "\"";
            dot_file << " [label=\"" << cluster.frame_id << "_" << cluster.local_time_cluster_id << " ("
                     << cluster.total_mass << ")\", __bounds=\"[" << cluster.bounding_box.Left() << ","
                     << cluster.bounding_box.Bottom() << "," << cluster.bounding_box.Back() << ","
                     << cluster.bounding_box.Right() << "," << cluster.bounding_box.Top() << ","
                     << cluster.bounding_box.Front() << "]\""
                     << ", __center_of_mass=\"[" << cluster.center_of_mass.X() << "," << cluster.center_of_mass.Y()
                     << "," << cluster.center_of_mass.Z() << "]\""
                     << ", __velocity=\"[" << cluster.velocity.X() << "," << cluster.velocity.Y() << ","
                     << cluster.velocity.Z() << "]\""
                     << ", __frame=" << cluster.frame_id << ", __local_id=" << cluster.local_time_cluster_id
                     << ", __total_mass=" << cluster.total_mass << "];" << std::endl;
        }
    }
    dot_file << "}" << std::endl;
    return true;
}
