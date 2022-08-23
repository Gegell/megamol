#include "ParticleClusterTracking.h"

#include "datatools/table/TableDataCall.h"
#include "geometry_calls/MultiParticleDataCall.h"
#include "mmcore/param/ButtonParam.h"
#include "mmcore/param/FilePathParam.h"
#include "mmcore/param/IntParam.h"

#include <fstream>

using namespace megamol::trialvolume;

ParticleClusterTracking::ParticleClusterTracking()
        : in_cluster_slot_("in_cluster", "The particle cluster call")
        , out_cluster_track_slot_("out_cluster_track", "The cluster track call")
        , start_button_("start", "Start tracking")
        , min_connection_count_(
              "min_connection_count", "Minimum particle count for two clusters to be considered connected")
        , dot_file_name_("dot_file_name", "The file name for the .dot file") {
    // Setup the input slot
    in_cluster_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_cluster_slot_);

    // Setup the output slot as table call for now
    out_cluster_track_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(0), &ParticleClusterTracking::getDataCallback);
    out_cluster_track_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(1), &ParticleClusterTracking::getExtentCallback);
    MakeSlotAvailable(&out_cluster_track_slot_);

    // Setup the minimum particle count
    min_connection_count_.SetParameter(new core::param::IntParam(1, 1));
    MakeSlotAvailable(&min_connection_count_);

    // Setup the dot file name
    dot_file_name_.SetParameter(
        new core::param::FilePathParam("out.dot", core::param::FilePathParam::Flag_File_ToBeCreated));
    MakeSlotAvailable(&dot_file_name_);

    // Setup manual start button
    start_button_.SetParameter(new core::param::ButtonParam());
    start_button_.SetUpdateCallback(&ParticleClusterTracking::buttonCallback);
    MakeSlotAvailable(&start_button_);
}

ParticleClusterTracking::~ParticleClusterTracking() {}

bool ParticleClusterTracking::create(void) {
    return true;
}

void ParticleClusterTracking::release(void) {}

bool ParticleClusterTracking::getExtentCallback(core::Call& call) {
    // TODO add the extents of the particle data call
    return true;
}

bool ParticleClusterTracking::getDataCallback(core::Call& call) {
    // TODO actually fill the table data call with the data
    return true;
}

bool ParticleClusterTracking::buttonCallback(core::param::ParamSlot& param) {
    computeTracks();
    generateDotFile();
    return true;
}

/**
 * Iterate over all time slices t and t-1 to figure out the relationship between
 * their clusters. For this we assume that the particles keep their ID between
 * frames.
 */
void ParticleClusterTracking::computeTracks(void) {
    auto* in_cluster_call = in_cluster_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_cluster_call == nullptr) {
        return;
    }
    // TODO check that we do not recompute the tracks if we do not have to

    // Before we do anything we need to call the input cluster call once to
    // get the needed metadata like frame count set as well.
    in_cluster_call->SetFrameID(0, true);
    if (!(*in_cluster_call)(1)) {
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[ParticleClusterTracking] Failed to get extents for frame 0");
        // If we fail here we cannot do anything as we do not know how many
        // frames we have to iterate over.
        return;
    }

    // First reset the current tracks
    cluster_metadata_.clear();

    // Keep track of the previous id accessor here
    typedef std::map<size_t, size_t> particle_cluster_id_map_t;
    std::vector<particle_cluster_id_map_t> cached_cluster_ids;

    // Iterate over every time step t
    for (size_t t = 0; t < in_cluster_call->FrameCount(); t++) {
        // Get the current cluster call
        in_cluster_call->SetFrameID(t);
        bool found_frame_data = true;
        do {
            if (!(*in_cluster_call)(1)) {
                megamol::core::utility::log::Log::DefaultLog.WriteError(
                    "[ParticleClusterTracking] Failed to get extents for frame t");
                found_frame_data = false;
            }
            if (!(*in_cluster_call)(0)) {
                megamol::core::utility::log::Log::DefaultLog.WriteError(
                    "[ParticleClusterTracking] Failed to get cluster call for frame %d, skipping.", t);
                found_frame_data = false;
            }
        } while (in_cluster_call->FrameID() != t);
        if (!found_frame_data) {
            continue;
        }

        // Check that the needed data is here)
        bool passed_checks = in_cluster_call->GetParticleListCount() > 0;
        for (size_t pl_idx = 0; pl_idx < in_cluster_call->GetParticleListCount(); pl_idx++) {
            auto& parts = in_cluster_call->AccessParticles(pl_idx);
            if (parts.GetColourDataType() != geocalls::SimpleSphericalParticles::COLDATA_FLOAT_I &&
                parts.GetColourDataType() != geocalls::SimpleSphericalParticles::COLDATA_DOUBLE_I) {
                core::utility::log::Log::DefaultLog.WriteError(
                    "[ParticleClusterTracking] Particle list %d, frame %d has no intensity colour data, skipping.",
                    pl_idx, t);
                passed_checks = false;
                break;
            }
            if (!parts.HasID()) {
                core::utility::log::Log::DefaultLog.WriteError(
                    "[ParticleClusterTracking] Particle list %d, frame %d has no ID data, skipping.", pl_idx, t);
                passed_checks = false;
                break;
            }
        }
        if (!passed_checks) {
            continue;
        }

        // 1. Find out how many clusters we have in the current time step
        auto min = std::numeric_limits<float>::max();
        auto max = 0.0f;
        for (size_t pl_idx = 0; pl_idx < in_cluster_call->GetParticleListCount(); pl_idx++) {
            auto& parts = in_cluster_call->AccessParticles(pl_idx);
            min = std::min(min, parts.GetMinColourIndexValue());
            max = std::max(max, parts.GetMaxColourIndexValue());
        }
        // 2. Reserve the amount of new clusters we have in the current time step
        auto new_cluster_count = (static_cast<size_t>(max) + 1) - 2;
        auto& current_cluster_list = cluster_metadata_.emplace_back(new_cluster_count);
        for (auto i = 0; i < new_cluster_count; i++) {
            current_cluster_list[i].local_time_cluster_id = i;
            current_cluster_list[i].frame_id = t;
        }

        // Note however that id 0 is reserved for unassigned particles
        // and id 1 is reserved for too sparse particles

        // 3. Iterate over all particles in the current time step
        for (size_t pl_idx = 0; pl_idx < in_cluster_call->GetParticleListCount(); pl_idx++) {
            auto& parts = in_cluster_call->AccessParticles(pl_idx);
            auto ps = parts.GetParticleStore();

            // 3.1. Extract the cluster ID (colour index)
            auto acc_cluster_id = ps.GetCRAcc();
            auto acc_particle_id = ps.GetIDAcc();

            auto acc_x = ps.GetXAcc();
            auto acc_y = ps.GetYAcc();
            auto acc_z = ps.GetZAcc();
            for (size_t p_idx = 0; p_idx < parts.GetCount(); p_idx++) {
                auto c_id = acc_cluster_id->Get_u32(p_idx);
                if (c_id <= 1) {
                    continue;
                }
                c_id -= 2;
                auto& cluster = current_cluster_list[c_id];

                // 3.2. Expand the cluster bounding box to include the particle
                if (cluster.bounding_box.IsEmpty()) {
                    cluster.bounding_box.Set(acc_x->Get_f(p_idx), acc_y->Get_f(p_idx), acc_z->Get_f(p_idx),
                        acc_x->Get_f(p_idx), acc_y->Get_f(p_idx), acc_z->Get_f(p_idx));
                    cluster.bounding_box.Grow(1e-4f);
                } else {
                    cluster.bounding_box.GrowToPoint(acc_x->Get_f(p_idx), acc_y->Get_f(p_idx), acc_z->Get_f(p_idx));
                }
                cluster.num_particles++;

                // 3.3. Mark the cluster as connected to the previous time step
                if (cached_cluster_ids.size() > 0) {
                    // TODO check if particle list size stays the same size between frames
                    auto const & map = cached_cluster_ids[pl_idx];
                    auto const p_id = acc_particle_id->Get_u32(p_idx);
                    auto const matched_id = map.find(p_id);
                    if (matched_id != map.end()) {
                        auto prev_cluster_id = matched_id->second;
                        cluster.parents[prev_cluster_id]++;
                    }
                }
            }
        }

        // 4. Store the previous id accessor for the next time step
        cached_cluster_ids.clear();
        for (size_t pl_idx = 0; pl_idx < in_cluster_call->GetParticleListCount(); pl_idx++) {
            auto& parts = in_cluster_call->AccessParticles(pl_idx);
            auto ps = parts.GetParticleStore();

            auto acc_id = ps.GetIDAcc();
            auto acc_cluster_id = ps.GetCRAcc();
            // Keep the previous id alive for the next time step
            particle_cluster_id_map_t particle_cluster_map;
            for (size_t p_idx = 0; p_idx < parts.GetCount(); p_idx++) {
                auto p_id = acc_id->Get_u32(p_idx);
                auto c_id = acc_cluster_id->Get_u32(p_idx);
                if (c_id <= 1) {
                    continue;
                }
                particle_cluster_map[p_id] = c_id - 2;
            }
            cached_cluster_ids.push_back(particle_cluster_map);
        }

        // 5. Report some statistics
        for (auto& cluster : current_cluster_list) {
            megamol::core::utility::log::Log::DefaultLog.WriteInfo(
                "[ParticleClusterTracking] Cluster (%d,%d) has %d parents, with %d particles.", cluster.frame_id,
                cluster.local_time_cluster_id, cluster.parents.size(), cluster.num_particles);
            for (auto& p : cluster.parents) {
                auto parent_cluster = cluster_metadata_[cluster_metadata_.size() - 2][p.first];
                megamol::core::utility::log::Log::DefaultLog.WriteInfo(
                    "[ParticleClusterTracking] Cluster (%d,%d) has parent (%d,%d) with %d connections.",
                    cluster.frame_id, cluster.local_time_cluster_id, parent_cluster.frame_id,
                    parent_cluster.local_time_cluster_id, p.second);
            }
        }

        // 6. Save the current cluster list, we have partial results should the process be interrupted
        generateDotFile(true);
    }
}

bool ParticleClusterTracking::generateDotFile(bool silent) {
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
                if (p.second < min_connection_count_.Param<core::param::IntParam>()->Value()) {
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
                     << cluster.num_particles << ")\", __bounds=\"[" << cluster.bounding_box.Left() << ","
                     << cluster.bounding_box.Bottom() << "," << cluster.bounding_box.Back() << ","
                     << cluster.bounding_box.Right() << "," << cluster.bounding_box.Top() << ","
                     << cluster.bounding_box.Front() << "]\", __frame=" << cluster.frame_id
                     << ", __local_id=" << cluster.local_time_cluster_id << ", __global_id=" << cluster.num_particles
                     << "];" << std::endl;
        }
    }
    dot_file << "}" << std::endl;
    return true;
}
