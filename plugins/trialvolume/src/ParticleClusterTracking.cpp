#include "ParticleClusterTracking.h"

#include "datatools/table/TableDataCall.h"
#include "geometry_calls/MultiParticleDataCall.h"
#include "mmcore/param/ButtonParam.h"

using namespace megamol::trialvolume;

ParticleClusterTracking::ParticleClusterTracking()
        : in_cluster_slot_("in_cluster", "The particle cluster call")
        , out_cluster_track_slot_("out_cluster_track", "The cluster track call")
        , start_button_("start", "Start tracking") {
    // Setup the input slot
    in_cluster_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_cluster_slot_);

    // Setup the output slot as table call for now
    out_cluster_track_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(0), &ParticleClusterTracking::getDataCallback);
    out_cluster_track_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(1), &ParticleClusterTracking::getExtentCallback);
    MakeSlotAvailable(&out_cluster_track_slot_);

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
    std::vector<std::shared_ptr<geocalls::Accessor>> acc_prev_cluster_id;

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
        auto new_cluster_count = static_cast<size_t>(max) + 1;
        auto cur_list_offset = cluster_metadata_.size();
        for (auto i = 0; i < new_cluster_count; i++) {
            cluster_metadata_.emplace_back();
            cluster_metadata_.back().frame_id = t;
            cluster_metadata_.back().local_time_cluster_id = i;
        }

        // Note however that id 0 is reserved for unassigned particles
        // and id 1 is reserved for too sparse particles

        // 3. Iterate over all particles in the current time step and
        for (size_t pl_idx = 0; pl_idx < in_cluster_call->GetParticleListCount(); pl_idx++) {
            auto& parts = in_cluster_call->AccessParticles(pl_idx);
            auto ps = parts.GetParticleStore();

            // 3.1. Extract the cluster ID (colour index)
            auto acc_cluster_id = ps.GetCRAcc();

            // 3.2. Expand the cluster bounding box to include the particle
            auto acc_x = ps.GetXAcc();
            auto acc_y = ps.GetYAcc();
            auto acc_z = ps.GetZAcc();
            for (size_t p_idx = 0; p_idx < parts.GetCount(); p_idx++) {
                auto id = acc_cluster_id->Get_u32(p_idx);
                if (id <= 1) {
                    continue;
                }
                auto& cluster = cluster_metadata_[cur_list_offset + id];
                // TODO this bbox still has the 0,0,0 point unjustly included
                if (cluster.bounding_box.IsEmpty()) {
                    cluster.bounding_box.Set(acc_x->Get_f(p_idx), acc_y->Get_f(p_idx), acc_z->Get_f(p_idx),
                        acc_x->Get_f(p_idx), acc_y->Get_f(p_idx), acc_z->Get_f(p_idx));
                    cluster.bounding_box.Grow(1e-4f);
                } else {
                    cluster.bounding_box.GrowToPoint(acc_x->Get_f(p_idx), acc_y->Get_f(p_idx), acc_z->Get_f(p_idx));
                }
                cluster.num_particles++;

                // 3.3. Mark the cluster as connected to the previous time step
                if (acc_prev_cluster_id.size() > 0) {
                    // TODO check if particle list size stays the same size between frames
                    auto prev_id = acc_prev_cluster_id[pl_idx]->Get_u32(p_idx);
                    if (prev_id > 1) {
                        cluster.parents.insert(cur_list_offset + prev_id);
                    }
                }
            }
        }

        // 4. Store the previous id accessor for the next time step
        acc_prev_cluster_id.clear();
        for (size_t pl_idx = 0; pl_idx < in_cluster_call->GetParticleListCount(); pl_idx++) {
            auto& parts = in_cluster_call->AccessParticles(pl_idx);
            auto ps = parts.GetParticleStore();
            acc_prev_cluster_id.push_back(ps.GetCRAcc());
        }

        // 5. Report some statistics
        for (auto c_idx = cur_list_offset; c_idx < cluster_metadata_.size(); c_idx++) {
            auto& cluster = cluster_metadata_[c_idx];
            if (cluster.parents.size() > 0) {
                megamol::core::utility::log::Log::DefaultLog.WriteInfo("[ParticleClusterTracking] Cluster %d has %d parents.", c_idx, cluster.parents.size());
            }
        }
    }
}
