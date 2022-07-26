#include "ParticleClusterTracking.h"

#include "datatools/table/TableDataCall.h"
#include "geometry_calls/MultiParticleDataCall.h"

using namespace megamol::trialvolume;

ParticleClusterTracking::ParticleClusterTracking()
        : in_cluster_slot_("in_cluster", "The particle cluster call")
        , out_cluster_track_slot_("out_cluster_track", "The cluster track call") {
    // Setup the input slot
    in_cluster_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_cluster_slot_);

    // Setup the output slot as table call for now
    out_cluster_track_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(0), &ParticleClusterTracking::getDataCallback);
    out_cluster_track_slot_.SetCallback(datatools::table::TableDataCall::ClassName(),
        datatools::table::TableDataCall::FunctionName(1), &ParticleClusterTracking::getExtentCallback);
    MakeSlotAvailable(&out_cluster_track_slot_);
}

ParticleClusterTracking::~ParticleClusterTracking() {}

bool ParticleClusterTracking::create(void) {
    return true;
}

void ParticleClusterTracking::release(void) {}

bool ParticleClusterTracking::getExtentCallback(core::Call& call) {
    computeTracks();
    return true;
}

bool ParticleClusterTracking::getDataCallback(core::Call& call) {
    computeTracks();
    return true;
}

/**
 * Iterate iver akk tune skuces t abd t-1 and figure out the relationship between
 * their clusters. For this we assume that the particles keep their ID between
 * frames.
 */
void ParticleClusterTracking::computeTracks(void) {
    auto* in_cluster_call = in_cluster_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (in_cluster_call == nullptr) {
        return;
    }
    // TODO check that we do not recompute the tracks if we do not have to

    // First reset the current tracks
    cluster_metadata_.clear();
    auto reverse_track_local_id_map = std::map<std::pair<size_t, size_t>, track_unique_id_t>();

    // Then we need to initialize the tracks with the first frame separately,
    // because we do not have a previous frame to compare to.
    in_cluster_call->SetFrameID(0);
    for (auto pl_idx = 0; pl_idx < in_cluster_call->GetParticleListCount(); pl_idx++) {
        auto const& pl = in_cluster_call->AccessParticles(pl_idx);
        auto const& ps = pl.GetParticleStore();
        auto curr_id_acc = ps.GetIDAcc();
        auto curr_cluster_id_acc = ps.GetCRAcc();
        auto curr_x_acc = ps.GetXAcc();
        auto curr_y_acc = ps.GetYAcc();
        auto curr_z_acc = ps.GetZAcc();
        for (auto p_idx = 0; p_idx < pl.GetCount(); p_idx++) {
            auto cluster_id = curr_cluster_id_acc->Get_u32(p_idx);
            auto global_id = reverse_track_local_id_map.find(std::make_pair(0, cluster_id));
            if (global_id == reverse_track_local_id_map.end()) {
                reverse_track_local_id_map[std::make_pair(0, cluster_id)] = cluster_metadata_.size();
                cluster_metadata_.emplace_back();
                cluster_metadata_.back().local_time_cluster_id = cluster_id;
                cluster_metadata_.back().frame_id = 0;
                cluster_metadata_.back().parents.clear();
                cluster_metadata_.back().bounding_box.Set(curr_x_acc->Get_f(p_idx), curr_y_acc->Get_f(p_idx),
                    curr_z_acc->Get_f(p_idx), curr_x_acc->Get_f(p_idx), curr_y_acc->Get_f(p_idx),
                    curr_z_acc->Get_f(p_idx));
            } else {
                auto metadata = &cluster_metadata_[global_id->second];
                metadata->bounding_box.GrowToPoint(
                    curr_x_acc->Get_f(p_idx), curr_y_acc->Get_f(p_idx), curr_z_acc->Get_f(p_idx));
            }
        }
    }

    // Now we need to compare the time steps t and t-1.
    auto prev_in_cluster_call = *in_cluster_call;
    for (unsigned int t = 1; t < in_cluster_call->FrameCount(); t++) {
        in_cluster_call->SetFrameID(t);
        if (!(*in_cluster_call)()) {
            megamol::core::utility::log::Log::DefaultLog.WriteError(
                "[ParticleClusterTracking] Failed to get data for frame %d, skipping it.", t);
            continue;
        }
        megamol::core::utility::log::Log::DefaultLog.WriteInfo(
            "[ParticleClusterTracking] Processing frame %d.", t);
        for (auto pl_idx = 0; pl_idx < prev_in_cluster_call.GetParticleListCount(); pl_idx++) {
            // For now hope that we keep the same amount of particle lists in each time slice
            auto prev_pl = prev_in_cluster_call.AccessParticles(pl_idx);
            auto prev_ps = prev_pl.GetParticleStore();
            auto prev_id_acc = prev_ps.GetIDAcc();
            auto prev_cluster_id_acc = prev_ps.GetCRAcc();

            auto curr_pl = in_cluster_call->AccessParticles(pl_idx);
            auto curr_ps = curr_pl.GetParticleStore();
            auto curr_id_acc = curr_ps.GetIDAcc();
            auto curr_cluster_id_acc = curr_ps.GetCRAcc();
            auto curr_x_acc = curr_ps.GetXAcc();
            auto curr_y_acc = curr_ps.GetYAcc();
            auto curr_z_acc = curr_ps.GetZAcc();
            // Iterate over all particles in the previous time slice and find the
            // corresponding ones in the current time slice
            for (auto i = 0; i < prev_pl.GetCount(); i++) {
                // TODO check that this way of accessing the particles is correct
                // Or should I use i instead of part_id here?
                auto part_id = prev_id_acc->Get_u32(i);
                auto prev_cluster_id = prev_cluster_id_acc->Get_u32(part_id);
                auto curr_cluster_id = curr_cluster_id_acc->Get_u32(part_id);

                // We have now found a link between a previous and current cluster.
                auto new_global_id = reverse_track_local_id_map.find(std::make_pair(t, curr_cluster_id));
                if (new_global_id == reverse_track_local_id_map.end()) {
                    // We have a new cluster
                    auto new_id = cluster_metadata_.size();
                    reverse_track_local_id_map[std::make_pair(t, curr_cluster_id)] = new_id;
                    cluster_metadata_.emplace_back();
                    cluster_metadata_.back().local_time_cluster_id = curr_cluster_id;
                    cluster_metadata_.back().frame_id = t;
                    cluster_metadata_.back().parents = {reverse_track_local_id_map[std::make_pair(t - 1, prev_cluster_id)]};
                    cluster_metadata_.back().bounding_box = vislib::math::Cuboid<float>(curr_x_acc->Get_f(part_id),
                        curr_y_acc->Get_f(part_id), curr_z_acc->Get_f(part_id), curr_x_acc->Get_f(part_id),
                        curr_y_acc->Get_f(part_id), curr_z_acc->Get_f(part_id));
                } else {
                    // The cluster was already seen, update the parents and bounding box if necessary
                    auto metadata = &cluster_metadata_[new_global_id->second];
                    metadata->parents.insert(reverse_track_local_id_map[std::make_pair(t - 1, prev_cluster_id)]);
                    metadata->bounding_box.GrowToPoint(curr_x_acc->Get_f(part_id), curr_y_acc->Get_f(part_id),
                        curr_z_acc->Get_f(part_id));
                }
            }
        }
    }
}
