#include "IdentityToICol.h"

megamol::trialvolume::IdentityToICol::IdentityToICol(void)
        : datatools::AbstractParticleManipulator("outData", "indata") {}

megamol::trialvolume::IdentityToICol::~IdentityToICol(void) {
    this->Release();
}

bool megamol::trialvolume::IdentityToICol::manipulateData(
    geocalls::MultiParticleDataCall& outData, geocalls::MultiParticleDataCall& inData) {

    // Only copy the id data if the hash of the data has changed
    if (inData.DataHash() == last_hash_) {
        return true;
    }
    last_hash_ = inData.DataHash();

    outData = inData;                   // also transfers the unlocker to 'outData'
    inData.SetUnlocker(nullptr, false); // keep original data locked
                                        // original data will be unlocked through outData

    // Determine the number of particles
    auto const particle_list_count = outData.GetParticleListCount();
    auto particle_count = 0u;
    for (size_t i = 0u; i < particle_list_count; ++i) {
        particle_count += outData.AccessParticles(i).GetCount();
    }

    // Allocate memory for the new particle color intensity data
    intensities_.clear();
    intensities_.reserve(particle_count);

    for (unsigned int i = 0; i < particle_list_count; ++i) {
        auto& parts = outData.AccessParticles(i);

        if (!parts.HasID()) {
            megamol::core::utility::log::Log::DefaultLog.WriteWarn("IdentityToICol: Particlelist %d has no ID\n", i);
            continue;
        }

        auto const idType = parts.GetIDDataType();
        auto const& idAcc = parts.GetParticleStore().GetIDAcc();
        auto const* address = intensities_.data() + intensities_.size();

        auto min_intensity = std::numeric_limits<float>::max();
        auto max_intensity = std::numeric_limits<float>::min();
        for (size_t pidx = 0; pidx < parts.GetCount(); ++pidx) {
            if (idType == geocalls::SimpleSphericalParticles::IDDATA_UINT32) {
                intensities_.push_back(static_cast<float>(idAcc->Get_u32(pidx)));
            } else {
                intensities_.push_back(static_cast<float>(idAcc->Get_u64(pidx)));
            }
            min_intensity = std::min(min_intensity, intensities_.back());
            max_intensity = std::max(max_intensity, intensities_.back());
        }

        parts.SetColourData(geocalls::SimpleSphericalParticles::COLDATA_FLOAT_I, address);
        parts.SetColourMapIndexValues(min_intensity, max_intensity);
    }
    return true;
}
