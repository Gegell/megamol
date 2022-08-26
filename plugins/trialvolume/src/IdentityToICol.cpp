#include "IdentityToICol.h"

megamol::trialvolume::IdentityToICol::IdentityToICol(void)
        : datatools::AbstractParticleManipulator("outData", "indata") {}

megamol::trialvolume::IdentityToICol::~IdentityToICol(void) {
    this->Release();
}

bool megamol::trialvolume::IdentityToICol::manipulateData(
    geocalls::MultiParticleDataCall& outData, geocalls::MultiParticleDataCall& inData) {

    outData = inData;                   // also transfers the unlocker to 'outData'
    inData.SetUnlocker(nullptr, false); // keep original data locked
                                        // original data will be unlocked through outData

    auto const plc = outData.GetParticleListCount();
    for (unsigned int i = 0; i < plc; ++i) {
        auto& p = outData.AccessParticles(i);

        if (!p.HasID()) {
            megamol::core::utility::log::Log::DefaultLog.WriteWarn("IdentityToICol: Particlelist %d has no ID\n", i);
            continue;
        }

        auto const idt = p.GetIDDataType();

        intensities.resize(p.GetCount());
        auto const basePtrOut = intensities.data();
        auto& idAcc = p.GetParticleStore().GetIDAcc();

        for (size_t pidx = 0; pidx < p.GetCount(); ++pidx) {
            if (idt == geocalls::SimpleSphericalParticles::IDDATA_UINT32) {
                basePtrOut[pidx] = static_cast<float>(idAcc->Get_u32(pidx));
            } else {
                basePtrOut[pidx] = static_cast<float>(idAcc->Get_u64(pidx));
            }
        }

        p.SetColourData(geocalls::SimpleSphericalParticles::COLDATA_FLOAT_I, intensities.data());
        p.SetColourMapIndexValues(0, 1);
    }
    return true;
}
