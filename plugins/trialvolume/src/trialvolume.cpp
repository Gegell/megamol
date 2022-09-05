/**
 * MegaMol
 * Copyright (c) 2021, MegaMol Dev Team
 * All rights reserved.
 */

#include "DualContouring.h"
#include "ParticleClusterTracking.h"
#include "ParticleToVolume.h"
#include "VtrFileReader.h"
#include "VolumeSegmentation.h"
#include "MeshSegmentationCall.h"
#include "MeshSegmentation.h"
#include "SegmentationAnalysis.h"
#include "IdentityToICol.h"

#include "mmcore/utility/plugins/AbstractPluginInstance.h"
#include "mmcore/utility/plugins/PluginRegister.h"

namespace megamol::trialvolume {
class TrialVolumePluginInstance : public megamol::core::utility::plugins::AbstractPluginInstance {
    REGISTERPLUGIN(TrialVolumePluginInstance)
public:
    TrialVolumePluginInstance()
            : megamol::core::utility::plugins::AbstractPluginInstance(
                  // machine-readable plugin assembly name
                  "trailvolume", // TODO: change to something reasonable later

                  // human-readable plugin description
                  "Initial trial plugin for bachelor project"){};

    ~TrialVolumePluginInstance() override = default;

    // Registers modules and calls
    void registerClasses() override {

        // register modules
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::ParticleToVolume>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::VtrFileReader>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::DualContouring>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::ParticleClusterTracking>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::IdentityToICol>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::VolumeSegmentation>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::MeshSegmentation>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::SegmentationAnalysis>();

        // register calls
        call_descriptions.RegisterAutoDescription<megamol::trialvolume::MeshSegmentationCall>();
    }
};
} // namespace megamol::trialvolume
