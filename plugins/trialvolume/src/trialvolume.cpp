/**
 * MegaMol
 * Copyright (c) 2021, MegaMol Dev Team
 * All rights reserved.
 */

#include "trialvolume/DualContouring.h"
#include "trialvolume/GraphCall.h"
#include "trialvolume/GraphDotWriter.h"
#include "trialvolume/IdentityToICol.h"
#include "trialvolume/MeshSegmentation.h"
#include "trialvolume/MeshSegmentationCall.h"
#include "trialvolume/ParticleClusterTracking.h"
#include "trialvolume/ParticleToVolume.h"
#include "trialvolume/SegmentationAnalysis.h"
#include "trialvolume/VolumeClusterTracking.h"
#include "trialvolume/VolumeSegmentation.h"
#include "trialvolume/VtrFileReader.h"

#include "mmcore/utility/plugins/AbstractPluginInstance.h"
#include "mmcore/utility/plugins/PluginRegister.h"

namespace megamol::trialvolume {
class TrialVolumePluginInstance : public megamol::core::utility::plugins::AbstractPluginInstance {
    REGISTERPLUGIN(TrialVolumePluginInstance)
public:
    TrialVolumePluginInstance()
            : megamol::core::utility::plugins::AbstractPluginInstance(
                  // machine-readable plugin assembly name
                  "trialvolume", // TODO: change to something reasonable later

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
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::VolumeClusterTracking>();
        module_descriptions.RegisterAutoDescription<megamol::trialvolume::GraphDotWriter>();

        // register calls
        call_descriptions.RegisterAutoDescription<megamol::trialvolume::MeshSegmentationCall>();
        call_descriptions.RegisterAutoDescription<megamol::trialvolume::GraphCall>();
    }
};
} // namespace megamol::trialvolume
