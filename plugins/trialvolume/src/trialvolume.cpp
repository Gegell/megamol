/**
 * MegaMol
 * Copyright (c) 2021, MegaMol Dev Team
 * All rights reserved.
 */

#include "mmcore/utility/plugins/AbstractPluginInstance.h"
#include "mmcore/utility/plugins/PluginRegister.h"

#include "ParticleToVolume.h"
#include "vtrFileReader.h"

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
        this->module_descriptions.RegisterAutoDescription<megamol::trialvolume::ParticleToVolume>();
        this->module_descriptions.RegisterAutoDescription<megamol::trialvolume::vtrFileReader>();
        //
        // TODO: Register your plugin's modules here:
        // this->module_descriptions.RegisterAutoDescription<megamol::trialvolume::MyModule1>();
        // this->module_descriptions.RegisterAutoDescription<megamol::trialvolume::MyModule2>();
        // ...
        //

        // register calls

        // TODO: Register your plugin's calls here:
        // this->call_descriptions.RegisterAutoDescription<megamol::trialvolume::MyCall1>();
        // this->call_descriptions.RegisterAutoDescription<megamol::trialvolume::MyCall2>();
        // ...
        //
    }
};
} // namespace megamol::trialvolume
