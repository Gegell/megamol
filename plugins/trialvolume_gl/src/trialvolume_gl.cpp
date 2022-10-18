/**
 * MegaMol
 * Copyright (c) 2021, MegaMol Dev Team
 * All rights reserved.
 */

#include "ParticleToVolumeGL.h"

#include "mmcore/utility/plugins/AbstractPluginInstance.h"
#include "mmcore/utility/plugins/PluginRegister.h"

namespace megamol::trialvolume_gl {
class TrialVolumeGLPluginInstance : public megamol::core::utility::plugins::AbstractPluginInstance {
    REGISTERPLUGIN(TrialVolumeGLPluginInstance)
public:
    TrialVolumeGLPluginInstance()
            : megamol::core::utility::plugins::AbstractPluginInstance(
                  // machine-readable plugin assembly name
                  "trialvolume_gl", // TODO: change to something reasonable later

                  // human-readable plugin description
                  "Initial trial plugin for bachelor project"){};

    ~TrialVolumeGLPluginInstance() override = default;

    // Registers modules and calls
    void registerClasses() override {

        // register modules
        module_descriptions.RegisterAutoDescription<megamol::trialvolume_gl::ParticleToVolumeGL>();

        // register calls
    }
};
} // namespace megamol::trialvolume_gl
