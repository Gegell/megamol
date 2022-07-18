#pragma once

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/CalleeSlot.h"
#include "mmcore/CallerSlot.h"
#include "mmcore/Module.h"
#include "mmcore/param/ParamSlot.h"

namespace megamol::trialvolume {

class VolumeSegmentation : public core::Module {
public:

    static const char* ClassName() {
        return "VolumeSegmentation";
    }

    static const char* Description() {
        return "Performs segmentation on a volume to generate a mask of connected parts.";
    }

    static bool IsAvailable() {
        return true;
    }

    VolumeSegmentation();

    ~VolumeSegmentation() override;

private:

    bool create() override;

    void release() override;

    inline bool anythingDirty() const {
        return false;
    }

    inline void resetDirtyFlags() {
    }

    bool getDataCallback(core::Call& call);

    bool getExtentCallback(core::Call& call);

    bool dummyCallback(core::Call& call);

    bool computeSegmentation(geocalls::VolumetricDataCall& call);

    /** The slot for the volume data */
    core::CallerSlot in_volume_data_slot_;

    /** The slot for the volume segmentation data */
    core::CalleeSlot out_volume_segment_data_slot_;

    /** The slot for the iso level */
    core::param::ParamSlot iso_level_slot_;

    /** The absolute iso level */
    float absolute_iso_level_;

    /** The hash of when the input data was last read */
    size_t input_data_hash_;
};

} // namespace megamol::trialvolume