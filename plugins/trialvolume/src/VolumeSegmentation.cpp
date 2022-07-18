#include "VolumeSegmentation.h"

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/param/FloatParam.h"

using namespace megamol::trialvolume;
using namespace megamol;

VolumeSegmentation::VolumeSegmentation()
        : iso_level_slot_("iso_level", "The relative iso level"),
        , out_volume_segment_data_slot_("out_volume_segment_data", "The volume segmentation data")
        , in_volume_data_slot_("in_volume_data", "The volume data") {
    // Setup the input volume data slot
    in_volume_data_slot_.SetCompatilbleCall<geocalls::VolumetricDataCallDescription>();
    MakeSlotAvailable(&in_volume_data_slot_);

    // Setup the output volume segmentation data slot
    out_volume_segment_data_slot_.SetCallback(geocalls::ClassName(), "GetData", &VolumeSegmentation::getData);
    out_volume_segment_data_slot_.SetCallback(geocalls::ClassName(), "GetExtent", &VolumeSegmentation::getExtent);
    MakeSlotAvailable(&out_volume_segment_data_slot_);

    // Setup the iso level slot
    iso_level_slot_.SetParameter(new core::param::FloatParam(0.5f, 0.0f, 1.0f));
    MakeSlotAvailable(&iso_level_slot_);
}

VolumeSegmentation::~VolumeSegmentation() {
    release();
}

bool VolumeSegmentation::create() {
    return true;
}

void VolumeSegmentation::release() {
    // TODO release any data here
}
