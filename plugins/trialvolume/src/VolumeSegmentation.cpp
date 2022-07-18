#include "VolumeSegmentation.h"

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/param/FloatParam.h"

using megamol::trialvolume::VolumeSegmentation;

VolumeSegmentation::VolumeSegmentation(void)
        : iso_level_slot_("iso_level", "The relative iso level")
        , out_volume_segment_data_slot_("out_volume_segment_data", "The volume segmentation data")
        , in_volume_data_slot_("in_volume_data", "The volume data") {
    // Setup the input volume data slot
    in_volume_data_slot_.SetCompatibleCall<geocalls::VolumetricDataCallDescription>();
    MakeSlotAvailable(&in_volume_data_slot_);

    // Setup the output volume segmentation data slot
    out_volume_segment_data_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_DATA),
        &VolumeSegmentation::getDataCallback);
    out_volume_segment_data_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &VolumeSegmentation::getExtentCallback);
    out_volume_segment_data_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &VolumeSegmentation::getExtentCallback);
    out_volume_segment_data_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &VolumeSegmentation::dummyCallback);
    out_volume_segment_data_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &VolumeSegmentation::dummyCallback);
    out_volume_segment_data_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &VolumeSegmentation::dummyCallback);
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

bool VolumeSegmentation::getDataCallback(core::Call& call) {
    return computeSegmentation(dynamic_cast<geocalls::VolumetricDataCall&>(call));
}

bool VolumeSegmentation::getExtentCallback(core::Call& call) {
    return computeSegmentation(dynamic_cast<geocalls::VolumetricDataCall&>(call));
}

bool VolumeSegmentation::dummyCallback(core::Call& call) {
    return true;
}

bool VolumeSegmentation::computeSegmentation(geocalls::VolumetricDataCall& call) {
    return false;
}