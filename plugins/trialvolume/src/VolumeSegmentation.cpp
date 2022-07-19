#include "VolumeSegmentation.h"

#include <stack>
#include <vector>

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
    auto* in_volume_data_call = in_volume_data_slot_.CallAs<geocalls::VolumetricDataCall>();
    if (in_volume_data_call == nullptr) {
        return false;
    }

    in_volume_data_call->SetFrameID(call.FrameID());
    if (!(*in_volume_data_call)(1)) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("[VolumeSegmentation] Failed to get volume extents");
        return false;
    }
    if (!(*in_volume_data_call)(0)) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("[VolumeSegmentation] Failed to get volume data");
        return false;
    }

    auto metadata = in_volume_data_call->GetMetadata();
    if (metadata->Components != 1 || metadata->ScalarType != geocalls::ScalarType_t::FLOATING_POINT) {
        // TODO this might be subject to change, if the volume also contains additional data such as the velocity field, and not only the density.
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "[VolumeSegmentation] Volume data must be single value, float density component");
        return false;
    }

    // Check if the volume data changed
    if (in_volume_data_call->DataHash() != input_data_hash_ || anythingDirty()) {
        input_data_hash_ = in_volume_data_call->DataHash();
        resetDirtyFlags();

        // Compute the absolute iso level from the relative iso level input
        float iso_level = iso_level_slot_.Param<core::param::FloatParam>()->Value();
        auto threshold = metadata->MinValues[0] + iso_level * (metadata->MaxValues[0] - metadata->MinValues[0]);
        // TODO might need to be replaced with only the absolute iso level instead of relative
        // as the min and max values are not necessarily the same for different frames

        // Fetch the volume data
        auto* volume_data = static_cast<float*>(in_volume_data_call->GetData());

        // Initialize the volume segmentation data like the volume data
        segment_ids_.resize(metadata->Resolution[0] * metadata->Resolution[1] * metadata->Resolution[2], 0);
        segment_count_ = 0;

        // Iterate over all indices and start region growing if the density is above the iso level
        for (size_t z = 0, index = 0; z < metadata->Resolution[2]; z++) {
            for (size_t y = 0; y < metadata->Resolution[1]; y++) {
                for (size_t x = 0; x < metadata->Resolution[0]; x++, index++) {
                    if (segment_ids_[index] == 0 && volume_data[index] > threshold) {
                        segment_count_++;
                        regionGrow(x, y, z, *metadata, volume_data, threshold, segment_count_);
                    }
                }
            }
        }

        megamol::core::utility::log::Log::DefaultLog.WriteInfo("[VolumeSegmentation] Found %lld segments", segment_count_);
    }

    // Copy most of the previous metadata attributes to the given call
    segment_metadata_ = metadata->Clone();
    segment_metadata_.MinValues[0] = 0.0f;
    segment_metadata_.MaxValues[0] = static_cast<double>(segment_count_);
    segment_metadata_.ScalarType = geocalls::ScalarType_t::UNSIGNED_INTEGER;
    segment_metadata_.ScalarLength = sizeof(uint16_t);
    
    // Set the segmentation data
    call.SetMetadata(&segment_metadata_);
    call.SetFrameID(in_volume_data_call->FrameID());
    call.SetData(segment_ids_.data());
    call.SetDataHash(input_data_hash_);

    return true;
}

void VolumeSegmentation::regionGrow(const size_t x, const size_t y, const size_t z,
    const geocalls::VolumetricMetadata_t& meta, const float* data, const float threshold, const uint16_t label) {
    typedef std::tuple<size_t, size_t, size_t> cellPos_t;

    auto to_check = std::stack<cellPos_t>();
    to_check.push(std::make_tuple(x, y, z));
    while (!to_check.empty()) {
        auto [x, y, z] = to_check.top();
        to_check.pop();
        size_t index = x + y * meta.Resolution[0] + z * meta.Resolution[0] * meta.Resolution[1];
        if (segment_ids_[index] != label && data[index] > threshold) {
            segment_ids_[index] = label;
            if (z > 0) {
                to_check.push(std::make_tuple(x, y, z - 1));
            }
            if (z < meta.Resolution[2] - 1) {
                to_check.push(std::make_tuple(x, y, z + 1));
            }
            if (y > 0) {
                to_check.push(std::make_tuple(x, y - 1, z));
            }
            if (y < meta.Resolution[1] - 1) {
                to_check.push(std::make_tuple(x, y + 1, z));
            }
            if (x > 0) {
                to_check.push(std::make_tuple(x - 1, y, z));
            }
            if (x < meta.Resolution[0] - 1) {
                to_check.push(std::make_tuple(x + 1, y, z));
            }
        }
    }
}
