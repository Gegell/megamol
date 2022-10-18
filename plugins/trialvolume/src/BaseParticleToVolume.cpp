
#include "trialvolume/BaseParticleToVolume.h"

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"

using namespace megamol;
using megamol::core::utility::log::Log;

trialvolume::BaseParticleToVolume::BaseParticleToVolume(void)
        : voxel_size_slot_("voxelSize", "Voxel size")
        , kernel_type_slot_("kernelType", "Kernel type")
        , kernel_metric_slot_("kernelMetric", "Kernel metric")
        , kernel_radius_slot_("kernelRadius", "Kernel radius")
        , kernel_boundary_slot_("kernelBoundary", "Kernel boundary handling")
        , out_density_slot_("outData", "Provides splatted particle volume")
        , out_velocity_slot_("outVelocity", "Provides splatted particle velocity volume")
        , in_particle_data_slot_("inParticleData", "Takes the particle data") {
    // Setup volumetric data output slot
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_DATA),
        &BaseParticleToVolume::getDataCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &BaseParticleToVolume::getExtentCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &BaseParticleToVolume::getExtentCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &BaseParticleToVolume::dummyCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &BaseParticleToVolume::dummyCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &BaseParticleToVolume::dummyCallback);
    MakeSlotAvailable(&out_density_slot_);

    // Setup velocity data output slot
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_DATA),
        &BaseParticleToVolume::getVelocityCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &BaseParticleToVolume::getExtentCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &BaseParticleToVolume::getExtentCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &BaseParticleToVolume::dummyCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &BaseParticleToVolume::dummyCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &BaseParticleToVolume::dummyCallback);
    MakeSlotAvailable(&out_velocity_slot_);


    // Setup particle data input slot
    in_particle_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_particle_data_slot_);

    // Setup voxel size slot
    voxel_size_slot_ << new core::param::FloatParam(
        1.0f, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max(), 0.1f);
    MakeSlotAvailable(&voxel_size_slot_);

    // Setup kernel type slot
    auto* ktp = new core::param::EnumParam(trialvolume::BaseParticleToVolume::KERNEL_TYPE_NEAREST);
    ktp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_TYPE_NEAREST, "Nearest");
    ktp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_TYPE_BUMP, "Bump");
    kernel_type_slot_ << ktp;
    MakeSlotAvailable(&kernel_type_slot_);

    // Setup kernel metric slot
    auto* kmp = new core::param::EnumParam(trialvolume::BaseParticleToVolume::KERNEL_METRIC_EUCLIDEAN);
    kmp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_METRIC_EUCLIDEAN, "Euclidean");
    kmp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_METRIC_MANHATTAN, "Manhattan");
    kmp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_METRIC_CHEBYSHEV, "Chebyshev");
    kernel_metric_slot_ << kmp;
    MakeSlotAvailable(&kernel_metric_slot_);

    // Setup kernel radius slot
    kernel_radius_slot_ << new core::param::FloatParam(
        1.0f, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max(), 0.1f);
    MakeSlotAvailable(&kernel_radius_slot_);

    // Setup kernel boundary slot
    auto* kbp = new core::param::EnumParam(trialvolume::BaseParticleToVolume::KERNEL_BOUNDARY_CLIP);
    kbp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_BOUNDARY_CLIP, "Clip");
    kbp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_BOUNDARY_CLAMP, "Clamp");
    kbp->SetTypePair(trialvolume::BaseParticleToVolume::KERNEL_BOUNDARY_WRAP, "Wrap");
    kernel_boundary_slot_ << kbp;
    MakeSlotAvailable(&kernel_boundary_slot_);
}

trialvolume::BaseParticleToVolume::~BaseParticleToVolume(void) {
    Release();
}

bool trialvolume::BaseParticleToVolume::dummyCallback(core::Call& caller) {
    return true;
}

bool trialvolume::BaseParticleToVolume::getDataCallback(core::Call& caller) {
    auto* outVolumetricDataCall = dynamic_cast<geocalls::VolumetricDataCall*>(&caller);

    if (outVolumetricDataCall == nullptr) {
        return false;
    }

    if (!assertData(*outVolumetricDataCall)) {
        return false;
    }

    outVolumetricDataCall->SetMetadata(&metadata_density_);

    // Set the data
    outVolumetricDataCall->SetFrameID(time_);
    outVolumetricDataCall->SetData(density_.data());
    outVolumetricDataCall->SetDataHash(data_hash_ * 2 + 0);
    outVolumetricDataCall->SetFrameCount(metadata_density_.NumberOfFrames);

    return true;
}

bool trialvolume::BaseParticleToVolume::getVelocityCallback(core::Call& caller) {
    auto* outVolumetricDataCall = dynamic_cast<geocalls::VolumetricDataCall*>(&caller);

    if (outVolumetricDataCall == nullptr) {
        return false;
    }

    if (!assertData(*outVolumetricDataCall)) {
        return false;
    }

    outVolumetricDataCall->SetMetadata(&metadata_velocity_);

    // Set the data
    outVolumetricDataCall->SetFrameID(time_);
    outVolumetricDataCall->SetData(velocity_.data());
    outVolumetricDataCall->SetDataHash(data_hash_ * 2 + 1);
    outVolumetricDataCall->SetFrameCount(metadata_velocity_.NumberOfFrames);

    return true;
}

bool trialvolume::BaseParticleToVolume::assertData(geocalls::VolumetricDataCall& caller) {
    auto* inMultiParticleDataCall = in_particle_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (inMultiParticleDataCall == nullptr) {
        return false;
    }

    // Update until current frame is reached? Copy pasted from ParticlesToDenisty.cpp
    auto frameID = caller.FrameID();
    do {
        inMultiParticleDataCall->SetFrameID(frameID, true);
        if (!(*inMultiParticleDataCall)(1)) {
            Log::DefaultLog.WriteError("%s: Unable to get extents.", ClassName());
            return false;
        }
        if (!(*inMultiParticleDataCall)(0)) {
            Log::DefaultLog.WriteError("%s: Unable to get data.", ClassName());
            return false;
        }
    } while (inMultiParticleDataCall->FrameID() != frameID);
    // Only update if hash is different / parameters changed
    if (time_ != inMultiParticleDataCall->FrameID() || inMultiParticleDataCall->DataHash() != in_data_hash_ ||
        anythingDirty()) {
        if (!createVolume(inMultiParticleDataCall)) {
            return false;
        }
        in_data_hash_ = inMultiParticleDataCall->DataHash();
        data_hash_++;
        time_ = inMultiParticleDataCall->FrameID();
        resetDirtyFlags();
    }
    // Figure out the volume metadata
    auto bbox = inMultiParticleDataCall->AccessBoundingBoxes().ObjectSpaceBBox();

    geocalls::VolumetricDataCall::Metadata temp_meta;

    temp_meta.GridType = geocalls::GridType_t::CARTESIAN;
    temp_meta.NumberOfFrames = inMultiParticleDataCall->FrameCount();
    temp_meta.Components = 0;
    temp_meta.ScalarLength = 0;

    auto const voxelSideLength = voxel_size_slot_.Param<core::param::FloatParam>()->Value();

    temp_meta.Extents[0] = bbox.Width();
    temp_meta.Extents[1] = bbox.Height();
    temp_meta.Extents[2] = bbox.Depth();

    temp_meta.Origin[0] = bbox.Left();
    temp_meta.Origin[1] = bbox.Bottom();
    temp_meta.Origin[2] = bbox.Back();

    temp_meta.Resolution[0] = static_cast<size_t>(std::ceil(temp_meta.Extents[0] / voxelSideLength));
    temp_meta.Resolution[1] = static_cast<size_t>(std::ceil(temp_meta.Extents[1] / voxelSideLength));
    temp_meta.Resolution[2] = static_cast<size_t>(std::ceil(temp_meta.Extents[2] / voxelSideLength));

    temp_meta.SliceDists[0] = new float[0];
    temp_meta.SliceDists[0][0] = temp_meta.Extents[0] / static_cast<float>(temp_meta.Resolution[0] - 1);
    temp_meta.SliceDists[1] = new float[0];
    temp_meta.SliceDists[1][0] = temp_meta.Extents[1] / static_cast<float>(temp_meta.Resolution[1] - 1);
    temp_meta.SliceDists[2] = new float[0];
    temp_meta.SliceDists[2][0] = temp_meta.Extents[2] / static_cast<float>(temp_meta.Resolution[2] - 1);

    temp_meta.IsUniform[0] = true;
    temp_meta.IsUniform[1] = true;
    temp_meta.IsUniform[2] = true;

    metadata_density_ = temp_meta.Clone();
    metadata_density_.Components = 1;
    metadata_density_.ScalarType = geocalls::ScalarType_t::FLOATING_POINT;
    metadata_density_.ScalarLength = sizeof(float);
    metadata_density_.MinValues = new double[1];
    metadata_density_.MinValues[0] = min_value_;
    metadata_density_.MaxValues = new double[1];
    metadata_density_.MaxValues[0] = max_value_;

    metadata_velocity_ = temp_meta.Clone();
    metadata_velocity_.Components = 3;
    metadata_velocity_.ScalarType = geocalls::ScalarType_t::FLOATING_POINT;
    metadata_velocity_.ScalarLength = sizeof(float);
    metadata_velocity_.MinValues = new double[3];
    std::copy(min_velocity_.begin(), min_velocity_.end(), metadata_velocity_.MinValues);
    metadata_velocity_.MaxValues = new double[3];
    std::copy(max_velocity_.begin(), max_velocity_.end(), metadata_velocity_.MaxValues);

    // TODO check that there are no memory leaks here... probably there are with the cloning :/

    return true;
}

bool trialvolume::BaseParticleToVolume::getExtentCallback(core::Call& caller) {
    auto* inMultiParticleDataCall = in_particle_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (inMultiParticleDataCall == nullptr) {
        return false;
    }

    auto* outVolumetricDataCall = dynamic_cast<geocalls::VolumetricDataCall*>(&caller);

    auto frameID = outVolumetricDataCall != nullptr     ? outVolumetricDataCall->FrameID()
                   : inMultiParticleDataCall != nullptr ? inMultiParticleDataCall->FrameID()
                                                        : 0;
    inMultiParticleDataCall->SetFrameID(frameID, true);
    if (!(*inMultiParticleDataCall)(1)) {
        Log::DefaultLog.WriteError("%s: could not get current frame extents (%u)", ClassName(), time_ - 1);
        return false;
    }

    if (outVolumetricDataCall != nullptr) {
        outVolumetricDataCall->AccessBoundingBoxes().SetObjectSpaceBBox(
            inMultiParticleDataCall->GetBoundingBoxes().ObjectSpaceBBox());
        outVolumetricDataCall->AccessBoundingBoxes().SetObjectSpaceClipBox(
            inMultiParticleDataCall->GetBoundingBoxes().ObjectSpaceClipBox());
        outVolumetricDataCall->AccessBoundingBoxes().MakeScaledWorld(1.0f);
        outVolumetricDataCall->SetFrameCount(inMultiParticleDataCall->FrameCount());
    }
    return true;
}

bool trialvolume::BaseParticleToVolume::createVolume(geocalls::MultiParticleDataCall* caller) {
    Log::DefaultLog.WriteInfo("%s: starting volume creation", ClassName());
    const auto startTime = std::chrono::high_resolution_clock::now();

    auto const bbox = caller->AccessBoundingBoxes().ObjectSpaceBBox();
    if (bbox.IsEmpty()) {
        Log::DefaultLog.WriteError("%s: bounding box is empty", ClassName());
        return false;
    }

    auto const voxelSideLength = voxel_size_slot_.Param<core::param::FloatParam>()->Value();

    x_cells_ = static_cast<size_t>(std::ceil(bbox.Width() / voxelSideLength));
    y_cells_ = static_cast<size_t>(std::ceil(bbox.Height() / voxelSideLength));
    z_cells_ = static_cast<size_t>(std::ceil(bbox.Depth() / voxelSideLength));

    density_.resize(x_cells_ * y_cells_ * z_cells_);
    std::fill(density_.begin(), density_.end(), 0.0f);

    velocity_.resize(x_cells_ * y_cells_ * z_cells_ * 3);
    std::fill(velocity_.begin(), velocity_.end(), 0.0f);

    if (!computeVolume(caller)) {
        Log::DefaultLog.WriteError("%s: could not create volume", ClassName());
    }

    auto [min, max] = std::minmax_element(density_.begin(), density_.end());
    min_value_ = *min;
    max_value_ = *max;

    auto [min_vel, max_vel] = std::minmax_element(velocity_.begin(), velocity_.end());
    min_velocity_ = {*min_vel, *min_vel, *min_vel}; // TODO change this to min/max of each component
    max_velocity_ = {*max_vel, *max_vel, *max_vel};

    size_t totalParticles = 0;
    for (size_t i = 0; i < caller->GetParticleListCount(); ++i) {
        auto const& ps = caller->AccessParticles(i);
        totalParticles += ps.GetCount();
    }

    const auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> diffMillis = endTime - startTime;
    Log::DefaultLog.WriteInfo("%s: creation of %u x %u x %u volume from %llu particles took %f ms.", ClassName(),
        x_cells_, y_cells_, z_cells_, totalParticles, diffMillis.count());
    return true;
}
