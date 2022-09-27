
#include "ParticleToVolume.h"

#include <functional>

#include <voro++.hh>

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/param/EnumParam.h"
#include "mmcore/param/FloatParam.h"

using namespace megamol;

bool trialvolume::ParticleToVolume::create(void) {
    return true;
}

void trialvolume::ParticleToVolume::release(void) {
    // TODO release any data here
}

trialvolume::ParticleToVolume::ParticleToVolume(void)
        : splatting_method_slot_("SplattingMethod", "The splatting method to use")
        , voxel_size_slot_("voxelSize", "Voxel size")
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
        &ParticleToVolume::getDataCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &ParticleToVolume::getExtentCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &ParticleToVolume::getExtentCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &ParticleToVolume::dummyCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &ParticleToVolume::dummyCallback);
    out_density_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &ParticleToVolume::dummyCallback);
    MakeSlotAvailable(&out_density_slot_);

    // Setup velocity data output slot
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_DATA),
        &ParticleToVolume::getVelocityCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &ParticleToVolume::getExtentCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &ParticleToVolume::getExtentCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &ParticleToVolume::dummyCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &ParticleToVolume::dummyCallback);
    out_velocity_slot_.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &ParticleToVolume::dummyCallback);
    MakeSlotAvailable(&out_velocity_slot_);


    // Setup particle data input slot
    in_particle_data_slot_.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    MakeSlotAvailable(&in_particle_data_slot_);

    // Setup voxel size slot
    voxel_size_slot_ << new core::param::FloatParam(
        1.0f, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max(), 0.1f);
    MakeSlotAvailable(&voxel_size_slot_);

    // Setup kernel type slot
    auto* ktp = new core::param::EnumParam(trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST);
    ktp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST, "Nearest");
    ktp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_TYPE_BUMP, "Bump");
    kernel_type_slot_ << ktp;
    MakeSlotAvailable(&kernel_type_slot_);

    // Setup kernel metric slot
    auto* kmp = new core::param::EnumParam(trialvolume::ParticleToVolume::KERNEL_METRIC_EUCLIDEAN);
    kmp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_METRIC_EUCLIDEAN, "Euclidean");
    kmp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_METRIC_MANHATTAN, "Manhattan");
    kmp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_METRIC_CHEBYSHEV, "Chebyshev");
    kernel_metric_slot_ << kmp;
    MakeSlotAvailable(&kernel_metric_slot_);

    // Setup kernel radius slot
    kernel_radius_slot_ << new core::param::FloatParam(
        1.0f, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max(), 0.1f);
    MakeSlotAvailable(&kernel_radius_slot_);

    // Setup kernel boundary slot
    auto* kbp = new core::param::EnumParam(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLIP);
    kbp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLIP, "Clip");
    kbp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLAMP, "Clamp");
    kbp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_WRAP, "Wrap");
    kernel_boundary_slot_ << kbp;
    MakeSlotAvailable(&kernel_boundary_slot_);

    // Setup splatting method slot
    auto* spm = new core::param::EnumParam(trialvolume::ParticleToVolume::SPLAT_METHOD_KERNEL);
    spm->SetTypePair(trialvolume::ParticleToVolume::SPLAT_METHOD_KERNEL, "Kernel");
    spm->SetTypePair(trialvolume::ParticleToVolume::SPLAT_METHOD_NATURAL_NEIGHBOR, "Natural Neighbor");
    splatting_method_slot_ << spm;
    MakeSlotAvailable(&splatting_method_slot_);
}

trialvolume::ParticleToVolume::~ParticleToVolume(void) {
    Release();
}

bool trialvolume::ParticleToVolume::dummyCallback(core::Call& caller) {
    return true;
}

bool trialvolume::ParticleToVolume::getDataCallback(core::Call& caller) {
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

bool trialvolume::ParticleToVolume::getVelocityCallback(core::Call& caller) {
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

bool trialvolume::ParticleToVolume::assertData(geocalls::VolumetricDataCall& caller) {
    auto* inMultiParticleDataCall = in_particle_data_slot_.CallAs<geocalls::MultiParticleDataCall>();
    if (inMultiParticleDataCall == nullptr) {
        return false;
    }

    // Update until current frame is reached? Copy pasted from ParticlesToDenisty.cpp
    auto frameID = caller.FrameID();
    do {
        inMultiParticleDataCall->SetFrameID(frameID, true);
        if (!(*inMultiParticleDataCall)(1)) {
            megamol::core::utility::log::Log::DefaultLog.WriteError("ParticleToVolume: Unable to get extents.");
            return false;
        }
        if (!(*inMultiParticleDataCall)(0)) {
            megamol::core::utility::log::Log::DefaultLog.WriteError("ParticleToVolume: Unable to get data.");
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

bool trialvolume::ParticleToVolume::getExtentCallback(core::Call& caller) {
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
        // What is this ? Why pass function <1>?
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "ParticlesToDensity: could not get current frame extents (%u)", time_ - 1);
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

bool trialvolume::ParticleToVolume::createVolume(geocalls::MultiParticleDataCall* caller) {
    megamol::core::utility::log::Log::DefaultLog.WriteInfo("ParticleToVolume: starting volume creation");
    const auto startTime = std::chrono::high_resolution_clock::now();
    size_t totalParticles = 0;

    auto const bbox = caller->AccessBoundingBoxes().ObjectSpaceBBox();
    if (bbox.IsEmpty()) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("ParticleToVolume: bounding box is empty");
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

    auto success = false;
    switch (splatting_method_slot_.Param<core::param::EnumParam>()->Value()) {
    case trialvolume::ParticleToVolume::SPLAT_METHOD_KERNEL:
        success = computeKernel(caller);
        break;
    case trialvolume::ParticleToVolume::SPLAT_METHOD_NATURAL_NEIGHBOR:
        success = computeNaturalNeighborhood(caller);
        break;
    }
    if (!success) {
        megamol::core::utility::log::Log::DefaultLog.WriteError("ParticleToVolume: could not create volume");
    }

    auto [min, max] = std::minmax_element(density_.begin(), density_.end());
    min_value_ = *min;
    max_value_ = *max;

    auto [min_vel, max_vel] = std::minmax_element(velocity_.begin(), velocity_.end());
    min_velocity_ = {*min_vel, *min_vel, *min_vel}; // TODO change this to min/max of each component
    max_velocity_ = {*max_vel, *max_vel, *max_vel};

    const auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> diffMillis = endTime - startTime;
    megamol::core::utility::log::Log::DefaultLog.WriteInfo(
        "ParticleToVolume: creation of %u x %u x %u volume from %llu particles took %f ms.", x_cells_, y_cells_,
        z_cells_, totalParticles, diffMillis.count());
    return true;
}

bool trialvolume::ParticleToVolume::computeKernel(geocalls::MultiParticleDataCall* caller) {
    auto const bbox = caller->AccessBoundingBoxes().ObjectSpaceBBox();

    auto const voxelSideLength = voxel_size_slot_.Param<core::param::FloatParam>()->Value();

    auto const kernelRadius = kernel_radius_slot_.Param<core::param::FloatParam>()->Value();
    auto const kernelCellSpan = static_cast<int>(std::ceil(kernelRadius / voxelSideLength));

    std::function<float(float, float, float)> lengthFunction;
    switch (kernel_metric_slot_.Param<core::param::EnumParam>()->Value()) {
    default:
    case trialvolume::ParticleToVolume::KERNEL_METRIC_EUCLIDEAN:
        lengthFunction = [](float const x, float const y, float const z) -> float {
            return std::sqrt(x * x + y * y + z * z);
        };
        break;
    case trialvolume::ParticleToVolume::KERNEL_METRIC_MANHATTAN:
        lengthFunction = [](float const x, float const y, float const z) -> float {
            return std::abs(x) + std::abs(y) + std::abs(z);
        };
        break;
    case trialvolume::ParticleToVolume::KERNEL_METRIC_CHEBYSHEV:
        lengthFunction = [](float const x, float const y, float const z) -> float {
            return std::max(std::abs(x), std::max(std::abs(y), std::abs(z)));
        };
        break;
    }

    std::function<float(float)> kernel;
    switch (kernel_type_slot_.Param<core::param::EnumParam>()->Value()) {
    default:
    case trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST:
        kernel = [](float const dist) -> float { return 0.0f; };
        break;
    case trialvolume::ParticleToVolume::KERNEL_TYPE_BUMP:
        kernel = [kernelRadius](float const dist) -> float {
            return dist <= kernelRadius ? std::exp(-1.0f / (1.0f - std::pow(dist / kernelRadius, 2.0f))) : 0.0f;
        };
        break;
    }

    std::function<float(float)> applyBoundary;
    switch (kernel_boundary_slot_.Param<core::param::EnumParam>()->Value()) {
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLAMP:
        applyBoundary = [](float const value) -> float { return std::max(0.0f, std::min(value, 1.0f)); };
        break;
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_WRAP:
        applyBoundary = [](float const value) -> float { return value - std::floor(value); };
        break;
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLIP:
        applyBoundary = [](float const value) -> float { return value; };
        break;
    }

    for (size_t i = 0; i < caller->GetParticleListCount(); i++) {
        auto& particleList = caller->AccessParticles(i);
        auto& ps = particleList.GetParticleStore();

        auto xAcc = ps.GetXAcc();
        auto yAcc = ps.GetYAcc();
        auto zAcc = ps.GetZAcc();

        auto const xDirAcc = ps.GetDXAcc();
        auto const yDirAcc = ps.GetDYAcc();
        auto const zDirAcc = ps.GetDZAcc();

        for (size_t j = 0; j < particleList.GetCount(); j++) {
            auto x = xAcc->Get_f(j);
            auto y = yAcc->Get_f(j);
            auto z = zAcc->Get_f(j);

            auto xNorm = (x - bbox.Left()) / bbox.Width();
            auto yNorm = (y - bbox.Bottom()) / bbox.Height();
            auto zNorm = (z - bbox.Back()) / bbox.Depth();

            auto isKernel = kernel_type_slot_.Param<core::param::EnumParam>()->Value() !=
                            trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST;
            if (!isKernel) {
                auto const xBounded = static_cast<size_t>(std::round(applyBoundary(xNorm) * x_cells_));
                auto const yBounded = static_cast<size_t>(std::round(applyBoundary(yNorm) * y_cells_));
                auto const zBounded = static_cast<size_t>(std::round(applyBoundary(zNorm) * z_cells_));

                // Check if we are inside the volume
                if (xBounded < 0 || xBounded >= x_cells_ || yBounded < 0 || yBounded >= y_cells_ || zBounded < 0 ||
                    zBounded >= z_cells_) {
                    continue;
                }

                auto const index = (zBounded * y_cells_ + yBounded) * x_cells_ + xBounded;
                density_[index] += 1.0f;

                velocity_[index * 3 + 0] += xDirAcc->Get_f(j);
                velocity_[index * 3 + 1] += yDirAcc->Get_f(j);
                velocity_[index * 3 + 2] += zDirAcc->Get_f(j);
            } else {
                for (auto dz = -kernelCellSpan; dz <= kernelCellSpan; ++dz) {
                    for (auto dy = -kernelCellSpan; dy <= kernelCellSpan; ++dy) {
                        for (auto dx = -kernelCellSpan; dx <= kernelCellSpan; ++dx) {
                            auto const xBounded = static_cast<size_t>(
                                std::round(applyBoundary(xNorm + static_cast<float>(dx) / x_cells_) * x_cells_));
                            auto const yBounded = static_cast<size_t>(
                                std::round(applyBoundary(yNorm + static_cast<float>(dy) / y_cells_) * y_cells_));
                            auto const zBounded = static_cast<size_t>(
                                std::round(applyBoundary(zNorm + static_cast<float>(dz) / z_cells_) * z_cells_));

                            // Check if we are inside the volume
                            if (xBounded < 0 || xBounded >= x_cells_ || yBounded < 0 || yBounded >= y_cells_ ||
                                zBounded < 0 || zBounded >= z_cells_) {
                                continue;
                            }

                            auto const index = (zBounded * y_cells_ + yBounded) * x_cells_ + xBounded;
                            // FIXME use offset to cell vertex
                            auto const dist =
                                lengthFunction(dx * voxelSideLength, dy * voxelSideLength, dz * voxelSideLength);
                            auto const weight = kernel(dist);

                            density_[index] += weight;

                            velocity_[index * 3 + 0] += xDirAcc->Get_f(j) * weight;
                            velocity_[index * 3 + 1] += yDirAcc->Get_f(j) * weight;
                            velocity_[index * 3 + 2] += zDirAcc->Get_f(j) * weight;
                        }
                    }
                }
            }
        }
    }

    // Normalize velocity
    for (size_t i = 0; i < x_cells_ * y_cells_ * z_cells_; i++) {
        if (density_[i] > 0.0f) {
            velocity_[i * 3 + 0] /= density_[i];
            velocity_[i * 3 + 1] /= density_[i];
            velocity_[i * 3 + 2] /= density_[i];
        }
    }

    return true;
}

bool trialvolume::ParticleToVolume::computeNaturalNeighborhood(geocalls::MultiParticleDataCall* caller) {
    auto const bbox = caller->AccessBoundingBoxes().ObjectSpaceBBox();

    auto isWrapping = false;
    switch (kernel_boundary_slot_.Param<core::param::EnumParam>()->Value()) {
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLAMP:
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "ParticleToVolume: Clamp boundary not supported for natural neighborhood");
        return false;
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_WRAP:
        isWrapping = true;
        break;
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLIP:
        isWrapping = false;
        break;
    }
    auto voroContainer = voro::container(bbox.Left(), bbox.Right(), bbox.Bottom(), bbox.Top(), bbox.Back(),
        bbox.Front(), 8, 8, 8, isWrapping, isWrapping, isWrapping, 8);
    auto particleId = 0;
    for (size_t i = 0; i < caller->GetParticleListCount(); i++) {
        auto& particleList = caller->AccessParticles(i);
        auto& ps = particleList.GetParticleStore();

        auto xAcc = ps.GetXAcc();
        auto yAcc = ps.GetYAcc();
        auto zAcc = ps.GetZAcc();

        for (size_t j = 0; j < particleList.GetCount(); j++, particleId++) {
            auto x = xAcc->Get_f(j);
            auto y = yAcc->Get_f(j);
            auto z = zAcc->Get_f(j);

            voroContainer.put(particleId, x, y, z);
        }
    }

    std::vector<int> neighbors;
    std::vector<double> weights;

    auto const kernelRadius = kernel_radius_slot_.Param<core::param::FloatParam>()->Value();
    std::function<float(float)> kernel;
    switch (kernel_type_slot_.Param<core::param::EnumParam>()->Value()) {
    case trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST:
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "ParticleToVolume: Nearest neighbor not supported for natural neighborhood");
        return false;
    default:
    case trialvolume::ParticleToVolume::KERNEL_TYPE_BUMP:
        kernel = [kernelRadius](float const dist) -> float {
            return dist <= kernelRadius ? std::exp(-1.0f / (1.0f - std::pow(dist / kernelRadius, 2.0f))) : 0.0f;
        };
        break;
    }

    for (auto z = 0; z < z_cells_; ++z)
        for (auto y = 0; y < y_cells_; ++y)
            for (auto x = 0; x < x_cells_; ++x) {

                auto xNorm = x / static_cast<double>(x_cells_ - 1);
                auto yNorm = y / static_cast<double>(y_cells_ - 1);
                auto zNorm = z / static_cast<double>(z_cells_ - 1);

                auto xLocal = xNorm * bbox.Width() + bbox.Left();
                auto yLocal = yNorm * bbox.Height() + bbox.Bottom();
                auto zLocal = zNorm * bbox.Depth() + bbox.Back();

                voro::voronoicell_neighbor cell(voroContainer);
                if (voroContainer.compute_ghost_cell(cell, xLocal, yLocal, zLocal)) {

                    auto const index = (z * y_cells_ + y) * x_cells_ + x;
                    cell.face_areas(weights);
                    cell.neighbors(neighbors);
                    auto weightSum = 0.0;
                    auto interpolated = 0.0;
                    auto interpolatedVelocity = std::array<double, 3>{0.0, 0.0, 0.0};

                    for (size_t i = 0; i < neighbors.size(); i++) {
                        auto neighborId = neighbors[i];
                        // Skip negative neighbors, e.g. if the neighbor is the boundary
                        if (neighborId < 0)
                            continue;

                        // Search for the particle list containing the neighbor
                        auto particleListId = 0;
                        for (; particleListId < caller->GetParticleListCount(); particleListId++) {
                            auto& particleList = caller->AccessParticles(particleListId);
                            if (neighborId >= particleList.GetCount()) {
                                neighborId -= particleList.GetCount();
                                continue;
                            } else {
                                break;
                            }
                        }
                        auto& particleList = caller->AccessParticles(particleListId);
                        auto& ps = particleList.GetParticleStore();

                        // Compute the distance between the two particles
                        auto xAcc = ps.GetXAcc();
                        auto yAcc = ps.GetYAcc();
                        auto zAcc = ps.GetZAcc();
                        auto const xNeighbor = xAcc->Get_f(neighborId);
                        auto const yNeighbor = yAcc->Get_f(neighborId);
                        auto const zNeighbor = zAcc->Get_f(neighborId);
                        auto const dx = xNeighbor - xLocal;
                        auto const dy = yNeighbor - yLocal;
                        auto const dz = zNeighbor - zLocal;
                        auto const dist = std::sqrt(dx * dx + dy * dy + dz * dz);

                        // Compute the un-normalized laplacian weight
                        weights[i] /= dist;
                        weightSum += weights[i];

                        // Interpolate the distance to the neighbor
                        // TODO change this to use actual point data (e.g. color, etc.)
                        interpolated += weights[i] * kernel(dist);

                        // Interpolate the velocity to the neighbor
                        auto const xDirAcc = ps.GetDXAcc();
                        auto const yDirAcc = ps.GetDYAcc();
                        auto const zDirAcc = ps.GetDZAcc();
                        auto const xDir = xDirAcc->Get_f(neighborId);
                        auto const yDir = yDirAcc->Get_f(neighborId);
                        auto const zDir = zDirAcc->Get_f(neighborId);
                        interpolatedVelocity[0] += weights[i] * xDir;
                        interpolatedVelocity[1] += weights[i] * yDir;
                        interpolatedVelocity[2] += weights[i] * zDir;
                    }
                    // Normalize the laplacian weights
                    interpolated /= weightSum;

                    // Normalize the velocity
                    interpolatedVelocity[0] /= weightSum;
                    interpolatedVelocity[1] /= weightSum;
                    interpolatedVelocity[2] /= weightSum;

                    // Apply the rbf
                    density_[index] = interpolated;

                    // Apply the velocity
                    velocity_[index * 3 + 0] = interpolatedVelocity[0];
                    velocity_[index * 3 + 1] = interpolatedVelocity[1];
                    velocity_[index * 3 + 2] = interpolatedVelocity[2];
                }
            }
    return true;
}
