
#include "ParticleToVolume.h"

#include "geometry_calls/VolumetricDataCall.h"
#include "mmcore/param/FloatParam.h"
#include "mmcore/param/EnumParam.h"

#include <functional>

using namespace megamol;

bool trialvolume::ParticleToVolume::create(void) {
    return true;
}

void trialvolume::ParticleToVolume::release(void) {
    // TODO release any data here
}

trialvolume::ParticleToVolume::ParticleToVolume(void) 
        : voxelSizeSlot("voxelSize", "Voxel size")
        , kernelTypeSlot("kernelType", "Kernel type")
        , kernelMetricSlot("kernelMetric", "Kernel metric")
        , kernelRadiusSlot("kernelRadius", "Kernel radius")
        , kernelBoundarySlot("kernelBoundary", "Kernel boundary handling")
        , outDataSlot("outData", "Provides splatted particle volume")
        , inParticleDataSlot("inParticleData", "Takes the particle data") {
    // Setup volumetric data output slot
    this->outDataSlot.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_DATA),
        &ParticleToVolume::getDataCallback);
    this->outDataSlot.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_EXTENTS),
        &ParticleToVolume::getExtentCallback);
    this->outDataSlot.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_GET_METADATA),
        &ParticleToVolume::getExtentCallback);
    this->outDataSlot.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_START_ASYNC),
        &ParticleToVolume::dummyCallback);
    this->outDataSlot.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_STOP_ASYNC),
        &ParticleToVolume::dummyCallback);
    this->outDataSlot.SetCallback(geocalls::VolumetricDataCall::ClassName(),
        geocalls::VolumetricDataCall::FunctionName(geocalls::VolumetricDataCall::IDX_TRY_GET_DATA),
        &ParticleToVolume::dummyCallback);
    this->MakeSlotAvailable(&this->outDataSlot);

    // Setup particle data input slot
    this->inParticleDataSlot.SetCompatibleCall<geocalls::MultiParticleDataCallDescription>();
    this->MakeSlotAvailable(&this->inParticleDataSlot);

    // Setup voxel size slot
    this->voxelSizeSlot << new core::param::FloatParam(
        1.0f, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max(), 0.1f);
    this->MakeSlotAvailable(&this->voxelSizeSlot);

    // Setup kernel type slot
    auto* ktp = new core::param::EnumParam(trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST);
    ktp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST, "Nearest");
    ktp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_TYPE_BUMP, "Bump");
    this->kernelTypeSlot << ktp;
    this->MakeSlotAvailable(&this->kernelTypeSlot);

    // Setup kernel metric slot
    auto* kmp = new core::param::EnumParam(trialvolume::ParticleToVolume::KERNEL_METRIC_EUCLIDEAN);
    kmp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_METRIC_EUCLIDEAN, "Euclidean");
    kmp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_METRIC_MANHATTAN, "Manhattan");
    kmp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_METRIC_CHEBYSHEV, "Chebyshev");
    this->kernelMetricSlot << kmp;
    this->MakeSlotAvailable(&this->kernelMetricSlot);

    // Setup kernel radius slot
    this->kernelRadiusSlot << new core::param::FloatParam(
        1.0f, std::numeric_limits<float>::epsilon(), std::numeric_limits<float>::max(), 0.1f);
    this->MakeSlotAvailable(&this->kernelRadiusSlot);

    // Setup kernel boundary slot
    auto* kbp = new core::param::EnumParam(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLIP);
    kbp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLIP, "Clip");
    kbp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLAMP, "Clamp");
    kbp->SetTypePair(trialvolume::ParticleToVolume::KERNEL_BOUNDARY_WRAP, "Wrap");
    this->kernelBoundarySlot << kbp;
    this->MakeSlotAvailable(&this->kernelBoundarySlot);

    // TODO: Extend to include splatting kernel selection, domain wrapping, etc.
}

trialvolume::ParticleToVolume::~ParticleToVolume(void) {
    this->Release();
}

bool trialvolume::ParticleToVolume::dummyCallback(core::Call& caller) {
    return true;
}

bool trialvolume::ParticleToVolume::getDataCallback(core::Call& caller) {
    auto* inMultiParticleDataCall = this->inParticleDataSlot.CallAs<geocalls::MultiParticleDataCall>();
    if (inMultiParticleDataCall == nullptr) {
        return false;
    }

    auto* outVolumetricDataCall = dynamic_cast<geocalls::VolumetricDataCall*>(&caller);
    if (outVolumetricDataCall != nullptr) { 
#if 1
        // Update until current frame is reached? Copy pasted from ParticlesToDenisty.cpp
        auto frameID = outVolumetricDataCall != nullptr ? outVolumetricDataCall->FrameID() : 0;
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
#endif
        // Only update if hash is different / parameters changed
        if (this->time != inMultiParticleDataCall->FrameID() 
            || inMultiParticleDataCall->DataHash() != this->inDataHash
            || this->anythingDirty()
            ) {
            if (!this->createVolume(inMultiParticleDataCall)) {
                return false;
            }
            this->inDataHash = inMultiParticleDataCall->DataHash();
            this->dataHash++;
            this->time = inMultiParticleDataCall->FrameID();
            this->resetDirtyFlags();
        }
        // Figure out the volume metadata
        auto bbox = inMultiParticleDataCall->AccessBoundingBoxes().ObjectSpaceBBox();

        metadata.Components = 1;
        metadata.GridType = geocalls::GridType_t::CARTESIAN;
        metadata.ScalarType = geocalls::ScalarType_t::FLOATING_POINT;
        metadata.ScalarLength = sizeof(float);

        metadata.NumberOfFrames = 1;

        metadata.MinValues = new double[1];
        metadata.MinValues[0] = this->minValue;
        metadata.MaxValues = new double[1];
        metadata.MaxValues[0] = this->maxValue;

        auto const voxelSideLength = this->voxelSizeSlot.Param<core::param::FloatParam>()->Value();

        metadata.Extents[0] = bbox.Width();
        metadata.Extents[1] = bbox.Height();
        metadata.Extents[2] = bbox.Depth();

        metadata.Origin[0] = bbox.Left();
        metadata.Origin[1] = bbox.Bottom();
        metadata.Origin[2] = bbox.Back();

        metadata.Resolution[0] = static_cast<size_t>(std::ceil(metadata.Extents[0] / voxelSideLength));
        metadata.Resolution[1] = static_cast<size_t>(std::ceil(metadata.Extents[1] / voxelSideLength));
        metadata.Resolution[2] = static_cast<size_t>(std::ceil(metadata.Extents[2] / voxelSideLength));

        metadata.SliceDists[0] = new float[0];
        metadata.SliceDists[0][0] = metadata.Extents[0] / static_cast<float>(metadata.Resolution[0] - 1);
        metadata.SliceDists[1] = new float[0];
        metadata.SliceDists[1][0] = metadata.Extents[1] / static_cast<float>(metadata.Resolution[1] - 1);
        metadata.SliceDists[2] = new float[0];
        metadata.SliceDists[2][0] = metadata.Extents[2] / static_cast<float>(metadata.Resolution[2] - 1);
        
        metadata.IsUniform[0] = true;
        metadata.IsUniform[1] = true;
        metadata.IsUniform[2] = true;

        outVolumetricDataCall->SetMetadata(&metadata);

        // Set the data
        outVolumetricDataCall->SetFrameID(this->time);
        outVolumetricDataCall->SetData(this->volume.data());
        outVolumetricDataCall->SetDataHash(this->dataHash);
    }

    return true;
}

bool trialvolume::ParticleToVolume::getExtentCallback(core::Call& caller) {
    auto* inMultiParticleDataCall = this->inParticleDataSlot.CallAs<geocalls::MultiParticleDataCall>();
    if (inMultiParticleDataCall == nullptr) {
        return false;
    }

    auto* outVolumetricDataCall = dynamic_cast<geocalls::VolumetricDataCall*>(&caller);

    auto frameID = outVolumetricDataCall != nullptr ? outVolumetricDataCall->FrameID()
                 : inMultiParticleDataCall != nullptr ? inMultiParticleDataCall->FrameID()
                 : 0;
    inMultiParticleDataCall->SetFrameID(frameID, true);
    if (!(*inMultiParticleDataCall)(1)) {
        // What is this ? Why pass function <1>?
        megamol::core::utility::log::Log::DefaultLog.WriteError(
            "ParticlesToDensity: could not get current frame extents (%u)", time - 1);
        return false;
    }

    if (outVolumetricDataCall != nullptr) {
        outVolumetricDataCall->AccessBoundingBoxes().SetObjectSpaceBBox(inMultiParticleDataCall->GetBoundingBoxes().ObjectSpaceBBox());
        outVolumetricDataCall->AccessBoundingBoxes().SetObjectSpaceClipBox(inMultiParticleDataCall->GetBoundingBoxes().ObjectSpaceClipBox());
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

    auto const voxelSideLength = this->voxelSizeSlot.Param<core::param::FloatParam>()->Value();

    auto const xCells = static_cast<size_t>(std::ceil(bbox.Width() / voxelSideLength));
    auto const yCells = static_cast<size_t>(std::ceil(bbox.Height() / voxelSideLength));
    auto const zCells = static_cast<size_t>(std::ceil(bbox.Depth() / voxelSideLength));

    auto const kernelRadius = this->kernelRadiusSlot.Param<core::param::FloatParam>()->Value();
    auto const kernelCellSpan = static_cast<int>(std::ceil(kernelRadius / voxelSideLength));

    std::function<float(float, float, float)> lengthFunction;
    switch (this->kernelMetricSlot.Param<core::param::EnumParam>()->Value())
    {
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
    switch (this->kernelTypeSlot.Param<core::param::EnumParam>()->Value())
    {
    default:
    case trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST:
        kernel = [](float const dist) -> float {
            return 0.0f;
        };
        break;
    case trialvolume::ParticleToVolume::KERNEL_TYPE_BUMP:
        kernel = [kernelRadius](float const dist) -> float {
            return dist <= kernelRadius ? std::exp(-1.0f / (1.0f - std::pow(dist / kernelRadius, 2.0f))) : 0.0f;
        };
        break;
    }

    std::function<float(float)> applyBoundary;
    switch (this->kernelBoundarySlot.Param<core::param::EnumParam>()->Value())
    {
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLAMP:
        applyBoundary = [](float const value) -> float {
            return std::max(0.0f, std::min(value, 1.0f));
        };
        break;
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_WRAP:
        applyBoundary = [](float const value) -> float {
            return value - std::floor(value);
        };
        break;
    case trialvolume::ParticleToVolume::KERNEL_BOUNDARY_CLIP:
        applyBoundary = [](float const value) -> float {
            return value;
        };
        break;
    }

    this->volume.resize(xCells * yCells * zCells);
    std::fill(this->volume.begin(), this->volume.end(), 0.0f);

    for (size_t i = 0; i < caller->GetParticleListCount(); i++) {
        auto& particleList = caller->AccessParticles(i);
        auto& ps = particleList.GetParticleStore();

        totalParticles += particleList.GetCount();

        auto xAcc = ps.GetXAcc();
        auto yAcc = ps.GetYAcc();
        auto zAcc = ps.GetZAcc();

        for (size_t j = 0; j < particleList.GetCount(); j++) {
            auto x = xAcc->Get_f(j);
            auto y = yAcc->Get_f(j);
            auto z = zAcc->Get_f(j);

            auto xNorm = (x - bbox.Left()) / bbox.Width();
            auto yNorm = (y - bbox.Bottom()) / bbox.Height();
            auto zNorm = (z - bbox.Back()) / bbox.Depth();

            auto xCell = static_cast<int>(std::round(xNorm * xCells));
            auto yCell = static_cast<int>(std::round(yNorm * yCells));
            auto zCell = static_cast<int>(std::round(zNorm * zCells));

            auto isKernel = this->kernelTypeSlot.Param<core::param::EnumParam>()->Value() != trialvolume::ParticleToVolume::KERNEL_TYPE_NEAREST;
            if (!isKernel) {
                auto index = (zCell * yCells + yCell) * xCells + xCell;
                this->volume[index] += 1.0f;
            } else {
                for (auto dz = -kernelCellSpan; dz <= kernelCellSpan; ++dz) {
                    for (auto dy = -kernelCellSpan; dy <= kernelCellSpan; ++dy) {
                        for (auto dx = -kernelCellSpan; dx <= kernelCellSpan; ++dx) {
                            auto const xBounded = std::round(applyBoundary(xNorm + static_cast<float>(dx) / xCells) * xCells);
                            auto const yBounded = std::round(applyBoundary(yNorm + static_cast<float>(dy) / yCells) * yCells);
                            auto const zBounded = std::round(applyBoundary(zNorm + static_cast<float>(dz) / zCells) * zCells);
                            
                            // Check if we are inside the volume
                            if (xBounded < 0 || xBounded >= xCells ||
                                yBounded < 0 || yBounded >= yCells ||
                                zBounded < 0 || zBounded >= zCells) {
                                continue;
                            }

                            auto const index = (static_cast<size_t>(zBounded) * yCells + static_cast<size_t>(yBounded)) * xCells + static_cast<size_t>(xBounded);
                            // FIXME use offset to cell vertex
                            auto const dist = lengthFunction(dx/voxelSideLength, dy/voxelSideLength, dz/voxelSideLength);
                            auto const weight = kernel(dist);

                            this->volume[index] += weight;
                        }
                    }
                }
            }
        }
    }

    auto [min, max] = std::minmax_element(this->volume.begin(), this->volume.end());
    this->minValue = *min;
    this->maxValue = *max;

    const auto endTime = std::chrono::high_resolution_clock::now();
    std::chrono::duration<float, std::milli> diffMillis = endTime - startTime;
    megamol::core::utility::log::Log::DefaultLog.WriteInfo(
        "ParticleToVolume: creation of %u x %u x %u volume from %llu particles took %f ms.", xCells, yCells, zCells,
        totalParticles, diffMillis.count());
    return true;
}
