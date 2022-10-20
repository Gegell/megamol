
#include "trialvolume/ParticleToVolume.h"

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
        : BaseParticleToVolume()
        , splatting_method_slot_("SplattingMethod", "The splatting method to use") {
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

bool trialvolume::ParticleToVolume::computeVolume(geocalls::MultiParticleDataCall* caller) {
    std::fill(density_.begin(), density_.end(), 0.0f);
    std::fill(velocity_.begin(), velocity_.end(), 0.0f);

    switch (splatting_method_slot_.Param<core::param::EnumParam>()->Value()) {
    case trialvolume::ParticleToVolume::SPLAT_METHOD_KERNEL:
        return computeKernel(caller);
    case trialvolume::ParticleToVolume::SPLAT_METHOD_NATURAL_NEIGHBOR:
        return computeNaturalNeighborhood(caller);
    default:
        return false;
    }
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
