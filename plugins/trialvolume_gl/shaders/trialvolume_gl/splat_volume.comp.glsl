#version 430
#extension GL_NV_shader_atomic_float : enable


layout(std430, binding = 0) buffer inPosition {
    vec3 partPos[];
};
layout(std430, binding = 1) buffer inVelocity {
    vec3 partVel[];
};

#ifdef GL_NV_shader_atomic_float
    layout(std430, binding = 2) buffer outDensity {
        float splatDensity[];
    };
    layout(std430, binding = 3) buffer outVelocity {
        vec3 splatVelocity[];
    };
#else
    layout(std430, binding = 2) buffer outDensity {
        uint splatDensity[];
    };
    layout(std430, binding = 3) buffer outVelocity {
        uvec3 splatVelocity[];
    };
#endif

uniform KernelData {
    uint type;
    uint boundary;
    uint metric;
    float radius;
} kernel;

uniform Dimension {
    uvec3 numCells;
    vec3 bboxMin;
    vec3 bboxMax;
    uint numParticles;
};

layout(local_size_x = 32, local_size_y = 1, local_size_z = 1) in;

#ifdef GL_NV_shader_atomic_float
    #define atomicAddFloat atomicAdd
#else
    #define atomicAddFloat(MEM, DATA) {                                       \
        uint expected_mem = floatBitsToUint(MEM);                             \
        uint result = floatBitsToUint(uintBitsToFloat(expected_mem) + DATA);  \
        uint prev_value = atomicCompSwap(MEM, expected_mem, result);          \
        while (prev_value != expected_mem) {                                  \
            expected_mem = prev_value;                                        \
            result = floatBitsToUint(uintBitsToFloat(expected_mem) + DATA);   \
            prev_value = atomicCompSwap(MEM, expected_mem, result);           \
        }                                                                     \
    }
#endif

float bumpKernel(float r, float eps) {
    return r*eps < 1.0 ? exp(-1.0/(1.0-(r*eps)*(r*eps))) : 0.0;
}

float gaussianKernel(float r, float eps) {
    return exp(-r*r/(2.0*eps*eps));
}

float metricLength(vec3 v) {
    if (kernel.metric == 0) {
        // Euclidean
        return length(v);
    } else if (kernel.metric == 1) {
        // Manhattan
        return abs(v.x) + abs(v.y) + abs(v.z);
    } else if (kernel.metric == 2) {
        // Chebyshev
        return max(max(abs(v.x), abs(v.y)), abs(v.z));
    } else {
        return 0.0;
    }
}

void main() {
    uint workGroupIndex = (gl_WorkGroupID.z * gl_NumWorkGroups.y + gl_WorkGroupID.y) * gl_NumWorkGroups.x + gl_WorkGroupID.x;
    uint workGroupSize = gl_WorkGroupSize.x * gl_WorkGroupSize.y * gl_WorkGroupSize.z;
    uint index = workGroupIndex * workGroupSize + gl_LocalInvocationID.x;

    if (index >= numParticles) {
        return;
    }

    vec3 pos = partPos[index];
    vec3 vel = partVel[index];

    vec3 cellSize = (bboxMax - bboxMin) / vec3(numCells);

    vec3 cellPos = pos / cellSize;
    ivec3 cellId = ivec3(floor(cellPos));
    cellPos -= vec3(cellId);

    ivec3 kernelSpan = ivec3(ceil(kernel.radius / cellSize));
    if (kernel.type == 0u) {
        // Nearest neighbor
        kernelSpan = ivec3(0);
    }

    float totalWeight = 0.0;
    // 2 Stages
    //   1st stage: compute total weight
    //   2nd stage: compute density and velocity
    for (int stage = 0; stage < 2; stage++) {
        ivec3 cellOffset;
        for (cellOffset.z = -kernelSpan.z; cellOffset.z <= kernelSpan.z; cellOffset.z++) {
            for (cellOffset.y = -kernelSpan.y; cellOffset.y <= kernelSpan.y; cellOffset.y++) {
                for (cellOffset.x = -kernelSpan.x; cellOffset.x <= kernelSpan.x; cellOffset.x++) {
                    vec3 r = (vec3(cellOffset) + 0.5 + cellPos) * cellSize;

                    float weight = 1.0;
                    if (kernel.type == 1) {
                        // Bump
                        weight = bumpKernel(metricLength(r), 1.0/kernel.radius);
                    } else if (kernel.type == 2) {
                        // Gaussian
                        weight = gaussianKernel(metricLength(r), 1.0/kernel.radius);
                    }

                    if (stage == 0) {
                        totalWeight += weight;
                    } else {
                        ivec3 cellIdOffset = cellId + cellOffset;

                        ivec3 clampedcellIdOffset = clamp(cellIdOffset, ivec3(0), ivec3(numCells) - 1);
                        if (kernel.boundary == 0u) {
                            // Clip
                            if (any(notEqual(clampedcellIdOffset, cellIdOffset))) {
                                continue;
                            }
                        } else if (kernel.boundary == 1u) {
                            // Clamp
                            cellIdOffset = clampedcellIdOffset;
                        } else if (kernel.boundary == 2u) {
                            // Wrap
                            cellIdOffset = cellIdOffset % ivec3(numCells);
                        }

                        uint cellIdOffsetFlat = (cellIdOffset.z * numCells.y + cellIdOffset.y) * numCells.x + cellIdOffset.x;

                        weight /= totalWeight;
                        #if 1
                            atomicAddFloat(splatDensity[cellIdOffsetFlat], weight);
                            atomicAddFloat(splatVelocity[cellIdOffsetFlat].x, weight * vel.x);
                            atomicAddFloat(splatVelocity[cellIdOffsetFlat].y, weight * vel.y);
                            atomicAddFloat(splatVelocity[cellIdOffsetFlat].z, weight * vel.z);
                        #else
                            splatDensity[cellIdOffsetFlat] = floatBitsToUint(uintBitsToFloat(splatDensity[cellIdOffsetFlat]) + weight);
                            splatVelocity[cellIdOffsetFlat]=floatBitsToUint(uintBitsToFloat(splatVelocity[cellIdOffsetFlat]) + weight * vel);
                        #endif
                    }
                }
            }
        }
    }
    splatDensity[index] = floatBitsToUint(1.0);
    splatVelocity[index] = floatBitsToUint(vel);
}
