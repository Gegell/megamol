#version 430
#extension GL_NV_shader_atomic_float : enable

struct Particle {
    vec3 position;
    vec3 velocity;
};

layout(std430, binding = 0) buffer InputBuffer {
    Particle particles[];
};

struct VoxelData {
#ifdef GL_NV_shader_atomic_float
    vec3 velocity;
    float density;
#else
    uvec3 velocity;
    uint density;
#endif
};

layout(std430, binding = 2) buffer OutputBuffer {
    VoxelData voxels[];
};

struct KernelInfo {
    uint type;
    uint boundary;
    uint metric;
    float radius;
};

uniform KernelInfo kernel;

uniform uvec3 numCells;
uniform vec3 bboxMin;
uniform vec3 bboxMax;
uniform uint numParticles;

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
    if (kernel.metric == 0u) {
        // Euclidean
        return length(v);
    } else if (kernel.metric == 1u) {
        // Manhattan
        return abs(v.x) + abs(v.y) + abs(v.z);
    } else if (kernel.metric == 2u) {
        // Chebyshev
        return max(max(abs(v.x), abs(v.y)), abs(v.z));
    } else {
        return 0.0;
    }
}

vec3 applyBoundary(vec3 pos) {
    switch (kernel.boundary) {
        default:
        case 0u:
            // Clip
            return pos;
        case 1u:
            // Clamp
            return clamp(pos, vec3(0.0), vec3(1.0));
        case 2u:
            // Wrap
            return fract(pos);
    }
}

float computeKernel(float dist) {
    switch (kernel.type) {
        default:
        case 0u:
            // Nearest Neighbor
            return 1.0;
        case 1u:
            // Bump
            return bumpKernel(dist, 1.0/kernel.radius);
        case 2u:
            // Gaussian
            return gaussianKernel(dist, 1.0/kernel.radius);
    }
}

void main() {
    uint index = gl_GlobalInvocationID.x;

    if (index >= numParticles) {
        return;
    }

    vec3 pos = particles[index].position;
    vec3 vel = particles[index].velocity;

    vec3 cellSize = (bboxMax - bboxMin) / vec3(numCells);

    vec3 unitPos = (pos - bboxMin) / (bboxMax - bboxMin);
    vec3 cellPos = fract(unitPos * vec3(numCells) + 0.5);

    ivec3 kernelSpan = ivec3(ceil(kernel.radius / cellSize));
    if (kernel.type == 0u) {
        // Nearest neighbor
        kernelSpan = ivec3(0);
    }

    float totalWeight = 0.0;
    ivec3 cellOffset;
    for (cellOffset.z = -kernelSpan.z; cellOffset.z <= kernelSpan.z; cellOffset.z++) {
        for (cellOffset.y = -kernelSpan.y; cellOffset.y <= kernelSpan.y; cellOffset.y++) {
            for (cellOffset.x = -kernelSpan.x; cellOffset.x <= kernelSpan.x; cellOffset.x++) {
                ivec3 bounded = ivec3(applyBoundary(unitPos + vec3(cellOffset) / vec3(numCells)) * vec3(numCells));
                if (any(notEqual(bounded, clamp(bounded, ivec3(0), ivec3(numCells) - 1)))) {
                    continue;
                }

                float dist = metricLength((cellOffset + cellPos) * cellSize);
                float weight = computeKernel(dist);
                uint cellIndex = (bounded.z * numCells.y + bounded.y) * numCells.x + bounded.x;

                atomicAddFloat(voxels[cellIndex].density, weight);
                atomicAddFloat(voxels[cellIndex].velocity.x, weight * vel.x);
                atomicAddFloat(voxels[cellIndex].velocity.y, weight * vel.y);
                atomicAddFloat(voxels[cellIndex].velocity.z, weight * vel.z);
            }
        }
    }
}
