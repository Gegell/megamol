#version 400

layout(location = 0) in vec3 bbox_min;
layout(location = 1) in vec3 bbox_max;
layout(location = 2) in vec3 center_of_mass;
layout(location = 3) in vec3 velocity;
layout(location = 4) in float total_mass;
layout(location = 5) in uint frame;
layout(location = 6) in uint frame_local_id;

uniform mat4 mvp;

out VS_OUT {
    vec3 bbox_min;
    float total_mass;
    vec3 bbox_max;
    uint frame;
} vs_out;

void main() {
    vec4 pos = vec4(center_of_mass.xyz, 1.0);
    gl_Position = mvp * pos;
    vs_out.bbox_min = bbox_min;
    vs_out.bbox_max = bbox_max;
    vs_out.frame = frame;
    vs_out.total_mass = total_mass;
}
