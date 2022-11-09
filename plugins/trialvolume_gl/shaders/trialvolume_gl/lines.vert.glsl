#version 400

#define COLORMODE_VELOCITY 0
#define COLORMODE_TOTAL_MASS 1
#define COLORMODE_LOCAL_ID 2
#define COLORMODE_FRAME 3

layout(location = 0) in vec3 bbox_min;
layout(location = 1) in vec3 bbox_max;
layout(location = 2) in vec3 center_of_mass;
layout(location = 3) in vec3 velocity;
layout(location = 4) in float total_mass;
layout(location = 5) in uint frame;
layout(location = 6) in uint frame_local_id;

uniform mat4 mvp;
uniform int color_mode;

uniform float max_mass;
uniform float max_frame;
uniform float max_frame_local_id;

out VS_OUT {
    vec4 color;
    float total_mass;
} vs_out;

void main() {
    vec4 pos = vec4(center_of_mass.xyz, 1.0);
    gl_Position = mvp * pos;
    vs_out.total_mass = total_mass;
    vs_out.color = vec4(1.0);
    switch (color_mode) {
        default:
        case COLORMODE_VELOCITY:
            vs_out.color.rgb = normalize(velocity) * .5 + .5;
            break;
        case COLORMODE_TOTAL_MASS:
            vs_out.color.rgb = vec3(total_mass) / max_mass;
            break;
        case COLORMODE_LOCAL_ID:
            vs_out.color.rgb = vec3(frame_local_id) / max_frame_local_id;
            break;
        case COLORMODE_FRAME:
            vs_out.color.rgb = vec3(frame) / max_frame;
            break;
    }
}
