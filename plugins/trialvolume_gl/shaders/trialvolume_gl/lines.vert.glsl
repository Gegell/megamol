#version 400

layout(location = 0) in vec3 bbox_min;
layout(location = 1) in vec3 bbox_max;
layout(location = 2) in vec3 center_of_mass;
layout(location = 3) in vec3 velocity;
layout(location = 4) in float total_mass;
layout(location = 5) in uint frame;
layout(location = 6) in uint frame_local_id;

uniform mat4 mvp;
out vec4 color;

void main() {
    vec4 pos = vec4(center_of_mass.xyz, 1.0);
    color = vec4(normalize(velocity), 1.0);
    gl_Position = mvp * pos;
}
