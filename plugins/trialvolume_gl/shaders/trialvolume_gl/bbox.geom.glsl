#version 400

layout(lines) in;
layout (line_strip, max_vertices = 24) out;

in VS_OUT {
    vec3 bbox_min;
    float total_mass;
    vec3 bbox_max;
    uint frame;
} gs_in[];

uniform mat4 mvp;
uniform float time;
uniform float min_mass;
out vec4 color;

const vec3 box_mults[8] = vec3[8](
    vec3(-1.0, -1.0, -1.0),
    vec3(-1.0, -1.0,  1.0),
    vec3(-1.0,  1.0, -1.0),
    vec3(-1.0,  1.0,  1.0),
    vec3( 1.0, -1.0, -1.0),
    vec3( 1.0, -1.0,  1.0),
    vec3( 1.0,  1.0, -1.0),
    vec3( 1.0,  1.0,  1.0)
);

const uint box_indices[24] = uint[](
    0, 1, 2, 3, 4, 5, 6, 7,
    0, 2, 1, 3, 4, 6, 5, 7,
    0, 4, 1, 5, 2, 6, 3, 7
);

void generate_cube_lines(vec3 box_min, vec3 box_max) {
    vec3 box_size = box_max - box_min;
    vec3 box_center = box_min + box_size / 2.0;
    vec3 box_half_size = box_size / 2.0;


    for (int i = 0; i < 24; i += 2) {
        for (int j = 0; j < 2; j++) {
            vec3 vert_mult = box_mults[box_indices[i + j]];
            vec3 vert_pos = box_center + vert_mult * box_half_size;
            gl_Position = mvp * vec4(vert_pos, 1.0);
            EmitVertex();
        }
        EndPrimitive();
    }
}

void main() {
    if (gs_in[0].frame != uint(floor(time))) {
        return;
    }
    float fTime = fract(time);

    float mass = mix(gs_in[0].total_mass, gs_in[1].total_mass, fTime);
    if (mass < min_mass) {
        return;
    }

    vec3 bbox_min = mix(gs_in[0].bbox_min, gs_in[1].bbox_min, fTime);
    vec3 bbox_max = mix(gs_in[0].bbox_max, gs_in[1].bbox_max, fTime);

    generate_cube_lines(bbox_min, bbox_max);
}
