#version 400

layout(lines) in;
layout (line_strip, max_vertices = 2) out;

in VS_OUT {
    vec4 color;
    float total_mass;
} gs_in[];

uniform float min_mass;
out vec4 color;

void main() {
    float mass = min(gs_in[0].total_mass, gs_in[1].total_mass);
    if (mass < min_mass) {
        return;
    }

    gl_Position = gl_in[0].gl_Position;
    color = gs_in[0].color;
    EmitVertex();

    gl_Position = gl_in[1].gl_Position;
    color = gs_in[1].color;
    EmitVertex();
}
