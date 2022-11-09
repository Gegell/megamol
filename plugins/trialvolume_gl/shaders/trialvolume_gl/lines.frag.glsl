#version 400

out layout(location = 0) vec4 frag_color;

in VS_OUT {
    vec4 color;
    float total_mass;
} fs_in;

uniform float min_mass;

void main() {
    if (fs_in.total_mass < min_mass) {
        discard;
    }
    frag_color = fs_in.color;
}
