#version 400

out layout(location = 0) vec4 frag_color;

in vec4 color;

void main() {
    frag_color = color;
}
