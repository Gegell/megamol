#version 400

out layout(location = 0) vec4 frag_color;

uniform vec4 line_color;

void main() {
    frag_color = line_color;
}
