uniform vec4 viewAttr;

uniform vec3 camIn;
uniform vec3 camUp;
uniform vec3 camRight;

uniform mat4 MVinv;
uniform mat4 MVP;
uniform mat4 MVPinv;
uniform mat4 MVPtransp;
uniform mat4 NormalM;

varying vec4 objPos;
varying vec4 camPos;
varying vec4 lightPos;

uniform vec2 planes;