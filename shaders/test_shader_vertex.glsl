#version 460 core

layout (location = 0) in vec3 pos;

out vec4 vertexColor;

uniform mat4 mvp;

void main()
{
    gl_Position = mvp * vec4(pos.x, pos.y, pos.z, 1.0);
    vertexColor = vec4(0.2, 0.7, 0.9, 1.0);
}