#version 460 core

layout (location = 0) in vec2 pos;
// layout (location = 1) in vec2 offset;
// layout (location = 2) in float radius;

out vec4 vertexColor;
// out vec2 center;
// out vec2 curr;
// out float radius_out;

void main()
{
    // gl_Position = vec4(pos.x+offset.x, pos.y+offset.y, 1.0, 1.0);
    gl_Position = vec4(pos.x, pos.y, 1.0, 1.0);
    vertexColor = vec4(0.2, 0.7, 0.9, 1.0);
    // radius_out = radius;
    // center = offset;
    // curr = vec2(pos.x+offset.x, pos.y+offset.y);
}