#version 460 core

layout (location = 0) in vec2 pos;
layout (location = 1) in vec2 offset;
layout (location = 2) in float radius;
layout (location = 3) in vec2 velocity;
layout (location = 4) in float mass;

out vec4 vertexColor;
out vec2 center;
out vec2 curr;
out float radius_out;

void main()
{
    gl_Position = vec4(pos.x+offset.x, pos.y+offset.y, 1.0, 1.0);

    float vel_mag = dot(velocity, velocity);
    vel_mag = clamp(vel_mag*10000, 0, 1);

    vec4 slow_color = vec4(0.2, 0.2, 0.9, 1.0);
    vec4 fast_color = vec4(0.9, 0.2, 0.2, 1.0);

    vertexColor = mix(slow_color, fast_color, vel_mag);
    
    radius_out = radius;
    center = offset;
    curr = vec2(pos.x+offset.x, pos.y+offset.y);
}