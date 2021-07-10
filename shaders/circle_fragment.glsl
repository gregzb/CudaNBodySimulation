#version 460 core

out vec4 FragColor;

uniform float inv_width_multiplier;

in vec4 vertexColor;
in vec2 center;
in vec2 curr;
in float radius_out;

void main()
{
    float dist = inv_width_multiplier*inv_width_multiplier*(curr.x-center.x)*(curr.x-center.x) + (curr.y-center.y)*(curr.y-center.y);
    if (dist > radius_out*radius_out) {
        discard;
    }
    FragColor = vertexColor;
} 