#version 460 core

out vec4 frag_color;

uniform vec2 screen_size;

uniform vec2 highlight_center;
uniform vec2 radii;

uniform vec4 highlight_color;

void main()
{
    vec2 center_screen_pos = (highlight_center + 1) / 2 * screen_size;
    vec2 screen_pos = gl_FragCoord.xy;
    float min_dim = min(screen_size.x, screen_size.y);
    vec2 pixel_radii = radii * min_dim;

    float dist = distance(screen_pos, center_screen_pos);
    if (dist < pixel_radii.x || dist > pixel_radii.y) discard;
    
    frag_color = highlight_color;
}