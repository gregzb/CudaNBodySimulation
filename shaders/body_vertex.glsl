#version 460 core

layout (location = 0) in vec3 pos;
layout (location = 1) in vec3 normal;
layout (location = 2) in vec3 offset;
layout (location = 3) in vec3 color;

out vec4 vertex_color;
out vec3 normal_out;
out vec3 frag_pos;

uniform mat4 model;
uniform mat4 view;
uniform mat4 projection;

uniform mat3 normal_model_view;

void main()
{
    vec4 model_transformed = model * vec4(pos, 1.0f) + vec4(offset, 0.0f);
    vec4 view_transformed = view * model_transformed;
    gl_Position = projection * view_transformed;
    vertex_color = vec4(color, 1.0f);
    normal_out = normal_model_view*normal;
    frag_pos = vec3(view_transformed);
}