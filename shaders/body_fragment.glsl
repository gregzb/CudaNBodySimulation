#version 460 core

out vec4 frag_color;

uniform vec3 light_pos;
uniform vec3 view_pos;

in vec4 vertex_color;
in vec3 frag_pos;
in vec3 normal_out;

void main()
{
    vec3 norm = normalize(normal_out);

    vec3 light_color = vec3(0.8, 0.8, 0.8);
    vec3 light_dir = normalize(light_pos - frag_pos);

    // frag_color = vec4(norm, 1.0f);

    float ambient_strength = 0.15;
    vec3 ambient = ambient_strength * light_color;

    float angle_of_incidence = max(dot(norm, light_dir), 0.0);
    vec3 diffuse = angle_of_incidence * light_color;

    float specular_strength = 0.8;
    vec3 view_dir = normalize(view_pos - frag_pos);
    vec3 reflect_dir = reflect(-light_dir, norm);
    float spec = pow(max(dot(view_dir, reflect_dir), 0.0), 64);
    vec3 specular = specular_strength * spec * light_color;

    vec3 result = (ambient+diffuse+specular) * vec3(vertex_color);
    frag_color = vec4(result, 1.0);
}