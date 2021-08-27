#pragma once

#include <glad/glad.h>
#include <vector>
#include <unordered_map>
#include <string>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "shader.hpp"

class shader_program
{
private:
    GLuint shader_program_id;
    std::unordered_map<std::string, GLint> uniforms;

public:
    shader_program(const std::vector<GLuint> &shader_ids) : shader_program_id(glCreateProgram())
    {
        for (auto shader_id : shader_ids)
        {
            glAttachShader(shader_program_id, shader_id);
        }
        glLinkProgram(shader_program_id);

        int success;
        char infoLog[512];
        glGetProgramiv(shader_program_id, GL_LINK_STATUS, &success);
        if (!success)
        {
            glGetProgramInfoLog(shader_program_id, 512, NULL, infoLog);
            std::cout << "ERROR::SHADER_PROGRAM::LINKING_FAILED\n"
                      << infoLog << std::endl;
        }
    }

    inline void use() const
    {
        glUseProgram(shader_program_id);
    }

    inline GLuint get_id() const {
        return shader_program_id;
    }

    inline void add_uniform(const std::string& uniform_name) {
        uniforms[uniform_name] = glGetUniformLocation(get_id(), uniform_name.c_str());
    }

    inline GLint get_uniform_location(const std::string& uniform_name) const {
        return uniforms.at(uniform_name);
    }

    inline void set_mvp(const glm::mat4 &model, const glm::mat4 &view, const glm::mat4 &projection) {
        glUniformMatrix4fv(get_uniform_location("model"), 1, GL_FALSE, &model[0][0]);
        glUniformMatrix4fv(get_uniform_location("view"), 1, GL_FALSE, &view[0][0]);
        glUniformMatrix4fv(get_uniform_location("projection"), 1, GL_FALSE, &projection[0][0]);
    }
};