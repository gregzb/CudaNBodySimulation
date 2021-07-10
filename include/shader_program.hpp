#pragma once

#include <glad/glad.h>
#include <vector>

#include "shader.hpp"

class shader_program
{
private:
    GLuint shader_program_id;

public:
    shader_program(std::vector<GLuint> shader_ids) : shader_program_id(glCreateProgram())
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

    void use()
    {
        glUseProgram(shader_program_id);
    }

    GLuint get_id() {
        return shader_program_id;
    }
};