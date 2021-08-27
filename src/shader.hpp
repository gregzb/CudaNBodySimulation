# pragma once

# include <glad/glad.h>

# include <string>

# include "utils.hpp"

template <GLenum shader_type>
class shader
{
private:
    GLuint shader_id;

public:
    inline shader(const std::string &file_path) {
        std::string contents = file_to_string(file_path);

        shader_id = glCreateShader(shader_type);

        const char* shader_source = contents.c_str();

        glShaderSource(shader_id, 1, &shader_source, NULL);
        glCompileShader(shader_id);

        int  success;
        char infoLog[512];
        glGetShaderiv(shader_id, GL_COMPILE_STATUS, &success);

        if(!success)
        {
            glGetShaderInfoLog(shader_id, 512, NULL, infoLog);
            std::cout << "ERROR::SHADER::" << file_path << "::COMPILATION_FAILED\n" << infoLog << std::endl;
        }
    }

    inline GLuint get_id() const {
        return shader_id;
    }

    inline void release() {
        glDeleteShader(shader_id);
    }
};