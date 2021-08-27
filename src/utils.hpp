#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

inline std::ostream& operator<<(std::ostream &stream, const glm::vec3 &vec) {
    stream <<  "{" << vec.x << " " << vec.y << " " << vec.z << "}";
    return stream;
}

inline std::string file_to_string(const std::string &file_path)
{
    std::ifstream file(file_path);
    std::string contents;
    file.seekg(0, std::ios::end);
    contents.resize(file.tellg());
    file.seekg(0, std::ios::beg);
    file.read(&contents[0], contents.size());
    file.close();
    return contents;
}