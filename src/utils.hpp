#pragma once

#include <iostream>
#include <fstream>
#include <string>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

std::ostream& operator<<(std::ostream &stream, const glm::vec3 &vec) {
    stream <<  "{" << vec.x << " " << vec.y << " " << vec.z << "}";
    return stream;
}

std::string fileToString(std::string file_path)
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