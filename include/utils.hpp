#pragma once

#include <iostream>
#include <fstream>
#include <string>

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