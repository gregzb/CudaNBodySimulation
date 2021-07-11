#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <array>

#pragma pack(push, 4)
class vertex {
private:
    std::array<glm::float32, 3> data;
public:
    vertex() : data({0, 0, 0}) {

    }

    vertex(const glm::vec3 &pos) : data({pos.x, pos.y, pos.z}) {

    }

    void set_pos(const glm::vec3 &pos) {
        data[0] = pos.x;
        data[1] = pos.y;
        data[2] = pos.z;
    }

    glm::vec3 get_pos() {
        return {data[0], data[1], data[2]};
    }
};
#pragma pack(pop)