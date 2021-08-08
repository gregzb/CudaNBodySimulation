#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include <vector>
#include <array>

#pragma pack(push, 1)
class vertex {
private:
    // pos[3], normal[3]
    std::array<glm::float32, 6> data;
public:
    vertex() : data({0, 0, 0}) {

    }

    vertex(const glm::vec3 &pos) : data({pos.x, pos.y, pos.z, 0, 0, 0}) {

    }

    vertex(const glm::vec3 &pos, const glm::vec3 &normal) : data({pos.x, pos.y, pos.z, normal.x, normal.y, normal.z}) {

    }

    void set_pos(const glm::vec3 &pos) {
        data[0] = pos.x;
        data[1] = pos.y;
        data[2] = pos.z;
    }

    glm::vec3 get_pos() {
        return {data[0], data[1], data[2]};
    }

    void set_normal(const glm::vec3 &normal) {
        data[3] = normal.x;
        data[4] = normal.y;
        data[5] = normal.z;
    }
};
#pragma pack(pop)