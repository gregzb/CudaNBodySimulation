#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "utils.hpp"

struct body {
    glm::vec3 pos;
    glm::vec3 vel;
    float mass;
};