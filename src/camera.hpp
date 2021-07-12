#pragma once

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

class camera {
    private:

    float fov;
    glm::vec3 pos;
    glm::vec3 target;

    public:
    camera(float fov_);
    camera(float fov_, const glm::vec3 &pos_);
    camera(float fov_, const glm::vec3 &pos_, const glm::vec3 &target_);

    glm::mat4 get_view_matrix();
    glm::mat4 get_projection_matrix();

    void set_pos(const glm::vec3 &pos_);
    glm::vec3 get_pos();

    void set_target(const glm::vec3 &target_);
    glm::vec3 get_target();

    void set_fov(float fov_);
    float get_fov();
};