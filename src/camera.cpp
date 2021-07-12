#include "camera.hpp"
#include "graphics_settings.hpp"

camera::camera(float fov_) : fov(fov_) {

}

camera::camera(float fov_, const glm::vec3 &pos_) : fov(fov_), pos(pos_) {

}

camera::camera(float fov_, const glm::vec3 &pos_, const glm::vec3 &target_) : fov(fov_), pos(pos_), target(target_) {

}

glm::mat4 camera::get_view_matrix() {
    glm::mat4 view = glm::lookAt(pos, target, {0, 1, 0}); // make sure to initialize matrix to identity matrix first
    return view;
}

glm::mat4 camera::get_projection_matrix() {
    glm::mat4 projection = glm::perspective(fov, (float)SCR_WIDTH / (float)SCR_HEIGHT, 0.1f, 100.0f);
    return projection;
}

void camera::set_pos(const glm::vec3 &pos_) {
    pos = pos_;
}

glm::vec3 camera::get_pos() {
    return pos;
}

void camera::set_target(const glm::vec3 &target_) {
    target = target_;
}

glm::vec3 camera::get_target() {
    return target;
}

void camera::set_fov(float fov_) {
    fov = fov_;
}

float camera::get_fov() {
    return fov;
}