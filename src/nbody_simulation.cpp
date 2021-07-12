#include "nbody_simulation.hpp"

nbody_simulation::nbody_simulation(float time_scale_) : time_scale(time_scale_) {

}

nbody_simulation::nbody_simulation(const std::vector<body> &bodies_, float time_scale_ ) : bodies(bodies_), time_scale(time_scale_) {

}

void nbody_simulation::calculcate_accelerations() {
    accelerations.assign(bodies.size(), glm::vec3());
}

void nbody_simulation::init() {

}

void nbody_simulation::step() {
    
    time_steps++;
}

void nbody_simulation::add_body(const body &body_) {
    bodies.push_back(body_);
}