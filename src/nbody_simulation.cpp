#include "nbody_simulation.hpp"

nbody_simulation::nbody_simulation(float time_scale_) : time_scale(time_scale_) {

}

nbody_simulation::nbody_simulation(const std::vector<body> &bodies_, float time_scale_ ) : bodies(bodies_), time_scale(time_scale_) {

}

void nbody_simulation::calculcate_accelerations() {
    const float G = 6.67430f*std::pow(10.0f, -11);
    const float epsilon = 0.0000001f;

    potential_energy = 0;

    accelerations.assign(bodies.size(), glm::vec3());
    for (unsigned i = 0; i < bodies.size(); i++) {
        for (unsigned j = 0; j < bodies.size(); j++) {
            if (i == j) continue;
            glm::vec3 r(bodies[j].pos-bodies[i].pos);
            glm::vec3 accel(glm::normalize(r)*G*bodies[j].mass/(glm::dot(r, r)+epsilon));
            accelerations[i] += accel;
            if (j > i) {
                potential_energy += (-G * bodies[i].mass * bodies[j].mass) / glm::length(r);
            }
        }
    }
}

void nbody_simulation::init() {
    calculcate_accelerations();
}

void nbody_simulation::step() {
    float half = 1.0f/2;

    kinetic_energy = 0;
    for (unsigned i = 0; i < bodies.size(); i++) {
        bodies[i].vel += accelerations[i]*time_scale*half;
        bodies[i].pos += bodies[i].vel*time_scale;
    }
    calculcate_accelerations();
    for (unsigned i = 0; i < bodies.size(); i++) {
        float mag = glm::length(bodies[i].vel);
        kinetic_energy += mag*mag*bodies[i].mass/2;
        bodies[i].vel += accelerations[i]*time_scale*half;
    }
    time_steps++;
}

void nbody_simulation::add_body(const body &body_) {
    bodies.push_back(body_);
}