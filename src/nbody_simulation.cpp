#include "nbody_simulation.hpp"

#include <x86intrin.h>

#include <exception>

nbody_simulation::nbody_simulation(float time_scale_) : time_scale(time_scale_) {

}

nbody_simulation::nbody_simulation(const std::vector<body> &bodies_, float time_scale_ ) : bodies(bodies_), time_scale(time_scale_) {

}

nbody_simulation::rect nbody_simulation::get_bounding_rect() const{
    float lx, ly, rx, ry;
    lx = std::numeric_limits<float>::max();
    ly = std::numeric_limits<float>::max();
    rx = -std::numeric_limits<float>::max();
    ry = -std::numeric_limits<float>::max();
    float eps = std::numeric_limits<float>::epsilon() * 10;
    for (auto & body : bodies) {
        lx = std::min(lx, body.pos.x);
        ly = std::min(ly, body.pos.y);
        rx = std::max(rx, body.pos.x);
        ry = std::max(ry, body.pos.y);
    }

    return rect{lx-eps, ly-eps, rx-lx+eps, ry-ly+eps};
}

uint32_t nbody_simulation::get_key (float x, float y) {
    uint32_t x_int = *(uint32_t*)&x;
    uint32_t y_int = *(uint32_t*)&y;

    constexpr uint32_t mask_x = 0b01010101010101010101010101010101u;
    constexpr uint32_t mask_y = 0b10101010101010101010101010101010u;

    uint32_t final_value = _pdep_u32(x_int >> (32-9-16), mask_x) | _pdep_u32(y_int >> (32-9-16), mask_y);
    return final_value;
};

void nbody_simulation::naive_cpu_calculcate_accelerations() {
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

void nbody_simulation::barnes_hut_cpu_calculcate_accelerations() {
    rect root_rect {get_bounding_rect()};

    std::vector<int> keys(bodies.size());
    for (unsigned i = 0; i < bodies.size(); i++) {
        const auto &constrained_pos = convert_xy(root_rect, bodies[i].pos.x, bodies[i].pos.y);
        keys[i] = get_key(constrained_pos.first, constrained_pos.second);
    }

    // thrust::device_vector<int> device_keys = keys;
    // thrust::device_vector<body> device_bodies = bodies;

    // thrust::sort_by_key(device_keys.begin(), device_keys.begin() + device_keys.size(), device_bodies.begin());
}

void nbody_simulation::calculcate_accelerations() {
    switch(calculation_backend) {
        case CalculationBackend::NAIVE_CPU:
            naive_cpu_calculcate_accelerations();
            break;
        case CalculationBackend::BARNES_HUT_CPU:
            barnes_hut_cpu_calculcate_accelerations();
            break;
        case CalculationBackend::NAIVE_GPU:
            naive_gpu_calculcate_accelerations();
            break;
        case CalculationBackend::BARNES_HUT_GPU:
            barnes_hut_gpu_calculcate_accelerations();
            break;
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