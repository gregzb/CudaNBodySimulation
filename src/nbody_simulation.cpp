#include "nbody_simulation.hpp"

#include <x86intrin.h>

#include <exception>
#include <algorithm>

nbody_simulation::nbody_simulation(float time_scale_) : time_scale(time_scale_) {

}

nbody_simulation::nbody_simulation(const std::vector<body> &bodies_, float time_scale_ ) : bodies(bodies_), time_scale(time_scale_), accelerations(bodies_.size()){

}

nbody_simulation::rect3d nbody_simulation::get_bounding_rect() const{
    float lx, ly, lz, rx, ry, rz;
    lx = std::numeric_limits<float>::max();
    ly = std::numeric_limits<float>::max();
    lz = std::numeric_limits<float>::max();
    rx = -std::numeric_limits<float>::max();
    ry = -std::numeric_limits<float>::max();
    rz = -std::numeric_limits<float>::max();
    float eps = std::numeric_limits<float>::epsilon() * 10;
    for (auto & body : bodies) {
        lx = std::min(lx, body.pos.x);
        ly = std::min(ly, body.pos.y);
        lz = std::min(lz, body.pos.z);
        rx = std::max(rx, body.pos.x);
        ry = std::max(ry, body.pos.y);
        rz = std::max(rz, body.pos.z);
    }

    return rect3d{lx-eps, ly-eps, lz-eps, rx-lx+2*eps, ry-ly+2*eps, rz-lz+2*eps};
}

uint64_t nbody_simulation::get_key (float x, float y, float z) {
    uint32_t x_int = *(uint32_t*)&x;
    uint32_t y_int = *(uint32_t*)&y;
    uint32_t z_int = *(uint32_t*)&z;

    constexpr uint64_t mask_x = 0b1001001001001001001001001001001001001001001001001001001001001001ull; // 22 ones
    constexpr uint64_t mask_y = 0b0100100100100100100100100100100100100100100100100100100100100100ull; // 21 ones
    constexpr uint64_t mask_z = 0b0010010010010010010010010010010010010010010010010010010010010010ull; // 21 ones

    uint64_t final_value = _pdep_u64(x_int >> (32-9-22), mask_x) | _pdep_u64(y_int >> (32-9-21), mask_y) | _pdep_u64(z_int >> (32-9-21), mask_z);
    return final_value;
};

void nbody_simulation::naive_cpu_calculate_accelerations() {
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

void nbody_simulation::barnes_hut_cpu_calculate_accelerations() {
    rect3d root_rect {get_bounding_rect()};

    int num_bodies = bodies.size();

    using key_body = std::pair<uint64_t, body>;
    std::vector<key_body> bodies_with_keys;
    for (int i = 0; i < num_bodies; i++) {
        const auto &constrained_pos = convert_xyz(root_rect, bodies[i].pos);
        auto key = get_key(constrained_pos.x, constrained_pos.y, constrained_pos.z);
        bodies_with_keys.emplace_back(key, bodies[i]);
    }

    std::sort(bodies_with_keys.begin(), bodies_with_keys.end(), [](const key_body &a, const key_body &b){
        return a.first < b.first;
    });

    std::vector<float> mass(num_bodies);
    std::vector<glm::vec3> com(num_bodies);

    std::vector<uint64_t> keys(num_bodies);
    for (int i = 0; i < num_bodies; i++) {
        keys[i] = bodies_with_keys[i].first;
        mass[i] = bodies_with_keys[i].second.mass;
        com[i] = bodies_with_keys[i].second.pos;
    }



    std::vector<node> tree;
    tree.emplace_back(makeNode<node, std::vector>(1));
    for (const body &body_ : bodies) {
        add_point_mass(tree[0], 0, body_.pos, body_.mass);
    }

    tree.emplace_back(makeNode<node, std::vector>(tree.back().start.size() * 8));

}

void nbody_simulation::calculate_accelerations() {
    switch(calculation_backend) {
        case CalculationBackend::NAIVE_CPU:
            naive_cpu_calculate_accelerations();
            break;
        case CalculationBackend::BARNES_HUT_CPU:
            barnes_hut_cpu_calculate_accelerations();
            break;
        case CalculationBackend::NAIVE_GPU:
            naive_gpu_calculate_accelerations();
            break;
        case CalculationBackend::BARNES_HUT_GPU:
            barnes_hut_gpu_calculate_accelerations();
            break;
    }
}

void nbody_simulation::init() {
    calculate_accelerations();
}

void nbody_simulation::step() {
    float half = 1.0f/2;

    kinetic_energy = 0;
    for (unsigned i = 0; i < bodies.size(); i++) {
        bodies[i].vel += accelerations[i]*time_scale*half;
        bodies[i].pos += bodies[i].vel*time_scale;
    }
    calculate_accelerations();
    for (unsigned i = 0; i < bodies.size(); i++) {
        float mag = glm::length(bodies[i].vel);
        kinetic_energy += mag*mag*bodies[i].mass/2;
        bodies[i].vel += accelerations[i]*time_scale*half;
    }
    time_steps++;
}

void nbody_simulation::add_body(const body &body_) {
    bodies.push_back(body_);
    accelerations.emplace_back();
}