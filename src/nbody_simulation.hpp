#pragma once

#include <iostream>
#include <vector>
#include <array>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "utils.hpp"
#include "body.hpp"

class nbody_simulation
{
public:
    enum class CalculationBackend {
        NAIVE_CPU = 0,
        BARNES_HUT_CPU = 1,
        NAIVE_GPU = 2,
        BARNES_HUT_GPU = 3
    };

private:
    std::vector<body> bodies;
    long long time_steps = 0;
    float time_scale = 1;

    float potential_energy = 0;
    float kinetic_energy = 0;

    std::vector<glm::vec3> accelerations;

    CalculationBackend calculation_backend = CalculationBackend::NAIVE_CPU;

    void naive_cpu_calculate_accelerations();
    void barnes_hut_cpu_calculate_accelerations();
    void naive_gpu_calculate_accelerations();
    void barnes_hut_gpu_calculate_accelerations();
    void calculate_accelerations();

    struct rect3d
    {
        float lx, ly, lz, sx, sy, sz;
    };

    struct node{
        std::vector<int> start_index;
        std::vector<std::array<int, 8>> children;
    };

    template <class NODE, template <class> class CONTAINER>
    static inline NODE makeNode(int items) {
        return NODE { CONTAINER<int>(items, -1), CONTAINER<std::array<int, 8>>(items, {-1, -1, -1, -1, -1, -1, -1, -1})};
    }

    template <class NODE>
    static inline void add_point_mass(NODE &node_, int idx, const glm::vec3 &com, float mass) {
        node_.com[idx] = (node_.com[idx] * node_.mass[idx] + com * mass) / (node_.mass[idx] + mass);
        node_.mass[idx] += mass;
    }

    template <class NODE>
    static inline glm::vec3 new_com(NODE &node_, int idx, const glm::vec3 &com, float mass) {
        return (node_[idx].com * node_[idx].mass + com * mass) / (node_[idx].mass + mass);
    }

public:

    nbody_simulation(float time_scale_ = 1.0f);
    nbody_simulation(const std::vector<body> &bodies_, float time_scale_ = 1.0f);

    void init();
    void step();
    void add_body(const body &body_);

    rect3d get_bounding_rect() const;

    inline glm::vec3 convert_xyz(const rect3d &root_rect, const glm::vec3 &pos) {
        //between 1 and 2
        return glm::vec3{(pos.x-root_rect.lx)/root_rect.sx+1, (pos.y-root_rect.ly)/root_rect.sy+1, (pos.z-root_rect.lz)/root_rect.sz+1};
    }

    uint64_t get_key (float x, float y, float z);

    inline void set_time_scale(float time_scale_)
    {
        time_scale = time_scale_;
    }

    inline float get_time_scale()
    {
        return time_scale;
    }

    inline long long get_time_steps()
    {
        return time_steps;
    }

    inline float get_energy() {
        return potential_energy + kinetic_energy;
    }

    inline const std::vector<body>& get_bodies() {
        return bodies;
    }

    inline CalculationBackend get_backend() {
        return calculation_backend;
    }

    inline void set_backend(CalculationBackend state) {
        calculation_backend = state;
    }
};