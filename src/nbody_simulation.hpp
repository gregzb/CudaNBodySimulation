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

    struct rect3d
            {
        float lx, ly, lz, sx, sy, sz;
            };

    struct node{
        std::vector<int> start_index;
        std::vector<int> end_index;
        std::vector<int> children_index;
    };

    struct com_mass {
        glm::dvec3 com;
        double mass;
    };

private:
    std::vector<body> bodies;
    long long time_steps = 0;
    float time_scale = 1;

    float potential_energy = 0;
    float kinetic_energy = 0;

    std::vector<glm::vec3> accelerations;

    CalculationBackend calculation_backend = CalculationBackend::NAIVE_CPU;

    float barnes_hut_factor = 0.7;
    int tree_depth = 8;

    void naive_cpu_calculate_accelerations();
    void barnes_hut_cpu_calculate_accelerations();
    void naive_gpu_calculate_accelerations();
    void barnes_hut_gpu_calculate_accelerations();
    void calculate_accelerations();

    template <class NODE, template <class> class CONTAINER>
    static inline NODE makeNode(int items) {
        return NODE { CONTAINER<int>(items, -1), CONTAINER<int>(items, -1), CONTAINER<int>(items, 0)};
    }

    void add_layer(std::vector<node> &tree, const std::vector<uint64_t> &keys);
    glm::vec3 calc_acceleration(const std::vector<node> &tree, const std::vector<body> &sorted_bodies, const std::vector<com_mass> &com_mass_, int focus_idx, const rect3d &curr_rect, int layer = 0, int node_idx = 0);

    inline static com_mass combine_com(const com_mass &a, const com_mass &b) {
        return {(a.com * a.mass + b.com * b.mass) / (a.mass + b.mass), a.mass + b.mass};
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

    inline float get_time_scale() const
    {
        return time_scale;
    }

    inline long long get_time_steps() const
    {
        return time_steps;
    }

    inline float get_energy() const{
        return potential_energy + kinetic_energy;
    }

    inline const std::vector<body>& get_bodies() {
        return bodies;
    }

    inline CalculationBackend get_backend() const{
        return calculation_backend;
    }

    inline void set_backend(CalculationBackend state) {
        calculation_backend = state;
    }

    inline float get_barnes_hut_factor() const {
        return barnes_hut_factor;
    }

    inline void set_barnes_hut_factor(float barnes_hut_factor_) {
        barnes_hut_factor = barnes_hut_factor_;
    }
};