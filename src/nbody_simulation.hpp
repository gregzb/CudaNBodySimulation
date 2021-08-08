#pragma once

#include <iostream>
#include <vector>

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

    void naive_cpu_calculcate_accelerations();
    void barnes_hut_cpu_calculcate_accelerations();
    void naive_gpu_calculcate_accelerations();
    void barnes_hut_gpu_calculcate_accelerations();
    void calculcate_accelerations();

public:

    nbody_simulation(float time_scale_ = 1.0f);
    nbody_simulation(const std::vector<body> &bodies_, float time_scale_ = 1.0f);

    void init();
    void step();
    void add_body(const body &body_);

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