#pragma once

#include <iostream>
#include <vector>

#include <glm/glm.hpp>
#include <glm/gtc/matrix_transform.hpp>
#include <glm/gtc/type_ptr.hpp>

#include "utils.hpp"
#include "body.hpp"

class nbody_simulation {
    private:
    std::vector<body> bodies;
    long long time_steps = 0;
    float time_scale = 1;

    std::vector<glm::vec3> accelerations;

    void calculcate_accelerations();

    public:
    nbody_simulation(float time_scale_ = 1.0f);
    nbody_simulation(const std::vector<body> &bodies_, float time_scale_ = 1.0f);

    void init();
    void step();
    void add_body(const body &body_);
    
    inline void set_time_scale(float time_scale_) {
        time_scale = time_scale_;
    }

    inline float get_time_scale() {
        return time_scale;
    }

    inline long long get_time_steps() {
        return time_steps;
    }
};