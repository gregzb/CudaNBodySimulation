#include "nbody_simulation.hpp"

#include <x86intrin.h>

#include <exception>
#include <algorithm>
#include <numeric>
#include <bitset>
#include <tuple>

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

    float max_dim = std::max({rx-lx, ry-ly, rz-lz});
    auto adjust = [max_dim](float &l, float &r) {
        float diff = (max_dim - (r-l)) / 2;
        l -= diff;
        r += diff;
    };
    adjust(lx, rx);
    adjust(ly, ry);
    adjust(lz, rz);

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
    const float epsilon = 0.00001f;

    potential_energy = 0;

    accelerations.assign(bodies.size(), glm::vec3());
    for (unsigned i = 0; i < bodies.size(); i++) {
        for (unsigned j = 0; j < bodies.size(); j++) {
            if (i == j) continue;
            glm::vec3 r(bodies[j].pos-bodies[i].pos);
            glm::vec3 accel(glm::normalize(r)*G*bodies[j].mass/(glm::dot(r, r) + epsilon));
            accelerations[i] += accel;
            if (j > i) {
                potential_energy += (-G * bodies[i].mass * bodies[j].mass) / glm::length(r);
            }
        }
    }
}

void nbody_simulation::add_layer(std::vector<node> &tree, const std::vector<uint64_t> &keys) {
    constexpr uint64_t all_ones_last = 0b1111111111111111111111111111111111111111111111111111111111111110ull;
    constexpr uint64_t three_mask = 0b111ull;

    int stage = tree.size();
    int sec_last = stage-1;

    std::partial_sum(tree[sec_last].children_index.begin(), tree[sec_last].children_index.end(), tree[sec_last].children_index.begin());

    tree.emplace_back(makeNode<node, std::vector>(tree.back().children_index.back() * 8));
    for (unsigned i = 0; i < tree[sec_last].start_index.size(); i++) {
        int first_idx = tree[sec_last].start_index[i];
        int last_idx = tree[sec_last].end_index[i];
        if (first_idx == -1) continue; // wont have any leaves

        uint64_t first_value = keys[first_idx];
        uint64_t base_mask = all_ones_last << (63 - sec_last*3);
        uint64_t base_value = first_value & base_mask;
        int shift_amt = 61 - sec_last*3;

        int prev_idx = -1;

        for (int j = 0; j < 8; j++) {
            uint64_t extension = static_cast<uint64_t>(j) << shift_amt;
            uint64_t extension_mask = three_mask << shift_amt;

            uint64_t search_value = base_value | extension;
            uint64_t search_mask = base_mask | extension_mask;

            auto found = std::lower_bound(keys.begin()+first_idx, keys.begin()+last_idx+1, search_value);
            int found_idx = found-keys.begin();
            if (found == keys.begin()+last_idx+1) {
//                exit(1);
                if (prev_idx != -1) {
                    int end_idx = found_idx-1;
                    int curr_stage_node_idx = (tree[sec_last].children_index[i]-1) * 8 + j;
                    tree[stage].start_index[curr_stage_node_idx-1] = prev_idx;
                    tree[stage].end_index[curr_stage_node_idx-1] = end_idx;
                    tree[stage].children_index[curr_stage_node_idx-1] = 1;
                }
                prev_idx = -1;
                continue;
            }

            uint64_t masked_find = (*found) & search_mask;
            if (prev_idx != -1) {
                int end_idx = found_idx-1;
                int curr_stage_node_idx = (tree[sec_last].children_index[i]-1) * 8 + j;
                tree[stage].start_index[curr_stage_node_idx-1] = prev_idx;
                tree[stage].end_index[curr_stage_node_idx-1] = end_idx;
                tree[stage].children_index[curr_stage_node_idx-1] = 1;
            }
            prev_idx = masked_find != search_value ? -1 : found_idx;
        }

        if (prev_idx != -1) {
            int curr_stage_node_idx = (tree[sec_last].children_index[i]-1) * 8 + 7;
            tree[stage].start_index[curr_stage_node_idx] = prev_idx;
            tree[stage].end_index[curr_stage_node_idx] = last_idx;
            tree[stage].children_index[curr_stage_node_idx] = 1;
        }
    }
}

glm::vec3 nbody_simulation::calc_acceleration(const std::vector<node> &tree, const std::vector<body> &sorted_bodies, const std::vector<com_mass> &com_mass_, int focus_idx, const rect3d &curr_rect, int layer, int node_idx) {
    const float G = 6.67430f*std::pow(10.0f, -11.0f);
    const float epsilon = 0.000001f;

    // inclusive on both sides
    auto query_com = [&com_mass_](int l, int r) {
        com_mass tmp = com_mass_[l];
        tmp.mass = -tmp.mass;
        return combine_com(com_mass_[r+1], tmp);
    };

    int start = tree[layer].start_index[node_idx];
    int end = tree[layer].end_index[node_idx];

    const body &curr_body = sorted_bodies[focus_idx];

    com_mass current_com_mass = query_com(start, end);

    float s2 = curr_rect.sx * curr_rect.sx;
    glm::vec3 diff = curr_body.pos - glm::vec3(current_com_mass.com);
    float d2 = glm::dot(diff, diff) + epsilon;

    bool inside = start <= focus_idx && focus_idx <= end;

    if (end-start+1 == (int)inside) return {0, 0, 0};

    if (s2/d2 < barnes_hut_factor*barnes_hut_factor) {
        com_mass other_com_mass = current_com_mass;
        if (inside) {
            com_mass tmp{curr_body.pos, -curr_body.mass};
            other_com_mass = combine_com(other_com_mass, tmp);
        }

        glm::vec3 r(glm::vec3(other_com_mass.com)-curr_body.pos);
        glm::vec3 accel(glm::normalize(r)*G*(float)other_com_mass.mass/(glm::dot(r, r) + epsilon));
        return accel;
    } else {
        glm::vec3 total{0, 0, 0};

        if (layer == tree_depth) {
            for (int i = tree[layer].start_index[node_idx]; i <= tree[layer].end_index[node_idx]; i++) {
                if (i == focus_idx) continue;
                glm::vec3 r(sorted_bodies[i].pos-curr_body.pos);
                glm::vec3 accel(glm::normalize(r)*G*sorted_bodies[i].mass/(glm::dot(r, r) + epsilon ));
                total += accel;
            }
            return total;
        }

        float nsx = curr_rect.sx / 2;
        float nsy = curr_rect.sy / 2;
        float nsz = curr_rect.sz / 2;

        float cx = curr_rect.lx + nsx;
        float cy = curr_rect.ly + nsy;
        float cz = curr_rect.lz + nsz;

        int children_base_idx = (tree[layer].children_index[node_idx]-1)*8;

        for (int rel_child = 0; rel_child<8; rel_child++) {
            int child_idx = children_base_idx + rel_child;

            bool x_greater = rel_child & 0b100;
            bool y_greater = rel_child & 0b010;
            bool z_greater = rel_child & 0b001;
            float nlx = x_greater ? cx : curr_rect.lx;
            float nly = y_greater ? cy : curr_rect.ly;
            float nlz = z_greater ? cz : curr_rect.lz;
            rect3d tmp_rect{nlx, nly, nlz, nsx, nsy, nsz};

            if (tree[layer+1].start_index[child_idx] != -1) {
                total += calc_acceleration(tree, sorted_bodies, com_mass_, focus_idx, tmp_rect, layer+1, child_idx);
            }
        }
        return total;
    }
}

void nbody_simulation::barnes_hut_cpu_calculate_accelerations() {
    rect3d root_rect {get_bounding_rect()};

    int num_bodies = bodies.size();

    if (num_bodies < 2) {
        return;
    }

    using key_body = std::tuple<uint64_t, body, int>;
    std::vector<key_body> bodies_with_keys;
    for (int i = 0; i < num_bodies; i++) {
        const auto &constrained_pos = convert_xyz(root_rect, bodies[i].pos);
        auto key = get_key(constrained_pos.x, constrained_pos.y, constrained_pos.z);
        bodies_with_keys.emplace_back(key, bodies[i], i);
    }

    std::sort(bodies_with_keys.begin(), bodies_with_keys.end(), [](const key_body &a, const key_body &b){
        return std::get<0>(a) < std::get<0>(b);
    });

    std::vector<com_mass> com_mass_(num_bodies+1, com_mass{{0, 0, 0}, 0});

    std::vector<uint64_t> keys(num_bodies);
    std::vector<body> sorted_bodies(num_bodies);
    for (int i = 0; i < num_bodies; i++) {
        keys[i] = std::get<0>(bodies_with_keys[i]);
        sorted_bodies[i] = std::get<1>(bodies_with_keys[i]);
        com_mass_[i+1] = {sorted_bodies[i].pos, sorted_bodies[i].mass};
    }

    std::partial_sum(com_mass_.begin(), com_mass_.end(), com_mass_.begin(), combine_com);


    std::vector<node> tree;
    tree.emplace_back(makeNode<node, std::vector>(1));
    tree[0].start_index[0] = 0;
    tree[0].end_index[0] = num_bodies-1;
    tree[0].children_index[0] = 1;

    for (int i = 0; i < tree_depth; i++) {
        add_layer(tree, keys);
    }

    accelerations.assign(bodies.size(), glm::vec3());
    for (int i = 0; i < num_bodies; i++) {
        int orig_idx = std::get<2>(bodies_with_keys[i]);
        accelerations[orig_idx] = calc_acceleration(tree, sorted_bodies, com_mass_, i, root_rect);
    }
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