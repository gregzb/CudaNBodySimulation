#include "nbody_simulation.hpp"

#include <algorithm>

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>
#include <thrust/binary_search.h>
#include <thrust/functional.h>
#include <thrust/device_ptr.h>
#include <thrust/iterator/zip_iterator.h>
#include <thrust/memory.h>
#include <thrust/device_new.h>

#include <chrono>

__global__ void calculate(glm::vec3 *accelerations, body *bodies, int num_bodies) {
    const float G = 6.67430f * std::pow(10.0f, -11);
    const float epsilon = 0.0000001f;

    uint32_t body_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (body_idx >= num_bodies) return;
    for (unsigned j = 0; j < num_bodies; j++) {
        if (j == body_idx) continue;
        glm::vec3 r(bodies[j].pos - bodies[body_idx].pos);
        glm::vec3 accel(glm::normalize(r) * G * bodies[j].mass / (glm::dot(r, r) + epsilon));
        accelerations[body_idx] += accel;
    }
}

void nbody_simulation::naive_gpu_calculate_accelerations() {
    int num_bodies = bodies.size();

    glm::vec3 *gpu_accelerations;
    cudaMallocManaged(&gpu_accelerations, num_bodies * sizeof(glm::vec3));
    std::memset(gpu_accelerations, 0, num_bodies * sizeof(glm::vec3));

    body *gpu_bodies;
    cudaMallocManaged(&gpu_bodies, num_bodies * sizeof(body));
    for (int i = 0; i < num_bodies; i++) {
        gpu_bodies[i] = bodies[i];
    }

    dim3 dim_block(256, 1, 1);
    uint32_t grid_size = (num_bodies + dim_block.x - 1) / dim_block.x;
    dim3 dim_grid(grid_size, 1, 1);

    calculate<<<dim_grid, dim_block>>>(gpu_accelerations, gpu_bodies, num_bodies);
    cudaDeviceSynchronize();

    accelerations.assign(num_bodies, glm::vec3());
    for (int i = 0; i < num_bodies; i++) {
        accelerations[i] = gpu_accelerations[i];
    }

    cudaFree(gpu_accelerations);
    cudaFree(gpu_bodies);
}

struct device_node {
    int *start_index;
    int *end_index;
    int *children_index;
    int n;
};

struct device_extract_com {
    __device__ nbody_simulation::com_mass operator()(const body &body_) const {
        return {body_.pos, body_.mass};
    }
};

struct device_combine_com {
    __device__ nbody_simulation::com_mass
    operator()(const nbody_simulation::com_mass &a, const nbody_simulation::com_mass &b) {
        return {(a.com * a.mass + b.com * b.mass) / (a.mass + b.mass), a.mass + b.mass};
    }
};

__device__ void make_device_node(device_node &d_node, int n) {
    d_node.start_index = new int[n];
    memset(d_node.start_index, -1, n * sizeof(int));

    d_node.end_index = new int[n];
    memset(d_node.end_index, -1, n * sizeof(int));

    d_node.children_index = new int[n];
    memset(d_node.children_index, 0, n * sizeof(int));

    d_node.n = n;
}

__host__ void initialize_device_node_from_host(device_node &d_node, int n) {
    cudaMalloc(&d_node.start_index, n * sizeof(int));
    cudaMemset(d_node.start_index, 0xFF, n * sizeof(int)); // this will set all the ints to -1

    cudaMalloc(&d_node.end_index, n * sizeof(int));
    cudaMemset(d_node.end_index, 0xFF, n * sizeof(int)); // this will set all the ints to -1

    cudaMalloc(&d_node.children_index, n * sizeof(int));
    cudaMemset(d_node.children_index, 0, n * sizeof(int)); // this will set all the ints to 0

    d_node.n = n;
}

__host__ void free_device_node_from_host(device_node &d_node) {
    cudaFree(d_node.start_index);
    cudaFree(d_node.end_index);
    cudaFree(d_node.children_index);
}

__device__ void free_device_node(device_node &d_node) {
    delete[] d_node.start_index;
    delete[] d_node.end_index;
    delete[] d_node.children_index;
}

struct device_initialize_tree {
    template<typename Tuple>
    __device__ void operator()(Tuple t) {
        device_node &layer = thrust::get<0>(t);
        int num_bodies = thrust::get<1>(t);
        make_device_node(layer, 1);
        layer.start_index[0] = 0;
        layer.end_index[0] = num_bodies - 1;
        layer.children_index[0] = 1;
    }
};

struct device_initialize_layer {
    template<typename Tuple>
    __device__ void operator()(Tuple t) {
        device_node &layer = thrust::get<0>(t);
        int new_num_indices = thrust::get<1>(t);
        make_device_node(layer, new_num_indices * 8);
    }
};

__device__ void print_binary(uint64_t val, char *out) {
    uint64_t mask = 0b1;
    mask <<= 63;
    for (int i = 0; i < 64; i++) {
        out[i] = (val & mask) ? '1' : '0';
        mask >>= 1;
    }
}

struct device_compute_layer {
    template<class Tuple>
    __device__ void operator()(Tuple t) {
        device_node *tree = thrust::get<0>(t);
        uint64_t *keys = thrust::get<1>(t);
        int data = thrust::get<2>(t);
        int stage = thrust::get<3>(t);
        int sec_last = stage - 1;

        constexpr uint64_t all_ones = 0b1111111111111111111111111111111111111111111111111111111111111111ull;
        constexpr uint64_t all_ones_last = 0b1111111111111111111111111111111111111111111111111111111111111110ull;
        constexpr uint64_t three_mask = 0b111ull;

        int idx = data >> 3;
        int extension_idx = data & 0b111;

        int first_idx = tree[sec_last].start_index[idx];
        int last_idx = tree[sec_last].end_index[idx];
        if (first_idx == -1) return;

        uint64_t first_value = keys[first_idx];
        uint64_t base_mask = all_ones_last << (63 - sec_last * 3);
        uint64_t base_value = first_value & base_mask;
        int shift_amt = 61 - sec_last * 3;

        uint64_t extension = static_cast<uint64_t>(extension_idx) << shift_amt;
        uint64_t extension_mask = three_mask << shift_amt;

        uint64_t search_value = base_value | extension;
        uint64_t search_mask = base_mask | extension_mask;

        uint64_t extension_upper = all_ones >> ((sec_last + 1) * 3);
        uint64_t search_value_upper = search_value | extension_upper;

        auto found_start = thrust::lower_bound(thrust::seq, keys + first_idx, keys + last_idx + 1, search_value);
        int found_start_idx = found_start - keys;

        if (found_start == keys + last_idx + 1) {
            return;
        };
        if ((((*found_start) & search_mask) != search_value)) {
            return;
        }

        auto found_end = thrust::upper_bound(thrust::seq, keys + first_idx, keys + last_idx + 1, search_value_upper);
        int found_end_idx = found_end - keys;

        int curr_stage_node_idx = (tree[sec_last].children_index[idx] - 1) * 8 + extension_idx;

        tree[stage].start_index[curr_stage_node_idx] = found_start_idx;
        tree[stage].end_index[curr_stage_node_idx] = found_end_idx - 1;
        tree[stage].children_index[curr_stage_node_idx] = 1;
    }
};

struct device_calc_acceleration {
    template<class Tuple>
    __device__ void operator()(Tuple t) {
        device_node *tree = thrust::get<0>(t);
        body *bodies = thrust::get<1>(t);
        nbody_simulation::com_mass *com_mass_ = thrust::get<2>(t);
        glm::vec3 *accelerations = thrust::get<3>(t);
        int *idxes = thrust::get<4>(t);
        int focus_idx = thrust::get<5>(t);
        float rect_size = thrust::get<6>(t);
        float barnes_hut_factor = thrust::get<7>(t);
        int tree_depth = thrust::get<8>(t);

        const float G = 6.67430f * std::pow(10.0f, -11.0f);
        const float epsilon = 0.01f;

        const body &curr_body = bodies[focus_idx];
        int write_idx = idxes[focus_idx];

        auto com_combiner = device_combine_com();

        auto query_com = [&com_mass_, &com_combiner](int l, int r) {
            nbody_simulation::com_mass tmp = com_mass_[l];
            tmp.mass = -tmp.mass;
            return com_combiner(com_mass_[r + 1], tmp);
        };

        struct stack_item {
            float size;
            int layer;
            int node_idx;
        };
        stack_item stack[64];
        int stack_size = 0;

        stack[stack_size] = {rect_size, 0, 0};
        stack_size++;

        glm::vec3 tmp_total_acceleration{0, 0, 0};
        while (stack_size > 0) {
            stack_size--;
            auto item = stack[stack_size];
            int start = tree[item.layer].start_index[item.node_idx];
            int end = tree[item.layer].end_index[item.node_idx];

            nbody_simulation::com_mass current_com_mass = query_com(start, end);
            float s2 = item.size * item.size;
            glm::vec3 diff = curr_body.pos - glm::vec3(current_com_mass.com);
            float d2 = glm::dot(diff, diff) + epsilon;

            bool inside = start <= focus_idx && focus_idx <= end;
            if (end - start + 1 == (int) inside) continue;

            if (s2 / d2 < barnes_hut_factor) {
                nbody_simulation::com_mass other_com_mass = current_com_mass;
                if (inside) {
                    nbody_simulation::com_mass tmp{curr_body.pos, -curr_body.mass};
                    other_com_mass = com_combiner(other_com_mass, tmp);
                }

                glm::vec3 r(glm::vec3(other_com_mass.com) - curr_body.pos);
                glm::vec3 accel(glm::normalize(r) * G * (float) other_com_mass.mass / (glm::dot(r, r) + epsilon));
                tmp_total_acceleration += accel;
            } else {
                if (item.layer == tree_depth) {
                    for (int i = tree[item.layer].start_index[item.node_idx];
                         i <= tree[item.layer].end_index[item.node_idx]; i++) {
                        if (i == focus_idx) continue;
                        glm::vec3 r(bodies[i].pos - curr_body.pos);
                        glm::vec3 accel(glm::normalize(r) * G * bodies[i].mass / (glm::dot(r, r) + epsilon));
                        tmp_total_acceleration += accel;
                    }
                } else {
                    int children_base_idx = (tree[item.layer].children_index[item.node_idx] - 1) * 8;

                    for (int rel_child = 0; rel_child < 8; rel_child++) {
                        int child_idx = children_base_idx + rel_child;

                        if (tree[item.layer + 1].start_index[child_idx] != -1) {
                            stack_item new_item{item.size / 2, item.layer + 1, child_idx};
                            stack[stack_size] = new_item;
                            stack_size++;
                        }
                    }
                }

            }
        }

        accelerations[write_idx] = tmp_total_acceleration;
    }
};

struct device_functor_free_node {
    __device__ void operator()(device_node &d_node) {
        free_device_node(d_node);
    }
};

void nbody_simulation::barnes_hut_gpu_calculate_accelerations() {
    rect3d root_rect{get_bounding_rect()};

    int num_bodies = bodies.size();

    if (num_bodies < 2) {
        return;
    }

    std::vector<uint64_t> keys(num_bodies);
    std::vector<int> idxes(num_bodies);
    for (int i = 0; i < num_bodies; i++) {
        const auto &constrained_pos = convert_xyz(root_rect, bodies[i].pos);
        keys[i] = get_key(constrained_pos.x, constrained_pos.y, constrained_pos.z);
        idxes[i] = i;
    }

    thrust::device_vector<uint64_t> device_keys = keys;
    thrust::device_vector<uint64_t> device_keys_copy = device_keys;

    thrust::device_vector<body> device_bodies = bodies;
    thrust::sort_by_key(device_keys.begin(), device_keys.end(), device_bodies.begin());

    thrust::device_vector<int> device_idxes = idxes;
    thrust::sort_by_key(device_keys_copy.begin(), device_keys_copy.end(), device_idxes.begin());

    thrust::device_vector<com_mass> com_mass_(num_bodies + 1, com_mass{{0, 0, 0}, 0});
    thrust::transform(device_bodies.begin(), device_bodies.end(), com_mass_.begin() + 1, device_extract_com());
    thrust::inclusive_scan(com_mass_.begin(), com_mass_.end(), com_mass_.begin(), device_combine_com());

    thrust::device_vector<device_node> tree(tree_depth + 1);
    thrust::constant_iterator<int> num_bodies_iterator(num_bodies - 1);

    device_node host_layer = tree[0];
    initialize_device_node_from_host(host_layer, 1);
    int *init_mem = (int *) malloc(sizeof(int) * 3);
    init_mem[0] = 0;
    init_mem[1] = num_bodies - 1;
    init_mem[2] = 1;
    cudaMemcpy(host_layer.start_index, init_mem + 0, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(host_layer.end_index, init_mem + 1, sizeof(int), cudaMemcpyHostToDevice);
    cudaMemcpy(host_layer.children_index, init_mem + 2, sizeof(int), cudaMemcpyHostToDevice);
    free(init_mem);
    tree[0] = host_layer;

    device_node *tree_ptr = thrust::raw_pointer_cast(tree.data());
    auto tree_iterator = thrust::make_constant_iterator(tree_ptr);

    for (int i = 0; i < tree_depth; i++) {
        int stage = i + 1;
        int sec_last = i;

        int prev_n = tree[sec_last].operator device_node().n;
        int *children_index = tree[sec_last].operator device_node().children_index;

        auto device_children_index = thrust::device_pointer_cast(children_index);

        thrust::inclusive_scan(device_children_index, device_children_index + prev_n, device_children_index);

        int *last_children_value_ptr = (int *) malloc(sizeof(int));
        cudaMemcpy(last_children_value_ptr, children_index + prev_n - 1, sizeof(int), cudaMemcpyDeviceToHost);
        int last_children_value = *last_children_value_ptr;
        free(last_children_value_ptr);

        device_node host_curr_layer = tree[stage];
        initialize_device_node_from_host(host_curr_layer, last_children_value * 8);
        tree[stage] = host_curr_layer;

        uint64_t *keys_ptr = thrust::raw_pointer_cast(device_keys.data());
        auto keys_iterator = thrust::make_constant_iterator(keys_ptr);
        auto idx = thrust::make_counting_iterator(0);
        auto stage_iterator = thrust::make_constant_iterator(stage);
        auto zipped_start = thrust::make_zip_iterator(
                thrust::make_tuple(tree_iterator, keys_iterator, idx, stage_iterator));
        auto zipped_end = thrust::make_zip_iterator(
                thrust::make_tuple(tree_iterator + prev_n * 8, keys_iterator + prev_n * 8, idx + prev_n * 8,
                                   stage_iterator + prev_n * 8));

        thrust::for_each(thrust::device, zipped_start, zipped_end, device_compute_layer());
    }

    body *bodies_ptr = thrust::raw_pointer_cast(device_bodies.data());
    auto bodies_iterator = thrust::make_constant_iterator(bodies_ptr);

    com_mass *com_masses = thrust::raw_pointer_cast(com_mass_.data());
    auto com_mass_iterator = thrust::make_constant_iterator(com_masses);

    accelerations.assign(bodies.size(), glm::vec3());
    thrust::device_vector<glm::vec3> device_accelerations(num_bodies, {0, 0, 0});

    glm::vec3 *accelerations_ptr = thrust::raw_pointer_cast(device_accelerations.data());
    auto accelerations_iterator = thrust::make_constant_iterator(accelerations_ptr);

    int *idxes_ptr = thrust::raw_pointer_cast(device_idxes.data());
    auto idxes_iterator = thrust::make_constant_iterator(idxes_ptr);

    auto focus_idx_iterator = thrust::make_counting_iterator(0);

    auto size_iterator = thrust::make_constant_iterator(root_rect.sx);

    auto barnes_hut_iterator = thrust::make_constant_iterator(barnes_hut_factor);

    auto tree_depth_iterator = thrust::make_constant_iterator(tree_depth);

    auto start_calculation_iter = thrust::make_zip_iterator(
            thrust::make_tuple(tree_iterator, bodies_iterator, com_mass_iterator, accelerations_iterator,
                               idxes_iterator, focus_idx_iterator, size_iterator, barnes_hut_iterator,
                               tree_depth_iterator));
    auto end_calculation_iter = thrust::make_zip_iterator(
            thrust::make_tuple(tree_iterator + num_bodies, bodies_iterator + num_bodies, com_mass_iterator + num_bodies,
                               accelerations_iterator + num_bodies, idxes_iterator + num_bodies,
                               focus_idx_iterator + num_bodies, size_iterator + num_bodies,
                               barnes_hut_iterator + num_bodies, tree_depth_iterator + num_bodies));
    thrust::for_each(thrust::device, start_calculation_iter, end_calculation_iter, device_calc_acceleration());

    thrust::host_vector<glm::vec3> host_accelerations;
    host_accelerations = device_accelerations;
    for (int i = 0; i < num_bodies; i++) {
        accelerations[i] = host_accelerations[i];
    }

    for (int i = 0; i < tree_depth + 1; i++) {
        device_node layer = tree[i];
        free_device_node_from_host(layer);
    }
}