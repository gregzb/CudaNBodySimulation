#include "nbody_simulation.hpp"

#include <thrust/host_vector.h>
#include <thrust/device_vector.h>
#include <thrust/sort.h>

__global__ void calculate(glm::vec3* accelerations, body* bodies, int num_bodies) {
    const float G = 6.67430f*std::pow(10.0f, -11);
    const float epsilon = 0.0000001f;

    int body_idx = blockIdx.x * blockDim.x + threadIdx.x;
    if (body_idx >= num_bodies) return;
    for (unsigned j = 0; j < num_bodies; j++) {
        if (j == body_idx) continue;
        glm::vec3 r(bodies[j].pos-bodies[body_idx].pos);
        glm::vec3 accel(glm::normalize(r)*G*bodies[j].mass/(glm::dot(r, r)+epsilon));
        accelerations[body_idx] += accel;
    }
}

void nbody_simulation::naive_gpu_calculcate_accelerations() {
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
    int grid_size = (num_bodies+dim_block.x-1)/dim_block.x;
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

void nbody_simulation::barnes_hut_gpu_calculcate_accelerations() {
    rect root_rect {get_bounding_rect()};

    std::vector<int> keys(bodies.size());
    for (unsigned i = 0; i < bodies.size(); i++) {
        const auto &constrained_pos = convert_xy(root_rect, bodies[i].pos.x, bodies[i].pos.y);
        keys[i] = get_key(constrained_pos.first, constrained_pos.second);
    }

    thrust::device_vector<int> device_keys = keys;
    thrust::device_vector<body> device_bodies = bodies;

    thrust::sort_by_key(device_keys.begin(), device_keys.begin() + device_keys.size(), device_bodies.begin());

    thrust::host_vector<int> sorted_keys = device_keys;
    for (auto val : sorted_keys) {
        std::cout << val << std::endl;
    }
}