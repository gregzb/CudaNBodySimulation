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


//device_node device_add_layer(device_node &tree, const thrust::device_vector<uint64_t> &keys) {
//    constexpr uint64_t all_ones_last = 0b1111111111111111111111111111111111111111111111111111111111111110ull;
//    constexpr uint64_t three_mask = 0b111ull;
//
//    int stage = tree.size();
//    int sec_last = stage-1;
//
//    thrust::inclusive_scan(tree.);
//}

//__device__ void device_add_layer(device_node* tree, uint64_t *keys, int curr_size) {
//    constexpr uint64_t all_ones = 0b1111111111111111111111111111111111111111111111111111111111111111ull;
//    constexpr uint64_t all_ones_last = 0b1111111111111111111111111111111111111111111111111111111111111110ull;
//    constexpr uint64_t three_mask = 0b111ull;
//
//    int stage = curr_size;
//    int sec_last = stage - 1;
//
//    int prev_num_indices = tree[sec_last].n;
//
//    thrust::inclusive_scan(thrust::device, tree[sec_last].children_index, tree[sec_last].children_index+prev_num_indices, tree[sec_last].children_index);
//    int new_num_indices = tree[sec_last].children_index[prev_num_indices-1]*8;
//    make_device_node(tree[stage], new_num_indices);
//    thrust::counting_iterator<int> start_idx(0);
//    thrust::counting_iterator<int> end_idx = start_idx + prev_num_indices*8;
//
//    // MUST BE COPY!
//    thrust::for_each(thrust::device, start_idx, end_idx, [=](int data) {
//        int idx = data >> 3;
//        int extension_idx = data & 0b111;
//
//        int first_idx = tree[sec_last].start_index[idx];
//        int last_idx = tree[sec_last].end_index[idx];
//        if (first_idx == -1) return;
//
//        uint64_t first_value = keys[first_idx];
//        uint64_t base_mask = all_ones_last << (63 - sec_last*3);
//        uint64_t base_value = first_value & base_mask;
//        int shift_amt = 61 - sec_last*3;
//
//        uint64_t extension = static_cast<uint64_t>(extension_idx) << shift_amt;
//        uint64_t extension_mask = three_mask << shift_amt;
//
//        uint64_t search_value = base_value | extension;
//        uint64_t search_mask = base_mask | extension_mask;
//
//        uint64_t extension_upper = all_ones >> ((sec_last+1)*3);
//        uint64_t search_value_upper = search_value | extension_upper;
//
//        auto found_start = thrust::lower_bound(thrust::device, keys+first_idx, keys+last_idx+1, search_value);
//        if (found_start == keys+last_idx+1 || (((*found_start) & search_mask) != search_value)) return;
//
//        auto found_end = thrust::upper_bound(thrust::device, keys+first_idx, keys+last_idx+1, search_value_upper);
//
//        int found_start_idx = found_start - keys;
//        int found_end_idx = found_end - keys;
//
//        int curr_stage_node_idx = (tree[sec_last].children_index[idx]-1) * 8 + extension_idx;
//        tree[stage].start_index[curr_stage_node_idx] = found_start_idx;
//        tree[stage].end_index[curr_stage_node_idx] = found_end_idx-1;
//        tree[stage].children_index[curr_stage_node_idx] = 1;
//
////        auto found = std::lower_bound(keys+first_idx, keys+last_idx+1, search_value);
//    });
//}

//__global__ void barnes_hut_gpu_function(glm::vec3* accelerations, body* bodies, uint64_t *keys, int n, int tree_depth, float barnes_hut_factor) {
////    thrust::sort_by_key(thrust::seq, keys, keys + n, bodies);
////    thrust::device_ptr<uint64_t> keys_ptr(keys);
////    thrust::sort(thrust::device, keys_ptr, keys_ptr + n);
////    glm::vec3 tot = {0, 0, 0};
////    uint64_t key_tot = 0;
////    for (int i = 0; i < n; i++) {
////        tot += bodies[i].pos;
////        key_tot += keys[i];
////    }
////    test_func<<<1, 1>>>(bodies, keys, n);
//
//    auto *com_mass_ = new nbody_simulation::com_mass[n+1];
//    nbody_simulation::com_mass tmp{{0, 0, 0}, 0};
//    com_mass_[0] = tmp;
//    thrust::transform(thrust::device, bodies, bodies + n, com_mass_+1, device_extract_com());
//    thrust::inclusive_scan(thrust::device, com_mass_, com_mass_+n+1, com_mass_, device_combine_com());
//
//    auto *tree = new device_node[tree_depth+1];
//
//    make_device_node(tree[0], 1);
//    tree[0].start_index[0] = 0;
//    tree[0].end_index[0] = n-1;
//    tree[0].children_index[0] = 1;
//
//    for (int i = 0; i < tree_depth; i++) {
//        device_add_layer(tree, keys, i+1);
//    }
//
//    for (int i = 0; i < tree_depth; i++) {
//        free_device_node(tree[i]);
//    }
//
//    delete[] tree;
//    delete[] com_mass_;
//}

struct device_node{
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
    __device__ nbody_simulation::com_mass operator()(const nbody_simulation::com_mass &a, const nbody_simulation::com_mass &b) {
        return {(a.com * a.mass + b.com * b.mass) / (a.mass + b.mass), a.mass + b.mass};
    }
};

__device__ void make_device_node(device_node &d_node, int n) {
    d_node.start_index = new int[n];
    memset(d_node.start_index, -1, n*sizeof(int));

    d_node.end_index = new int[n];
    memset(d_node.end_index, -1, n*sizeof(int));

    d_node.children_index = new int[n];
    memset(d_node.children_index, 0, n*sizeof(int));

    d_node.n = n;
}

struct device_initialize_tree {
    template <typename Tuple>
    __device__ void operator() (Tuple t) {
        device_node &layer = thrust::get<0>(t);
        int num_bodies = thrust::get<1>(t);
        make_device_node(layer, 3);
        layer.start_index[0] = 0;
        layer.end_index[0] = num_bodies-1;
        layer.children_index[0] = 1;
        layer.children_index[1] = 0;
        layer.children_index[2] = 1;
    }
};

struct device_initialize_layer {
    template <typename Tuple>
    __device__ void operator() (Tuple t) {
        device_node &layer = thrust::get<0>(t);
        int new_num_indices = thrust::get<1>(t);
        make_device_node(layer, new_num_indices);
    }
};

struct device_children_iterator {
    template <class Tuple>
    __device__ int* operator() (Tuple t) {
        return (thrust::get<0>(t) + thrust::get<1>(t))->children_index;
    }
};

//struct device_prefix_sum_layer {
//    __device__ void operator() (int idx) {
//        device_node &layer = thrust::get<0>(t);
//        int num_bodies = thrust::get<1>(t);
//        make_device_node(layer, 1);
//        layer.start_index[0] = 0;
//        layer.end_index[0] = num_bodies-1;
//        layer.children_index[0] = 1;
//    }
//};

__device__ void free_device_node(device_node &d_node) {
    delete[] d_node.start_index;
    delete[] d_node.end_index;
    delete[] d_node.children_index;
}

//__global__ void debug_kernel(int* dat, int n) {
//    for (int i = 0; i < n; i++) {
//        printf("%d: %d\n", i, dat[i]);
//    }
//}

void nbody_simulation::barnes_hut_gpu_calculate_accelerations() {
    rect3d root_rect {get_bounding_rect()};

//    size_t free, total;
//    printf("\n");
//    cudaMemGetInfo(&free,&total);
//    printf("%d MB free of total %d MB\n",free/1024/1024,total/1024/1024);

    static bool called = false;
    if (!called) {
        cudaDeviceSetLimit(cudaLimitMallocHeapSize, 1024ll*1024ll*1024ll*1ll); //set limit to 1gb
    }
    called = true;

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

    thrust::device_vector<com_mass> com_mass_(num_bodies+1, com_mass{{0, 0, 0}, 0});
    thrust::transform(device_bodies.begin(), device_bodies.end(), com_mass_.begin()+1, device_extract_com());
    thrust::inclusive_scan(com_mass_.begin(), com_mass_.end(), com_mass_.begin(), device_combine_com());

    thrust::device_vector<device_node> tree(tree_depth+1);
    thrust::constant_iterator<int> num_bodies_iterator(num_bodies-1);

    thrust::for_each_n(thrust::make_zip_iterator(thrust::make_tuple(tree.begin(), num_bodies_iterator)),
                     1,
                     device_initialize_tree());

    for (int i = 0; i < tree_depth; i++) {
        int stage = i + 1;
        int sec_last = i;

        int prev_n = tree[sec_last].operator device_node().n;
        int *children_index = tree[sec_last].operator device_node().children_index;

        auto device_children_index = thrust::device_pointer_cast(children_index);

        thrust::inclusive_scan(device_children_index, device_children_index+prev_n, device_children_index);
//        thrust::device_ptr<int> last_children_idx = thrust::device_new<int>(1);
//        thrust::transform(thrust::device, children_index+prev_n-1, children_index+prev_n, last_children_idx, );
        thrust::for_each_n(thrust::make_zip_iterator(thrust::make_tuple(tree.begin()+stage, device_children_index+prev_n-1)),
                           1,
                           device_initialize_layer());

//        auto prev_tree_level_it = tree.begin() + sec_last;
//        auto transformed_iterator = thrust::make_transform_iterator(prev_tree_level, );

//        // combine into one, fix function above, considering extending from unary_function, can this even be done?
//        thrust::constant_iterator<thrust::device_vector<device_node>::iterator> tree_iterator(tree.begin());
//        thrust::constant_iterator<int> sec_last_iterator(sec_last);
//        auto transform_start_iterator = thrust::make_transform_iterator(thrust::make_zip_iterator(thrust::make_tuple(tree_iterator, sec_last_iterator)), device_children_iterator());
//        auto transform_end_iterator = thrust::make_transform_iterator(thrust::make_zip_iterator(thrust::make_tuple(tree_iterator, sec_last_iterator)), device_children_iterator());
//        thrust::inclusive_scan(thrust::device, transform_start_iterator, transform_end_iterator, transform_start_iterator);

//        thrust::inclusive_scan(thrust::device, );
    }

//    thrust::for_each(thrust::make_zip_iterator(thrust::make_tuple(tree.begin(), num_bodies_iterator)), thrust::make_zip_iterator(thrust::make_tuple(tree.begin()+1, num_bodies_iterator+1)), device_initialize_tree);



//    thrust::host_vector<uint64_t> host_keys = device_keys;



//    std::sort(keys.begin(), keys.end());

//    body *device_bodies;
//    cudaMallocHost(&device_bodies, num_bodies * sizeof(body));
//
//    cudaMemcpy(device_bodies, bodies.data(), num_bodies, cudaMemcpyHostToDevice);
//
////    std::cout << cudaGetErrorName(err) << std::endl;
//
//    uint64_t *device_keys;
//    cudaMallocHost(&device_keys, num_bodies * sizeof(uint64_t));
//
//    cudaMemcpy(device_keys, keys.data(), num_bodies, cudaMemcpyHostToDevice);
//
//    glm::vec3 *device_accelerations;
//    cudaMallocHost(&device_accelerations, num_bodies * sizeof(glm::vec3));
//
////    for (int i = 0; i < num_bodies; i++) {
////        device_bodies[i] = bodies[i];
////        device_keys[i] = keys[i];
////        device_accelerations[i] = {0, 0, 0};
////    }
//
////    std::cout << "bruh" << std::endl;
////    thrust::sort_by_key(thrust::device, device_keys, device_keys + num_bodies, device_bodies);
////    std::cout << "eet" << std::endl;
////    std::cout << "here" << std::endl;
//    thrust::device_ptr<uint64_t> device_keys_ptr(device_keys);
//    thrust::device_ptr<body> device_bodies_ptr(device_bodies);
//    thrust::sort_by_key(device_keys_ptr, device_keys_ptr + num_bodies, device_bodies_ptr);
////    std::cout << "oop" << std::endl;
//    barnes_hut_gpu_function<<<1, 1>>>(device_accelerations, device_bodies, device_keys, num_bodies, tree_depth, barnes_hut_factor);
//    cudaDeviceSynchronize();
////    std::cout << "beep" << std::endl;
//
//    auto kernel_err = cudaGetLastError();
//
////    std::cout << "kernel: " << cudaGetErrorName(kernel_err) << std::endl;
//
//    cudaFreeHost(device_bodies);
//    cudaFreeHost(device_keys);
//    cudaFreeHost(device_accelerations);

//    exit(0);

//    thrust::device_vector<uint64_t> device_keys = keys;
//    thrust::device_vector<body> device_bodies = bodies;
//
//    thrust::sort_by_key(device_keys.begin(), device_keys.begin() + num_bodies, device_bodies.begin());
//    thrust::device_vector<com_mass> com_mass_(num_bodies+1, com_mass{{0, 0, 0}, 0});
//    thrust::transform(device_bodies.begin(), device_bodies.end(), com_mass_.begin() + 1, device_extract_com());
//    thrust::inclusive_scan(com_mass_.begin(), com_mass_.end(), com_mass_.begin(), device_combine_com());
//
//    thrust::device_vector<device_node> tree;
//    tree.push_back(device_node{thrust::device_vector<int>(1, -1), thrust::device_vector<int>(1, -1),thrust::device_vector<int>(1, 0)});
//
//    for (int i = 0; i < tree_depth; i++) {
//        tree.push_back(device_add_layer(tree.back().operator device_node(), device_keys));
//    }
}