#include <iostream>
#include <math.h>

#include <chrono>
using namespace std::chrono;

// Kernel function to add the elements of two arrays
__global__ void add(int n, float *x, float *y)
{
    int chunk_size = n/(gridDim.x*blockDim.x);
    int idx = (blockIdx.x * blockDim.x + threadIdx.x) * chunk_size;
    for (int i = idx; i < idx+chunk_size; i++)
        y[i] = x[i] + y[i];
}

void addCpu(int n, float *x, float *y)
{
    for (int i = 0; i < n; i++)
        y[i] = x[i] + y[i];
}

int main(void)
{
    int N = 1 << 24;
    float *x, *y;

    // Allocate Unified Memory – accessible from CPU or GPU
    cudaMallocManaged(&x, N * sizeof(float));
    cudaMallocManaged(&y, N * sizeof(float));

    // initialize x and y arrays on the host
    for (int i = 0; i < N; i++)
    {
        x[i] = 1.0f;
        y[i] = 2.0f;
    }

    auto start = high_resolution_clock::now();
    add<<<16384, 128>>>(N, x, y);
    cudaDeviceSynchronize();
    // addCpu(N, x, y);
    auto stop = high_resolution_clock::now();
    auto duration = duration_cast<microseconds>(stop - start);
  
    // To get the value of duration use the count()
    // member function on the duration object
    std::cout << duration.count() << std::endl;

    // Wait for GPU to finish before accessing on host

    // Check for errors (all values should be 3.0f)
    float maxError = 0.0f;
    for (int i = 0; i < N; i++) {
        // if (fmax(maxError, fabs(y[i] - 3.0f)) != maxError) {
        //     std::cout << fmax(maxError, fabs(y[i] - 3.0f)) << " " << i << std::endl;
        // }
        // std::cout << i << std::endl;
        maxError = fmax(maxError, fabs(y[i] - 3.0f));
    }
    std::cout << "Max error: " << maxError << std::endl;

    // Free memory
    cudaFree(x);
    cudaFree(y);

    return 0;
}