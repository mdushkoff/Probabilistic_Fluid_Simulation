/*
 * Perform a probabilistic fluid simulation
 * which computes over a fixed number of timesteps.
 */

// System Includes
#include <chrono>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <unistd.h>  // getopt
#include <utility>

// Conditional Includes
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif // USE_CUDA

// Local Includes
#include "../includes/fluid.hpp"
#include "../includes/utils.hpp"

// Definitions
#define NUM_CHANNELS (4)

// Conditional Definitions
#ifdef USE_CUDA
#define NUM_STREAMS (8)
#endif // USE_CUDA

// Helper Definitions
#define DURATION_MS(x) (std::chrono::duration_cast<std::chrono::milliseconds>(x).count())
#define DURATION_US(x) (std::chrono::duration_cast<std::chrono::microseconds>(x).count())

/*
 * Print the program usage
 */
void usage(char *prog){
    std::cerr << "Usage: " << prog << " <n_timesteps> <delta_t> <viscosity> <initial_field> <output_dir>" << std::endl;
}

int main(int argc, char *argv[]){
    // Timing variables
    std::chrono::high_resolution_clock::time_point time_start;
    std::chrono::high_resolution_clock::time_point time_end;

    // PNG inputs
    int status;
    png_image input_image;

    // Image data
    vp_field vp;  // 2D velocity + pressure field (RG,B)

#ifdef USE_CUDA
    // Device data
    vp_field d_vp;  // 2D Velocity + pressure field (RG,B) on device memory

    // Streams
    cudaStream_t streams[NUM_STREAMS];
#endif

    // Parse arguments
    if (argc != 6){
        usage(argv[0]);
        return 1;
    }
    int n_timesteps = atoi(argv[1]);
    float delta_t = atof(argv[2]);
    float viscosity = atof(argv[3]);
    
    // Validate arguments
    if (n_timesteps <= 0){
        std::cerr << "Timesteps must be greater than 0." << std::endl;
        return 1;
    }
    if (delta_t <= 0){
        std::cerr << "Delta T must be greater than 0." << std::endl;
        return 1;
    }

    // Try to load the initial field
    status = read_png_to_array(&input_image, argv[4], &(vp.data));
    if (status != 0){
        std::cerr << "Something went wrong reading the initial input field..." << std::endl;
        return 1;
    }

    // Set internal properties
    vp.x = input_image.width;
    vp.y = input_image.height;
    vp.z = NUM_CHANNELS;
#ifdef USE_CUDA
    d_vp.x = input_image.width;
    d_vp.y = input_image.height;
    d_vp.z = NUM_CHANNELS;
#endif // USE_CUDA

#ifdef USE_CUDA
    // Set up streams
    for (int i = 0; i < NUM_STREAMS; ++i){
        cudaStreamCreate(&streams[i]);
    }

    // Calculate device size
    size_t image_bytes = sizeof(float)*input_image.height*input_image.width*NUM_CHANNELS;

    // Allocate device memory
    cudaMalloc((void**)&(d_vp.data), image_bytes);

    // Copy memory to the device
    cudaMemcpy(d_vp.data, vp.data, image_bytes, cudaMemcpyHostToDevice);
#endif // USE_CUDA

    // Debug TODO: condition on debug output
    std::cout << "Simulating [" << input_image.height << " x " << input_image.width << "] domain for "
        << n_timesteps << " timesteps at dt=" << delta_t << "..." << std::endl;

    // Loop through all timesteps
    time_start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<n_timesteps; i++){
#ifdef USE_CUDA
        // Simulate the fluid for a single timestep on the device
        simulate_fluid_step(&d_vp, delta_t, viscosity);

        // Copy memory back to the host asynchronously  (TODO: Validate that this doesn't corrupt things)
        cudaMemcpyAsync(vp.data, d_vp.data, image_bytes, cudaMemcpyDeviceToHost, streams[i % NUM_STREAMS]);
#else
        // Simulate the fluid for a single timestep
        simulate_fluid_step(&vp, delta_t, viscosity);
#endif // USE_CUDA

        // TODO: Write data to a new output file (don't measure this)

    }
    time_end = std::chrono::high_resolution_clock::now();
    
    // Calculate the total time
    std::cout << n_timesteps << " timesteps took " << DURATION_US(time_end - time_start) << " us." << std::endl;

#ifdef USE_CUDA
    // Synchronize the device operations
    cudaDeviceSynchronize();
#endif // USE_CUDA

    // Clean up PNG data
    png_image_free(&input_image);

#ifdef USE_CUDA
    // Free pinned host memory
    cudaFreeHost(vp.data);

    // Free CUDA memory
    cudaFree(d_vp.data);
#else
    // Free memory
    free(vp.data);
#endif // USE_CUDA

    return 0;
}
