/*
 * Perform a probabilistic fluid simulation
 * which computes over a fixed number of timesteps.
 */

// System Includes
#include <boost/filesystem.hpp>
#include <chrono>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
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
#define VP_RANGE (10.0)

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

/*
 * Write an output to the results directory.
 *
 * Inputs:
 *     vp -  A 3D array representing the current velocity/pressure values
 *     output_image - A reference to a png_image structure
 *     i - The current iteration number
 *     dir - The output directory
 */
void write_to_output(vp_field *vp, png_image *output_image, int i, char *dir){
    // Copy data to a temporary buffer and normalize it
    int w = vp->x;
    int h = vp->y;
    int c = vp->z;
    size_t image_bytes = sizeof(float)*h*w*c;
    float *out = (float*)malloc(image_bytes);
    memcpy(out, vp->data, image_bytes);

    // Create output image structure
    /*png_image output_image;
    output_image.width = w;
    output_image.height = h;*/
    //output_image.format = PNG_FORMAT_RGBA;

    // Normalize the output array
    for (int x=0; x<w; x++){
        for (int y=0; y<h; y++){
            for (int z=0; z<c; z++){
                int idx = y*h*c+x*c+z;
                out[idx] = (out[idx] + VP_RANGE/2.0) / VP_RANGE;
            }
        }
    }

    // Write data to a png file
    std::string ext(".png");
    boost::filesystem::path base(dir);
    boost::filesystem::path fn(std::to_string(i) + ext);
    boost::filesystem::path outpath = base / fn;
    std::cout << "[" << i << "] Writing to : " << outpath.c_str() << std::endl;
    write_png_from_array(output_image, outpath.c_str(), &out);

    // Free memory
    free(out);
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
    vp_field tmp; // Temporary buffer

#ifdef USE_CUDA
    // Device data
    vp_field d_vp;  // 2D Velocity + pressure field (RG,B) on device memory
    vp_field d_tmp; // Temporary buffer on device memory

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
    tmp.x = input_image.width;
    tmp.y = input_image.height;
    tmp.z = NUM_CHANNELS;
#ifdef USE_CUDA
    d_vp.x = input_image.width;
    d_vp.y = input_image.height;
    d_vp.z = NUM_CHANNELS;
#endif // USE_CUDA

    // Rescale data to be within the range
    for (int x=0; x<vp.x; x++){
        for (int y=0; y<vp.y; y++){
            for (int z=0; z<vp.z; z++){
                int idx = x*vp.y*vp.z+y*vp.z+z;
                float v = vp.data[idx];
                vp.data[idx] = (v*VP_RANGE) - VP_RANGE/2.0;
            }
        }
    }

    // Calculate device size
    size_t image_bytes = sizeof(float)*input_image.height*input_image.width*NUM_CHANNELS;

    // Initialize temporary buffer
    tmp.data = (float*)malloc(image_bytes);
#ifdef USE_CUDA
    cudaMalloc((void**)&(d_tmp.data), image_bytes);

    // Set up streams
    for (int i = 0; i < NUM_STREAMS; ++i){
        cudaStreamCreate(&streams[i]);
    }

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
        simulate_fluid_step(&d_vp, &d_tmp, delta_t, viscosity);

        // Copy memory back to the host asynchronously  (TODO: Validate that this doesn't corrupt things)
        cudaMemcpyAsync(vp.data, d_vp.data, image_bytes, cudaMemcpyDeviceToHost, streams[i % NUM_STREAMS]);
#else
        // Simulate the fluid for a single timestep
        simulate_fluid_step(&vp, &tmp, delta_t, viscosity);
#endif // USE_CUDA

        // TODO: Write data to a new output file (don't measure this)
        write_to_output(&vp, &input_image, i, argv[5]);
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
    cudaFree(d_tmp.data);
#else
    // Free memory
    free(vp.data);
    free(tmp.data);
#endif // USE_CUDA

    return 0;
}
