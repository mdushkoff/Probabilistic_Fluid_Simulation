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
#define VP_RANGE (2.0)

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
    std::cerr << "Usage: " << prog << " <n_timesteps> <delta_t> <viscosity> <input_image> <velocity_field> [output_dir]" << std::endl;
}

/*
 * Write an output to the results directory.
 *
 * Inputs:
 *     image - Image data structure
 *     img_handle - A reference to a png_image structure
 *     i - The current iteration number
 *     dir - The output directory
 */
void write_to_output(vp_field *image, png_image *img_handle, int i, char *dir){
    // Write data to a png file
    std::string ext(".png");
    boost::filesystem::path base(dir);
    boost::filesystem::path fn(std::to_string(i) + ext);
    boost::filesystem::path outpath = base / fn;
    std::cout << "[" << i << "] Writing to : " << outpath.c_str() << std::endl;
    write_png_from_array(img_handle, outpath.c_str(), &(image->data));
}

int main(int argc, char *argv[]){
    // Timing variables
    std::chrono::high_resolution_clock::time_point time_start;
    std::chrono::high_resolution_clock::time_point time_end;

    // PNG velocitys
    int status;
    png_image input_image;
    png_image velocity_image;

    // Flags
    int flags = 0;

    // Image data
    vp_field image;   // Input image
    vp_field vp;      // 2D velocity + pressure field (RG,B)
    vp_field itmp;    // Temporary image buffer
    vp_field vtmp;    // Temporary velocity buffer

#ifdef USE_CUDA
    // Device data
    vp_field d_image; // Input image on device memory
    vp_field d_vp;    // 2D Velocity + pressure field (RG,B) on device memory
    vp_field d_itmp;  // Temporary image buffer on device memory
    vp_field d_vtmp;  // Temporary velocity buffer on device memory

    // Streams
    cudaStream_t streams[NUM_STREAMS];
#endif

    // Parse arguments
    if (argc != 6 && argc != 7){
        usage(argv[0]);
        return 1;
    }
    if (argc == 6){
        // Run in timing mode with no file writing
        flags |= 1;
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

    // Try to load the velocity image
    status = read_png_to_array(&input_image, argv[4], &(image.data));
    if (status != 0){
        std::cerr << "Something went wrong reading the input image..." << std::endl;
        return 1;
    }

    // Try to load the initial velocity field
    status = read_png_to_array(&velocity_image, argv[5], &(vp.data));
    if (status != 0){
        std::cerr << "Something went wrong reading the initial velocity field..." << std::endl;
        return 1;
    }

    // Set internal properties
    image.x = input_image.width;
    image.y = input_image.height;
    image.z = NUM_CHANNELS;
    vp.x = velocity_image.width;
    vp.y = velocity_image.height;
    vp.z = NUM_CHANNELS;
    itmp.x = input_image.width;
    itmp.y = input_image.height;
    itmp.z = NUM_CHANNELS;
    vtmp.x = velocity_image.width;
    vtmp.y = velocity_image.height;
    vtmp.z = NUM_CHANNELS;
#ifdef USE_CUDA
    d_image.x = input_image.width;
    d_image.y = input_image.height;
    d_image.z = NUM_CHANNELS;
    d_vp.x = velocity_image.width;
    d_vp.y = velocity_image.height;
    d_vp.z = NUM_CHANNELS;
    d_itmp.x = input_image.width;
    d_itmp.y = input_image.height;
    d_itmp.z = NUM_CHANNELS;
    d_vtmp.x = velocity_image.width;
    d_vtmp.y = velocity_image.height;
    d_vtmp.z = NUM_CHANNELS;
#endif // USE_CUDA

    // Rescale velocity data to be within the range
    for (int x=0; x<vp.x; x++){
        for (int y=0; y<vp.y; y++){
            for (int z=0; z<vp.z; z++){
                int idx = x*vp.y*vp.z+y*vp.z+z;
                float v = vp.data[idx];
                vp.data[idx] = (v*VP_RANGE) - VP_RANGE/2.0;
            }
        }
    }

    // Calculate buffer sizes
    size_t input_bytes = sizeof(float)*input_image.height*input_image.width*NUM_CHANNELS;
    size_t velocity_bytes = sizeof(float)*velocity_image.height*velocity_image.width*NUM_CHANNELS;

    // Initialize temporary buffers
    itmp.data = (float*)malloc(input_bytes);
    vtmp.data = (float*)malloc(velocity_bytes);
    for (int i=0; i<(velocity_image.height*velocity_image.width*NUM_CHANNELS); i++){
        if ((i %4 ) == 3){
            vtmp.data[i] = 1.0;
        }
        else {
            vtmp.data[i] = -1.0;
        }
    }
#ifdef USE_CUDA
    cudaMalloc((void**)&(d_itmp.data), input_bytes);
    cudaMalloc((void**)&(d_vtmp.data), velocity_bytes);

    // Set up streams
    for (int i = 0; i < NUM_STREAMS; ++i){
        cudaStreamCreate(&streams[i]);
    }

    // Allocate device memory
    cudaMalloc((void**)&(d_vp.data), velocity_bytes);

    // Copy memory to the device
    cudaMemcpy(d_vp.data, vp.data, velocity_bytes, cudaMemcpyHostToDevice);
#endif // USE_CUDA

    // Debug TODO: condition on debug output
    std::cout << "Simulating [" << velocity_image.height << " x " << velocity_image.width << "] domain for "
        << n_timesteps << " timesteps at dt=" << delta_t << "..." << std::endl;

    // Loop through all timesteps
    time_start = std::chrono::high_resolution_clock::now();
    for (int i=0; i<n_timesteps; i++){
#ifdef USE_CUDA
        // Simulate the fluid for a single timestep on the device
        simulate_fluid_step(&d_vp, &d_vtmp, delta_t, viscosity);

        // Run advection on the color data
        advect_color_step(&d_image, &d_itmp, &d_vp, delta_t);

        // Copy memory back to the host asynchronously  (TODO: Validate that this doesn't corrupt things)
        cudaMemcpyAsync(vp.data, d_vp.data, velocity_bytes, cudaMemcpyDeviceToHost, streams[i % NUM_STREAMS]);
#else
        // Simulate the fluid for a single timestep
        //std::cout << "STEP[" << i << "]: (" << vp.x << "," << vp.y << "," << vp.z << ") x (" << vtmp.x << "," << vtmp.y << "," << vtmp.z << ")" << std::endl; 
        simulate_fluid_step(&vp, &vtmp, delta_t, viscosity);

        // Run advection on the color data
        advect_color_step(&image, &itmp, &vp, delta_t);
#endif // USE_CUDA

        // Write data to a new output file (don't measure this)
        if (flags == 0){
            write_to_output(&image, &input_image, i, argv[6]);
        }
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
    png_image_free(&velocity_image);

#ifdef USE_CUDA
    // Free pinned host memory
    cudaFreeHost(image.data);
    cudaFreeHost(vp.data);

    // Free CUDA memory
    cudaFree(d_image.data);
    cudaFree(d_vp.data);
    cudaFree(d_itmp.data);
    cudaFree(d_vtmp.data);
#else
    // Free memory
    free(image.data);
    free(vp.data);
    free(itmp.data);
    free(vtmp.data);
#endif // USE_CUDA

    return 0;
}
