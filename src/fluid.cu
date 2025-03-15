/*
 * CUDA kernels for computing a fluid simulation
 * on GPU devices.
 */

// System Includes
#include <iostream>

// Local Includes
#include "../includes/fluid.hpp"

// Definitions
#define BLOCK_SIZE (256)

__device__ int asIdx(int i, int j, int k, int width, int channels) {
    return (j * width * channels) + (i * channels) + k;
}

__global__ void advect(vp_field *vp, vp_field *vp_out, double dt){
    // TODO: Perform advection
    int width = vp->x, height = vp->y, depth = vp->z;
    float *field = vp->data;
    float *new_field = vp_out->data;

    float fwidth = (float)width;
    float fheight = (float)height;

    int index_x = blockDim.x * blockIdx.x + threadIdx.x;
    int index_y = blockDim.y * blockIdx.y + threadIdx.y;

    int stride_x = blockDim.x * gridDim.x;
    int stride_y = blockDim.y * gridDim.y;


    if (index_x >= width || index_y >= height) return;

    for (int j = index_y; j < height; j += stride_y) {
        for (int i = index_x; i < width; i += stride_x) { 
            int x_idx = asIdx(i, j, 0, width, depth);
            int y_idx = asIdx(i, j, 1, width, depth);

            float x_prev = (float)i - dt * field[x_idx] / fwidth;
            float y_prev = (float)j - dt * field[y_idx] / fheight;

            float xrem = fmodf(fmodf(x_prev, fwidth) + fwidth, fwidth);
            float yrem = fmodf(fmodf(y_prev, fheight) + fheight, fheight);
            x_prev = xrem;
            y_prev = yrem;

            int i0 = int(x_prev), j0 = int(y_prev);
            int i1 = (i0 + 1) % width, j1 = (j0 + 1) % height;
            float sx = x_prev - (float)i0, sy = y_prev - (float)j0;

            for (int vel_dir = 0; vel_dir < 2; vel_dir++) {
                new_field[asIdx(i, j, vel_dir, width, depth)] =
                    (1 - sx) * (1 - sy) * field[asIdx(i0, j0, vel_dir, width, depth)] +
                    sx * (1 - sy) * field[asIdx(i1, j0, vel_dir, width, depth)] +
                    (1 - sx) * sy * field[asIdx(i0, j1, vel_dir, width, depth)] +
                    sx * sy * field[asIdx(i1, j1, vel_dir, width, depth)];
            }
        }
    }
}

__global__ void advect_color(vp_field *image, vp_field *out, vp_field *vp, float dt){
    // TODO: Perform color advection
    // Get dimensions
}

__global__ void diffuse(vp_field *vp, vp_field *vp_out, float viscosity, float dt){
    // TODO: Perform diffusion
}

__global__ void addForces(vp_field *vp, float *forces){
    // TODO: Perform force addition
}

__global__ void computePressure(vp_field *vp, vp_field *vp_out, float dt){
    // TODO: Perform pressure computation
}

__global__ void subtractPressureGradient(vp_field *vp, vp_field *vp_out, float dt){
    // TODO: Perform pressure gradient subtraction
}

void simulate_fluid_step(vp_field *vp, vp_field *tmp, float dt, float viscosity){
    //advect<<<,BLOCK_SIZE>>>(vp);
    //diffuse<<<,BLOCK_SIZE>>>(vp, viscosity);
    ////addForces<<<,BLOCK_SIZE>>>(vp_field, forces);  // TODO: eventually add forces
    //computePressure<<<,BLOCK_SIZE>>>(vp);
    //subtractPressureGradient<<<,BLOCK_SIZE>>>(vp);
}

void advect_color_step(vp_field *image, vp_field *itmp, vp_field *vp, float dt){
    // Run advection on the color data
    //advect_color<<<,BLOCK_SIZE>>>(image, itmp, vp, dt);
    
    // Swap buffers

}
