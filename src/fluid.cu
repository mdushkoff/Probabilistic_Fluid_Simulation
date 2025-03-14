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
#define GRID_SIZE (16)

namespace {
    __device__ int asIdx(int i, int j, int k, int width, int channels) {
        return (j * width * channels) + (i * channels) + k;
    }

    __device__ float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
        return (xl+xr+xt+xb+alpha*b)/beta;
    }

    void swapBuffers(vp_field *vp, vp_field *tmp) {
        float* swap = vp->data;
        vp->data = tmp->data;
        tmp->data = swap;
    }
}

__global__ void advect(vp_field *vp, vp_field *vp_out){
    // TODO: Perform advection
}

__global__ void advect_color(vp_field *image, vp_field *out, vp_field *vp, float dt){
    // TODO: Perform color advection
}

__global__ void diffuse(vp_field *vp, vp_field *vp_out, float viscosity, float dt){
    float alpha = viscosity * dt;
    float beta = 1.0 + 4.0 * alpha;

    int w = vp->x, h = vp->y, c = vp->z;
    float *data = vp->data;
    float *data_out = vp_out->data;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= w || j >= h) return;

    // For x and y components of velocity (channels 0 and 1)
    for (int k = 0; k < 2; k++) {
        int iminus = ( (i - 1 + w) % w );
        int iplus  = ( (i + 1)     % w );
        int jminus = ( (j - 1 + h) % h );
        int jplus  = ( (j + 1)     % h );

        // Get top, bottom, left, right elements
        float left = data[asIdx(iminus, j, k, w, c)];
        float right = data[asIdx(iplus, j, k, w, c)];
        float top = data[asIdx(i, jminus, k, w, c)];
        float bottom = data[asIdx(i, jplus, k, w, c)];
        float u_n = data[asIdx(i, j, k, w, c)];

        data_out[asIdx(i, j, k, w, c)] = jacobi(
            alpha * left,
            alpha * right,
            alpha * top,
            alpha * bottom,
            1.0f,
            beta,
            u_n
        );
    }
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
    dim3 blockDim(GRID_SIZE, GRID_SIZE);
    dim3 gridDim((vp->x + blockDim.x - 1)/blockDim.x,
                 (vp->y + blockDim.y - 1)/blockDim.y );

    //advect<<<,BLOCK_SIZE>>>(vp);

    // For diffusion, do Jacobian iteration and swaps on host
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        diffuse<<<gridDim, blockDim>>>(vp, tmp, viscosity, dt);
        if (iter < NUM_JACOBI_ITERS - 1)
            swapBuffers(vp, tmp);
    }

    ////addForces<<<,BLOCK_SIZE>>>(vp_field, forces);  // TODO: eventually add forces
    //computePressure<<<,BLOCK_SIZE>>>(vp);
    //subtractPressureGradient<<<,BLOCK_SIZE>>>(vp);
}

void advect_color_step(vp_field *image, vp_field *itmp, vp_field *vp, float dt){
    // Run advection on the color data
    //advect_color<<<,BLOCK_SIZE>>>(image, itmp, vp, dt);
    
    // Swap buffers

}