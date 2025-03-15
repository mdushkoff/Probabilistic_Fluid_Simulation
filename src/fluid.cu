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

__global__ void computeDivergence(vp_field *vp, vp_field *vp_out, float dt) {
    int w = vp->x, h = vp->y, c = vp->z;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= w || j >= h) return;

    float *data_in = vp->data;
    float *data_out = vp_out->data;

    float gamma = -1.0f / dt;  // Divergence scaling factor

    int iminus = (i - 1 + w) % w;
    int iplus  = (i + 1) % w;
    int jminus = (j - 1 + h) % h;
    int jplus  = (j + 1) % h;

    float uR = data_in[asIdx(iplus, j, 0, w, c)] - data_in[asIdx(iminus, j, 0, w, c)];
    float vT = data_in[asIdx(i, jplus, 1, w, c)] - data_in[asIdx(i, jminus, 1, w, c)];

    float divergence = gamma * (uR + vT);

    data_in[asIdx(i, j, 3, w, c)] = divergence;
    data_out[asIdx(i, j, 3, w, c)] = divergence;
}


__global__ void advect(vp_field *vp, vp_field *vp_out, double dt){
    // TODO: Perform advection
}

__global__ void advect_color(vp_field *image, vp_field *itmp, vp_field *vp, float dt){
    // TODO: Perform color advection
}

__global__ void diffuse(vp_field *vp, vp_field *vp_out, float viscosity, float dt){
    // TODO: Perform diffusion
}

__global__ void addForces(vp_field *vp, float *forces){
    // TODO: Perform force addition
}

__global__ void computePressure(vp_field *vp, vp_field *vp_out, float dt){
    float alpha = 1.0f;
    float beta  = 4.0f;

    int w = vp->x, h = vp->y, c = vp->z;
    float *data = vp->data;
    float *data_out = vp_out->data;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= w || j >= h) return;

    int iminus = (i - 1 + w) % w;
    int iplus  = (i + 1) % w;
    int jminus = (j - 1 + h) % h;
    int jplus  = (j + 1) %h;

    float pL = data[asIdx(iminus, j, 2, w, c)];
    float pR = data[asIdx(iplus,  j, 2, w, c)];
    float pT = data[asIdx(i, jminus, 2, w, c)];
    float pB = data[asIdx(i, jplus,  2, w, c)];

    float b = data[asIdx(i, j, 3, w, c)];

    data_out[asIdx(i, j, 2, w, c)] = jacobi(pL, pR, pT, pB, alpha, beta, b);
}

__global__ void subtractPressureGradient(vp_field *vp, vp_field *vp_out, float dt){
    int w = vp->x, h = vp->y, c = vp->z;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= w || j >= h) return;

    float *data_in = vp->data;
    float *data_out = vp_out->data;

    int iminus = (i - 1 + w) % w;
    int iplus  = (i + 1) % w;
    int jminus = (j - 1 + h) % h;
    int jplus  = (j + 1) % h;

    float pL = data_in[asIdx(iminus, j, 2, w, c)];
    float pR = data_in[asIdx(iplus, j, 2, w, c)];
    float pT = data_in[asIdx(i, jminus, 2, w, c)];
    float pB = data_in[asIdx(i, jplus, 2, w, c)];

    float gradX = (pR - pL) * dt * 0.5f;
    float gradY = (pB - pT) * dt * 0.5f;

    data_out[asIdx(i, j, 0, w, c)] = data_in[asIdx(i, j, 0, w, c)] - gradX; // x velocity
    data_out[asIdx(i, j, 1, w, c)] = data_in[asIdx(i, j, 1, w, c)] - gradY; // y velocity
}

void simulate_fluid_step(vp_field *vp, vp_field *tmp, float dt, float viscosity){
    dim3 blockDim(16, 16);
    dim3 gridDim((vp->x + blockDim.x - 1) / blockDim.x,
                 (vp->y + blockDim.y - 1) / blockDim.y);
    // advect<<<gridDim,blockDim>>>(vp, tmp, dt);
    //diffuse<<<,BLOCK_SIZE>>>(vp, viscosity);
    ////addForces<<<,BLOCK_SIZE>>>(vp_field, forces);  // TODO: eventually add forces
    // computePressuureGradient<<<,BLOCK_SIZE>>>(tmp, vp, dt);
    computeDivergence<<<gridDim,blockDim>>>(vp, tmp, dt);
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        computePressure<<<gridDim, blockDim>>>(vp, tmp, dt);
        if (iter < NUM_JACOBI_ITERS - 1)
            swapBuffers(vp, tmp);
    }
    subtractPressureGradient<<<gridDim,blockDim>>>(tmp, vp, dt);
}

void advect_color_step(vp_field *image, vp_field *itmp, vp_field *vp, float dt){

}