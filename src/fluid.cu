/*
 * CUDA kernels for computing a fluid simulation
 * on GPU devices.
 */

// System Includes
#include <iostream>
#include <stdio.h>

// Local Includes
#include "../includes/fluid.hpp"

// Definitions
#define BLOCK_SIZE_X (16)
#define BLOCK_SIZE_Y (16)

#define IX(i,j,k,w,c) ((j*w*c)+(i*c)+k)

namespace {
    __device__ int asIdx(int i, int j, int k, int width, int channels) {
        return (j * width * channels) + (i * channels) + k;
    }

    __device__ float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
        return (xl+xr+xt+xb+alpha*b)/beta;
    }

    void swapBuffers(float **vp, float **tmp) {
        /*float* swap = vp->data;
        vp->data = tmp->data;
        tmp->data = swap;*/
        float *swap = *vp;
        (*vp) = (*tmp);
        (*tmp) = swap;
    }
}

__global__ void computeDivergence(float *vp, float *vp_out, float dt, int vx, int vy, int vz) {
    int w = vx, h = vy, c = vz;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= w || j >= h) return;

    float *data_in = (vp);
    float *data_out = (vp_out);

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

__global__ void advect(float *vp, float *vp_out, float dt, int vx, int vy, int vz){
    // Perform advection
    int width = vx, height = vy, depth = vz;
    float *field = (vp);
    float *new_field = (vp_out);

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

__global__ void advect_color(float *image, float *itmp, float *vp, float dt, int ix, int iy, int iz, int vx, int vy, int vz){
    // Perform color advection
    // Get dimensions
    int iwidth = ix, iheight = iy, idepth = iz;
    int vwidth = vx, vheight = vy, vdepth = vz;

    // Float constants
    float fiwidth = (float)iwidth;
    float fiheight = (float)iheight;
    float fvwidth = (float)vwidth;
    float fvheight = (float)vheight;
    float viw = (fvwidth/fiwidth); // Fractional velocity width to image width
    float vih = (fvheight/fiheight); // Fractional velocity height to image height

    int index_x = blockDim.x * blockIdx.x + threadIdx.x;
    int index_y = blockDim.y * blockIdx.y + threadIdx.y;

    int stride_x = blockDim.x * gridDim.x;
    int stride_y = blockDim.y * gridDim.y;


    if (index_x >= iwidth || index_y >= iheight) return;

    // Loop through all pixels and advect
    for (int j = index_y; j < iheight; j += stride_y) {
        for (int i = index_x; i < iwidth; i += stride_x) {
            // Get velocity adjusted coordinates
            int vi = (int)((float)i * viw);
            int vj = (int)((float)j * vih);

            // Get velocity indices
            int vx_idx = asIdx(vi, vj, 0, vwidth, vdepth); // Cell index for x
            int vy_idx = asIdx(vi, vj, 1, vwidth, vdepth); // Cell index for y

            // Calculate previous coordinate
            float x_prev = (float)i - (dt/viw) * (vp)[vx_idx] / fiwidth; // Compute the where the color came from at the previous time stamp // Divided by h so can convert velocity into grid units
            float y_prev = (float)j - (dt/vih) * (vp)[vy_idx] / fiheight;

            // Clamp mode
            //x_prev = std::max(0.5f, std::min<float>(fiwidth - 1.5, x_prev)); // Get the center of cooridinates and prevent out-of-bounds errors
            //y_prev = std::max(0.5f, std::min<float>(fiheight - 1.5, y_prev));

            // Wrapping mode
            x_prev = fmod(fmod(x_prev, fiwidth)+fiwidth, fiwidth);
            y_prev = fmod(fmod(y_prev, fiheight)+fiheight, fiheight);
            //float xrem = fmod(x_prev, fiwidth);
            //float yrem = fmod(y_prev, fiheight);
            //x_prev = (x_prev >= 0) ? xrem : xrem + fiwidth;
            //y_prev = (y_prev >= 0) ? yrem : yrem + fiheight;

            int i0 = int(x_prev), j0 = int(y_prev);
            int i1 = (i0 + 1) % iwidth, j1 = (j0 + 1) % iheight; // Find the integer grid coordinates surrounding (x, y)
            float sx = x_prev - (float)i0, sy = y_prev - (float)j0; // Represent the fractional distances from (i0, j0)

            // Apply velocity to each channel
            for (int k=0; k<idepth; k++){
                // Bilinear interpolation
                (itmp)[asIdx(i, j, k, iwidth, idepth)] = 
                    (1 - sx) * (1 - sy) * (image)[asIdx(i0, j0, k, iwidth, idepth)] +
                    sx * (1 - sy) * (image)[asIdx(i1, j0, k, iwidth, idepth)] +
                    (1 - sx) * sy * (image)[asIdx(i0, j1, k, iwidth, idepth)] +
                    sx * sy * (image)[asIdx(i1, j1, k, iwidth, idepth)];
            }
        }
    }
}

__global__ void diffuse(float *vp, float *vp_out, float viscosity, float dt, int vx, int vy, int vz){
    float alpha = viscosity * dt;
    float beta = 1.0 + 4.0 * alpha;

    int w = vx, h = vy, c = vz;
    float *data = (vp);
    float *data_out = (vp_out);

    int tid_i = blockIdx.x * blockDim.x + threadIdx.x;
    int tid_j = blockIdx.y * blockDim.y + threadIdx.y;

    int stride_i = blockDim.x * gridDim.x;
    int stride_j = blockDim.y * gridDim.y;

    for (int i = tid_i; i < w; i += stride_i) {
        for (int j = tid_j; j < h; j += stride_j) {
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
    }
}

__global__ void addForces(float *vp, float *forces, int vx, int vy, int vz){
    // TODO: Perform force addition
}

__global__ void computePressure(float *vp, float *vp_out, float dt, int vx, int vy, int vz){
    float alpha = 1.0f;
    float beta  = 4.0f;

    int w = vx, h = vy, c = vz;
    float *data = (vp);
    float *data_out = (vp_out);

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

__global__ void subtractPressureGradient(float *vp, float *vp_out, float dt, int vx, int vy, int vz){
    int w = vx, h = vy, c = vz;

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= w || j >= h) return;

    float *data_in = (vp);
    float *data_out = (vp_out);

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

void simulate_fluid_step(float **vp, float **tmp, float dt, float viscosity, int vx, int vy, int vz){
    dim3 blockDim(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    dim3 gridDim((vx + blockDim.x - 1) / blockDim.x,
                 (vy + blockDim.y - 1) / blockDim.y);

    advect<<<gridDim,blockDim>>>(*vp, *tmp, dt, vx, vy, vz);
    swapBuffers(vp, tmp);

    // For diffusion, do Jacobian iteration and swaps on host
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        diffuse<<<gridDim, blockDim>>>(*vp, *tmp, viscosity, dt, vx, vy, vz);
        if (iter < NUM_JACOBI_ITERS - 1)
            swapBuffers(vp, tmp);
    }

    // If odd number of iters, swap so result is in vp
    if (NUM_JACOBI_ITERS % 2)
        swapBuffers(vp, tmp);

    //addForces<<<gridDim,blockDim>>>(vp_field, forces);  // TODO: eventually add forces
    computeDivergence<<<gridDim,blockDim>>>(*vp, *tmp, dt, vx, vy, vz);
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        computePressure<<<gridDim, blockDim>>>(*vp, *tmp, dt, vx, vy, vz);
        if (iter < NUM_JACOBI_ITERS - 1)
            swapBuffers(vp, tmp);
    }
    subtractPressureGradient<<<gridDim,blockDim>>>(*tmp, *vp, dt, vx, vy, vz);
}

void advect_color_step(float **image, float **itmp, float **vp, float dt, int ix, int iy, int iz, int vx, int vy, int vz){
    dim3 blockDim(BLOCK_SIZE_X, BLOCK_SIZE_Y);
    dim3 gridDim((ix + blockDim.x - 1) / blockDim.x,
                 (iy + blockDim.y - 1) / blockDim.y);
    advect_color<<<gridDim,blockDim>>>(*image, *itmp, *vp, dt, ix, iy, iz, vx, vy, vz);
    swapBuffers(image, itmp);
}
