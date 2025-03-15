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

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= width || j >= height) return;

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

__global__ void advect_color(vp_field *image, vp_field *out, vp_field *vp, float dt){
    // TODO: Perform color advection
    // Get dimensions
    int iwidth = image->x, iheight = image->y, idepth = image->z;
    int vwidth = vp->x, vheight = vp->y, vdepth = vp->z;

    // Float constants
    float fiwidth = (float)iwidth;
    float fiheight = (float)iheight;
    float fvwidth = (float)vwidth;
    float fvheight = (float)vheight;
    float viw = (fvwidth/fiwidth); // Fractional velocity width to image width
    float vih = (fvheight/fiheight); // Fractional velocity height to image height

    int i = blockIdx.x * blockDim.x + threadIdx.x;
    int j = blockIdx.y * blockDim.y + threadIdx.y;

    if (i >= iwidth || j >= iheight) return; // Out-of-bounds check

    int vi = (int)((float)i * viw);
    int vj = (int)((float)j * vih);

    // Get velocity indices
    int vx_idx = asIdx(vi, vj, 0, vwidth, vdepth); // Cell index for x
    int vy_idx = asIdx(vi, vj, 1, vwidth, vdepth); // Cell index for y

    // Calculate previous coordinate
    float x_prev = (float)i - (dt/viw) * vp->data[vx_idx] / fiwidth; // Compute the where the color came from at the previous time stamp // Divided by h so can convert velocity into grid units
    float y_prev = (float)j - (dt/vih) * vp->data[vy_idx] / fiheight;

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
        itmp->data[asIdx(i, j, k, iwidth, idepth)] = 
            (1 - sx) * (1 - sy) * image->data[asIdx(i0, j0, k, iwidth, idepth)] +
            sx * (1 - sy) * image->data[asIdx(i1, j0, k, iwidth, idepth)] +
            (1 - sx) * sy * image->data[asIdx(i0, j1, k, iwidth, idepth)] +
            sx * sy * image->data[asIdx(i1, j1, k, iwidth, idepth)];
    }
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