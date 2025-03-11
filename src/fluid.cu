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

__global__ void advect(vp_field *vp, vp_field *vp_out){
    // TODO: Perform advection
}

__global__ void advect_color(vp_field *image, vp_field *out, vp_field *vp, float dt){
    // TODO: Perform color advection
}

__global__ void diffuse(vp_field *vp, vp_field *vp_out, float viscosity, float dt){
    // TODO: Perform diffusion
}

__global__ void addForces(vp_field *vp, float *forces){
    // TODO: Perform force addition
}

__global__ void computePressure(vp_field *vp, vp_field *vp_out){
    // TODO: Perform pressure computation
}

__global__ void subtractPressureGradient(vp_field *vp, vp_field *vp_out){
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