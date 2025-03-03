/*
 * 
 */

// System Includes
#include <iostream>

// Local Includes
#include "../includes/fluid.hpp"

void advect(vp_field *vp){
    // TODO: Perform advection
}

void diffuse(vp_field *vp, float viscosity){
    // TODO: Perform diffusion
}

void addForces(vp_field *vp, float *forces){
    // TODO: Perform force addition
}

void computePressure(vp_field *vp){
    // TODO: Perform pressure computation
}

void subtractPressureGradient(vp_field *vp){
    // TODO: Perform pressure gradient subtraction
}

void simulate_fluid_step(vp_field *vp, float dt, float viscosity){
    // Execute operators in order
    advect(vp);
    diffuse(vp, viscosity);
    //addForces(vp_field, forces);  // TODO: eventually add forces
    computePressure(vp);
    subtractPressureGradient(vp);
}