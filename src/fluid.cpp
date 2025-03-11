/*
 * 
 */

// System Includes
#include <iostream>

// Local Includes
#include "../includes/fluid.hpp"

inline float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
    return (xl+xr+xt+xb+alpha*b)/beta;
}

inline int asIdx(int i, int j, int k, int height, int channels){
    return (j * height * channels) + (i * channels) + k;
}

void advect(vp_field *vp, vp_field *vp_out, double dt){
    // TODO: Perform advection
    int width = vp->x, height = vp->y, depth = vp->z;
    float *field = vp->data;
    float *new_field = vp_out->data;

    for (int j = 1; j < height - 1; j++) {
        for (int i = 1; i < width - 1; i++) { // Iterate through the inner part in the grid
            int x_idx = asIdx(i, j, 0, height, depth); // Cell index for x
            int y_idx = asIdx(i, j, 1, height, depth); // Cell index for y

            float x_prev = i - dt * field[x_idx] / width; // Compute the where the fluid came from at the previous time stamp // Divided by h so can convert velocity into grid units
            // Offset index to get velocity in y direction
            float y_prev = j - dt * field[y_idx] / height;

            x_prev = std::max(0.5f, std::min<float>(width - 1.5, x_prev)); // Get the center of cooridinates and prevent out-of-bounds errors
            y_prev = std::max(0.5f, std::min<float>(height - 1.5, y_prev));

            int i0 = int(x_prev), j0 = int(y_prev);
            int i1 = i0 + 1, j1 = j0 + 1; // Find the integer grid coordinates surrounding (x, y)
            float sx = x_prev - i0, sy = y_prev - j0; // Represent the fractional distances from (i0, j0)

            // Bilinear Interpolation
            for (int vel_dir = 0; vel_dir < 2; vel_dir++) {
                new_field[asIdx(i, j, vel_dir, height, depth)] = 
                    (1 - sx) * (1 - sy) * field[asIdx(i0, j0, vel_dir, height, depth)] +
                    sx * (1 - sy) * field[asIdx(i1, j0, vel_dir, height, depth)] +
                    (1 - sx) * sy * field[asIdx(i0, j1, vel_dir, height, depth)] +
                    sx * sy * field[asIdx(i1, j1, vel_dir, height, depth)];
            }
        }
    }
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
    // diffuse(vp, viscosity);
    //addForces(vp_field, forces);  // TODO: eventually add forces
    // computePressure(vp);
    // subtractPressureGradient(vp);
}