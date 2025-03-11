/*
 * 
 */

// System Includes
#include <iostream>

// Local Includes
#include "../includes/fluid.hpp"

namespace {
    //int NUM_CHANNELS = 4;

    inline int asIdx(int i, int j, int k, int height, int channels) {
        return (j * height * channels) + (i * channels) + k;
    }

    inline float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
        return (xl+xr+xt+xb+alpha*b)/beta;
    }
}

void advect(vp_field *vp, vp_field *vp_out){
    // Perform advection
    int w = vp->x, h = vp->y, c = vp->z;
    for (int x=0; x<w; x++){
        for (int y=0; y<h; y++){
            for (int z=0; z<c; z++){
                vp_out->data[asIdx(x, y, z, h, c)] = vp->data[asIdx(x, y, z, h, c)];  // TODO: Perform advection (This is just identity)
            }
        }
    }
}

void diffuse(vp_field *vp, vp_field *vp_out, float viscosity, float dt){
    // THEORY
    // del_u / del_t = v_const * nabla^2(u)
    // u_n+1 - v_const * delta_t * nabla^2(u_n+1) = u_n
    // u_n+1 - alpha * nabla^2(u_n+1) = u_n
    //     --> alpha = v_const * delta_t
    // nabla^2(u) \approx sum(neighbors) - 4(u_i,j)
    // (1 + 4\alpha) * u_n+1 - alpha * sum(neighbors_n+1) = u_n
    // --------------------------------------------
    // SOLUTION
    // u_n+1 = (u_n + alpha * (sum(neighbors_n)) / (1 + 4 * alpha)
    //     --> alpha = v_const * delta(t)
    //     --> solve the system of equations of all cells with jacobi iteration
    //         (spatial dependencies require nudging toward solution)

    float alpha = viscosity * dt;
    float beta = 1.0 + 4.0 * alpha;

    int w = vp->x, h = vp->y, c = vp->z;
    float *data = vp->data, *data_out = vp_out->data;
    
    // Jacobian iteration on outside to update all spatial dependencies
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++)
    {
        // Iterate over 2D grid
        for (int i = 0; i < w; i++)
        {
            for (int j = 0; j < h; j++)
            {
                // For x and y components of velocity (channels 0 and 1)
                for (int k = 0; k < 2; k++)
                {
                    // Get top, bottom, left, right elements
                    float left, right, top, bottom;
                    if (i != 0){
                        left = data[asIdx(i - 1, j, k, h, c)];
                    }
                    else {
                        left = 0.0; // TODO: Figure out boundary
                    }
                    if (i != (w-1)){
                        right = data[asIdx(i + 1, j, k, h, c)];
                    }
                    else {
                        right = 0.0; // TODO: Figure out boundary
                    }
                    if (j != 0){
                        top = data[asIdx(i, j - 1, k, h, c)];
                    }
                    else {
                        top = 0.0; // TODO: Figure out boundary
                    }
                    if (j != (h-1)){
                        bottom = data[asIdx(i, j + 1, k, h, c)];
                    }
                    else {
                        bottom = 0.0;  // TODO: Figure out boundary
                    }

                    // Solve Laplacian
                    data_out[asIdx(i, j, k, h, c)] = jacobi(
                        left,
                        right,
                        top,
                        bottom,
                        alpha, beta, 1.0);
                }
            }
        }

        // Swap buffers
        if (iter != (NUM_JACOBI_ITERS-1)){
            float *tp = vp_out->data;
            vp_out->data = vp->data;
            vp->data = tp;
            data = vp->data;
            data_out = vp_out->data;
        }
    }
}

void addForces(vp_field *vp, float *forces){
    // TODO: Perform force addition
    int w = vp->x, h = vp->y, c = vp->z;
    for (int x=0; x<w; x++){
        for (int y=0; y<h; y++){
            for (int z=0; z<c; z++){
                
            }
        }
    }
}

void computePressure(vp_field *vp, vp_field *vp_out){
    // Perform pressure computation
    int w = vp->x, h = vp->y, c = vp->z;
    for (int x=0; x<w; x++){
        for (int y=0; y<h; y++){
            for (int z=0; z<c; z++){
                vp_out->data[asIdx(x, y, z, h, c)] = vp->data[asIdx(x, y, z, h, c)]; // TODO: Perform pressure computation (This is just identity)
            }
        }
    }
}

void subtractPressureGradient(vp_field *vp, vp_field *vp_out){
    // Perform pressure gradient subtraction
    int w = vp->x, h = vp->y, c = vp->z;
    for (int x=0; x<w; x++){
        for (int y=0; y<h; y++){
            for (int z=0; z<c; z++){
                vp_out->data[asIdx(x, y, z, h, c)] = vp->data[asIdx(x, y, z, h, c)]; // TODO: Perform pressure gradient subtraction (This is just identity)
            }
        }
    }
}

void simulate_fluid_step(vp_field *vp, vp_field *tmp, float dt, float viscosity){
    // Execute operators in order (swapping buffers each step)
    advect(vp, tmp);
    diffuse(tmp, vp, viscosity, dt);
    //addForces(vp_field, forces);  // TODO: eventually add forces
    computePressure(vp, tmp);
    subtractPressureGradient(tmp, vp);
}