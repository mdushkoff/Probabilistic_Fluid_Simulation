/*
 * 
 */

// System Includes
#include <iostream>

// Local Includes
#include "../includes/fluid.hpp"

namespace {
    inline int asIdx(int i, int j, int k, int w, int h) {
        return (i * w) + (j * w) * h + k;
    }

    inline float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
        return (xl+xr+xt+xb+alpha*b)/beta;
    }
}

void advect(vp_field *vp){
    // TODO: Perform advection
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
    float beta = 1 + 4 * alpha;

    int w = vp->x, h = vp->y;
    float *data = vp->data, *data_out = vp_out->data;
    
    // For x and y components of velocity
    for (int k = 0; k < 2; k++)
    {
        for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++)
        {
            // Iterate over 2D grid
            for (int j = 0; j < h; j++)
            {
                for (int i = 0; i < w; i++)
                {
                    // Solve Laplacian
                    data_out[asIdx(i, j, k, w, h)] = jacobi(
                        data[asIdx(i - 1, j, k, w, h)],
                        data[asIdx(i + 1, j, k, w, h)],
                        data[asIdx(i, j - 1, k, w, h)],
                        data[asIdx(i, j + 1, k, w, h)],
                        alpha, beta, 1.0);
                }
            }
        }
    }
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
    // Create temporary buffer for results
    vp_field vp_out;
    vp_out.x = vp->x;
    vp_out.y = vp->y;
    vp_out.z = vp->z;
    vp_out.data = (float *)malloc(sizeof(float) * vp_out.x * vp_out.y * vp_out.z);

    // Execute operators in order (swapping buffers each step)
    advect(vp);
    diffuse(&vp_out, vp, viscosity, dt);
    //addForces(vp_field, forces);  // TODO: eventually add forces
    computePressure(vp);
    subtractPressureGradient(vp);

    free(vp_out.data);
}