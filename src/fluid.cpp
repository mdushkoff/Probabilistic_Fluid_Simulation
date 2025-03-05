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

void advect(vp_field *vp){
    // TODO: Perform advection
}

void diffuse(vp_field *vp, float viscosity){
    // !!! Needs delta_t and a results buffer !!!

    /*
     * THEORY
     * del_u / del_t = v_const * nabla^2(u)
     * u_n+1 - v_const * delta_t * nabla^2(u_n+1) = u_n
     * u_n+1 - alpha * nabla^2(u_n+1) = u_n
     * --> alpha = v_const * delta_t
     * nabla^2(u) \approx sum(neighbors) - 4(u_i,j)
     * (1 + 4\alpha) * u_n+1 - alpha * sum(neighbors_n+1) = u_n
     * --------------------------------------------
     * SOLUTION
     * u_n+1 = (u_n + alpha * (sum(neighbors_n)) / (1 + 4 * alpha)
     * --> alpha = v_const * delta(t)
     * --> solve the system of equations of all cells with jacobi iteration
     *     (spatial dependencies require nudging toward solution)
     */

    //float alpha = viscosity * delta_t;
    //float beta = 1 + 4 * alpha;
    int width = vp->x, height = vp->y, depth = vp->z;
    float *data = vp->data;
    
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++)
    {
        for (int k = 0; k < depth; k++)
        {
            for (int j = 0; j < height; j++)
            {
                for (int i = 0; i < width; i++)
                {
                    /*data_new[i, j, k] = jacobi(
                        data[i - 1, j, k], data[i + 1, j, k],
                        data[i, j - 1, k], data[i, j + 1, k],
                        alpha, beta, 1.0);*/
                    ;
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
    // Execute operators in order
    advect(vp);
    diffuse(vp, viscosity);
    //addForces(vp_field, forces);  // TODO: eventually add forces
    computePressure(vp);
    subtractPressureGradient(vp);
}