/*
 * 
 */

// System Includes
#include <iostream>

// Local Includes
#include "../includes/fluid.hpp"

namespace {
    //int NUM_CHANNELS = 4;

    inline int asIdx(int i, int j, int k, int width, int channels) {
        return (j * width * channels) + (i * channels) + k;
    }

    inline float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
        return (xl+xr+xt+xb+alpha*b)/beta;
    }
}

void advect(vp_field *vp, vp_field *vp_out, double dt){
    // TODO: Perform advection
    int width = vp->x, height = vp->y, depth = vp->z;
    float *field = vp->data;
    float *new_field = vp_out->data;

    // Store floating point constants
    float fwidth = (float)width;
    float fheight = (float)height;

    for (int j = 1; j < height - 1; j++) {
        for (int i = 1; i < width - 1; i++) { // Iterate through the inner part in the grid
            int x_idx = asIdx(i, j, 0, width, depth); // Cell index for x
            int y_idx = asIdx(i, j, 1, width, depth); // Cell index for y

            float x_prev = (float)i - dt * field[x_idx] / fwidth; // Compute the where the fluid came from at the previous time stamp // Divided by h so can convert velocity into grid units
            // Offset index to get velocity in y direction
            float y_prev = (float)j - dt * field[y_idx] / fheight;

            x_prev = std::max(0.5f, std::min<float>(fwidth - 1.5, x_prev)); // Get the center of cooridinates and prevent out-of-bounds errors
            y_prev = std::max(0.5f, std::min<float>(fheight - 1.5, y_prev));

            int i0 = int(x_prev), j0 = int(y_prev);
            int i1 = i0 + 1, j1 = j0 + 1; // Find the integer grid coordinates surrounding (x, y)
            float sx = x_prev - (float)i0, sy = y_prev - (float)j0; // Represent the fractional distances from (i0, j0)

            // Bilinear Interpolation
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
                    left = (i != 0) ? data[asIdx(i - 1, j, k, w, c)] : 0.0; // TODO: Figure out boundary
                    right = (i != (w-1)) ? data[asIdx(i + 1, j, k, w, c)] : 0.0; // TODO: Figure out boundary
                    top = (j != 0) ? data[asIdx(i, j - 1, k, w, c)] : 0.0; // TODO: Figure out boundary
                    bottom = (j != (h-1)) ? data[asIdx(i, j + 1, k, w, c)] : 0.0; // TODO: Figure out boundary
                    float u_n = data[asIdx(i, j, k, w, c)];

                    // Solve Laplacian
                    data_out[asIdx(i, j, k, w, c)] = jacobi(
                        alpha * left,
                        alpha * right,
                        alpha * top,
                        alpha * bottom,
                        1.0f,
                        beta,
                        u_n);
                }
            }
        }

        // Swap buffers
        float *tp = vp_out->data;
        vp_out->data = vp->data;
        vp->data = tp;
        data = vp->data;
        data_out = vp_out->data;
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
    int w = vp->x, h = vp->y, d = vp->z;
    float *data_in = vp->data;
    
    float alpha = -1.0f;
    float beta = 4.0f;
    

    // Compute divergence and store it in D
    // for (int i = 0; i < w; i++) {
    //     for (int j = 0; j < h; j++) {            
    //         float uR = (i < w - 1 ? data_in[asIdx(i + 1, j, 0, w, d)] : 0.0f) - (i > 0 ? data_in[asIdx(i - 1, j, 0, w, d)] : 0.0f);
    //         float vT = (j < h - 1 ? data_in[asIdx(i, j + 1, 1, w, d)] : 0.0f) - (j > 0 ? data_in[asIdx(i, j - 1, 1, w, d)] : 0.0f);
            
    //         data_in[asIdx(i, j, 3, w, d)] = 0.5f * (uR + vT);
    //     }
    // }
    
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                // Compute divergence
                float uR = (i < w - 1 ? data_in[asIdx(i + 1, j, 0, w, d)] : 0.0f) - (i > 0 ? data_in[asIdx(i - 1, j, 0, w, d)] : 0.0f);
                float vT = (j < h - 1 ? data_in[asIdx(i, j + 1, 1, w, d)] : 0.0f) - (j > 0 ? data_in[asIdx(i, j - 1, 1, w, d)] : 0.0f);
            
                // Pressure Values
                float pL = (i > 0     ? data_in[asIdx(i - 1, j, 2, w, d)] : 0.0f);
                float pR = (i < w - 1 ? data_in[asIdx(i + 1, j, 2, w, d)] : 0.0f);
                float pT = (j > 0     ? data_in[asIdx(i, j - 1, 2, w, d)] : 0.0f);
                float pB = (j < h - 1 ? data_in[asIdx(i, j + 1, 2, w, d)] : 0.0f);
                float b = 0.5f * (uR + vT);;

                vp_out->data[asIdx(i, j, 2, w, d)] = jacobi(pL, pR, pT, pB, alpha, beta, b);
            }
        }

        // Swap buffers
        float* tp = vp_out->data;
        vp_out->data = vp->data;
        vp->data = tp;
        data_in = vp->data;
    }
}

void subtractPressureGradient(vp_field *vp, vp_field *vp_out) {
    int w = vp->x, h = vp->y, d = vp->z;
    float *data_in = vp->data;

    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            // Compute pressure differences (gradient)
            float pR = (i < w - 1) ? data_in[asIdx(i + 1, j, 2, w, d)] : 0.0f;
            float pL = (i > 0) ? data_in[asIdx(i - 1, j, 2, w, d)] : 0.0f;
            float pB = (j < h - 1) ? data_in[asIdx(i, j + 1, 2, w, d)] : 0.0f;
            float pT = (j > 0) ? data_in[asIdx(i, j - 1, 2, w, d)] : 0.0f;

            // Compute gradient components
            float gradX = (pR - pL)/2.0f;
            float gradY = (pB - pT)/2.0f;

            // Subtract the gradient from velocity
            vp_out->data[asIdx(i, j, 0, w, d)] = data_in[asIdx(i, j, 0, w, d)] - gradX;
            vp_out->data[asIdx(i, j, 1, w, d)] = data_in[asIdx(i, j, 1, w, d)] - gradY;
        }
    }
}

void simulate_fluid_step(vp_field *vp, vp_field *tmp, float dt, float viscosity){
    // Execute operators in order (swapping buffers each step)
    advect(vp, tmp, dt);
    diffuse(tmp, vp, viscosity, dt);
    //addForces(vp_field, forces);  // TODO: eventually add forces
    computePressure(vp, tmp);
    subtractPressureGradient(tmp, vp);
    // memcpy(tmp->data, vp->data, sizeof(float) * vp->x * vp->y * vp->z);
    vp->data = tmp->data;
}