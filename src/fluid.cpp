/*
 * 
 */

// System Includes
#include <iostream>
#include <cassert>

// Local Includes
#include "../includes/fluid.hpp"

namespace {
    int NUM_CHANNELS = 4;

    inline int asIdx(int i, int j, int k, int height) {
        return (i * height * NUM_CHANNELS) + (j * NUM_CHANNELS) + k;
    }

    inline float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
        return (xl+xr+xt+xb+alpha*b)/beta;
    }
}

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
    int w = vp->x, h = vp->y, d = vp->z;
    float *data = vp->data;
    
    float alpha = -1.0f;
    float beta = 4.0f;
    
    // Temporary buffer to store new pressure values
    vp_field *newPressure = new vp_field;
    newPressure->x = w;
    newPressure->y = h;
    newPressure->z = d;
    newPressure->data = (float *)malloc(sizeof(float) * newPressure->x * newPressure->y * newPressure->z);

    // Compute divergence and store it in D
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {            
            float uR = (i < w - 1 ? data[asIdx(i + 1, j, 0, h)] : 1.0f) - (i > 0 ? data[asIdx(i - 1, j, 0, h)] : 0.0f);
            float vT = (j < h - 1 ? data[asIdx(i, j + 1, 1, h)] : 1.0f) - (j > 0 ? data[asIdx(i, j - 1, 1, h)] : 0.0f);
            
            data[asIdx(i, j, 3, h)] = 0.5f * (uR + vT);
        }
    }
    
    // Check
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                // Pressure Values
                float pL =  (i > 0 ? data[asIdx(i - 1, j, 2, h)] : 0.0f);
                float pR = (i < w - 1 ? data[asIdx(i + 1, j, 2, h)] : 0.0f);
                float pT = (j > 0 ? data[asIdx(i, j - 1, 2, h)] : 0.0f);
                float pB = (j < h - 1 ? data[asIdx(i, j + 1, 2, h)] : 0.0f);
                float b = data[asIdx(i, j, 3, h)];

                newPressure->data[asIdx(i, j, 2, h)] = jacobi(pL, pR, pT, pB, alpha, beta, b);
            }
        } 
        // Copy back new pressure values to the original data array
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                data[asIdx(i, j, 2, h)] = newPressure->data[asIdx(i, j, 2, h)];
            }
        }
    }

    
    delete newPressure;
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

void test_computePressure() {
    // Define grid size
    vp_field test_field;
    test_field.x = 3; 
    test_field.y = 3; 
    test_field.z = 3;
    int size = test_field.x * test_field.y * 4; // 4 channels per grid cell
    
    test_field.data = new float[size];

    // Initialize all values to zero
    for (int i = 0; i < size; i++) {
        test_field.data[i] = 0.0f;
    }

    // Set initial velocity field (channels 0 and 1)
    for (int i = 0; i < test_field.x; i++) {
        for (int j = 0; j < test_field.y; j++) {
            test_field.data[asIdx(i, j, 0, test_field.y)] = 1.0f;
            test_field.data[asIdx(i, j, 1, test_field.y)] = 1.0f;
            test_field.data[asIdx(i, j, 2, test_field.y)] = 1.0f;
        }
    }

    computePressure(&test_field);

    // Print pressure values in matrix format
    std::cout << "Computed Pressure Values:" << std::endl;
    for (int j = 0; j < test_field.y; ++j) {
        for (int i = 0; i < test_field.x; ++i) {
            float pressure = test_field.data[asIdx(i, j, 2, test_field.y)];
            std::cout << pressure << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Clean up memory
    delete[] test_field.data;
}