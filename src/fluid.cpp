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
        return (j * height * NUM_CHANNELS) + (i * NUM_CHANNELS) + k;
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

void computePressure(vp_field *vp, vp_field *vp_out){
    // TODO: Perform pressure computation
    int w = vp->x, h = vp->y, d = vp->z;
    float *data_in = vp->data;
    
    float alpha = -1.0f;
    float beta = 4.0f;
    

    // Compute divergence and store it in D
    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {            
            float uR = (i < w - 1 ? data_in[asIdx(i + 1, j, 0, w)] : 0.0f) - (i > 0 ? data_in[asIdx(i - 1, j, 0, w)] : 0.0f);
            float vT = (j < h - 1 ? data_in[asIdx(i, j + 1, 1, w)] : 0.0f) - (j > 0 ? data_in[asIdx(i, j - 1, 1, w)] : 0.0f);
            
            data_in[asIdx(i, j, 3, w)] = 0.5f * (uR + vT);
        }
    }
    
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        for (int i = 0; i < w; i++) {
            for (int j = 0; j < h; j++) {
                // Pressure Values
                float pL = (i > 0     ? data_in[asIdx(i - 1, j, 2, w)] : 0.0f);
                float pR = (i < w - 1 ? data_in[asIdx(i + 1, j, 2, w)] : 0.0f);
                float pT = (j > 0     ? data_in[asIdx(i, j - 1, 2, w)] : 0.0f);
                float pB = (j < h - 1 ? data_in[asIdx(i, j + 1, 2, w)] : 0.0f);
                float b = data_in[asIdx(i, j, 3, w)];

                vp_out->data[asIdx(i, j, 2, w)] = jacobi(pL, pR, pT, pB, alpha, beta, b);
            }

        } 

        if (iter < NUM_JACOBI_ITERS - 1) {
            float* tp = vp_out->data;
            vp_out->data = data_in;
            data_in = tp;
        }
    }
}

void subtractPressureGradient(vp_field *vp, vp_field *vp_out) {
    int w = vp->x, h = vp->y, d = vp->z;
    float *data_in = vp->data;

    memcpy(vp_out->data, data_in, sizeof(float) * w * h * d);

    for (int i = 0; i < w; i++) {
        for (int j = 0; j < h; j++) {
            // Compute pressure differences (gradient)
            float pR = (i < w - 1) ? data_in[asIdx(i + 1, j, 2, w)] : 0.0f;
            float pL = (i > 0) ? data_in[asIdx(i - 1, j, 2, w)] : 0.0f;
            float pB = (j < h - 1) ? data_in[asIdx(i, j + 1, 2, w)] : 0.0f;
            float pT = (j > 0) ? data_in[asIdx(i, j - 1, 2, w)] : 0.0f;

            // Compute gradient components
            float gradX = (pR - pL)/2.0f;
            float gradY = (pB - pT)/2.0f;

            // Subtract the gradient from velocity
            vp_out->data[asIdx(i, j, 0, w)] = *data_in - gradX;
            vp_out->data[asIdx(i, j, 1, w)] = *data_in - gradY;
        }
    }
}

void simulate_fluid_step(vp_field *vp, vp_field *vp_out, float dt, float viscosity){
    // Execute operators in order
    advect(vp);
    diffuse(vp, viscosity);
    //addForces(vp_field, forces);  // TODO: eventually add forces
    computePressure(vp, vp_out);
    subtractPressureGradient(vp, vp_out);
}

void test_computePressure() {
    // Define grid size
    vp_field test_field;
    test_field.x = 3; 
    test_field.y = 3; 
    test_field.z = 4;
    int size = test_field.x * test_field.y * test_field.z;
    test_field.data = new float[size];

    vp_field *newPressure = new vp_field;
    newPressure->x = test_field.x;
    newPressure->y = test_field.y;
    newPressure->z = test_field.z;
    newPressure->data = new float[size];
    

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

    computePressure(&test_field, newPressure);

    std::cout << "Computed Pressure Values:" << std::endl;
    for (int j = 0; j < newPressure->y; ++j) {
        for (int i = 0; i < newPressure->x; ++i) {
            float pressure = newPressure->data[asIdx(i, j, 2, newPressure->y)];
            std::cout << pressure << " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;
}

void test_subtractPressureGradient() {
    // Define grid size
    vp_field test_field;
    test_field.x = 2; 
    test_field.y = 2; 
    test_field.z = 4;
    int size = test_field.x * test_field.y * test_field.z;
    test_field.data = new float[size];

    vp_field *res_field = new vp_field;
    res_field->x = test_field.x;
    res_field->y = test_field.y;
    res_field->z = test_field.z;
    res_field->data = new float[size];
    

    // Initialize all values to zero
    for (int i = 0; i < size; i++) {
        test_field.data[i] = 0.0f;
    }

    // Set initial velocity field (channels 0 and 1)
    for (int j = 0; j < test_field.y; j++) {
        for (int i = 0; i < test_field.x; i++) {
            test_field.data[asIdx(i, j, 0, test_field.x)] = 1.0f;
            test_field.data[asIdx(i, j, 1, test_field.x)] = 1.0f;
            test_field.data[asIdx(i, j, 2, test_field.x)] = 1.0f;
        }
    }

    subtractPressureGradient(&test_field, res_field);

    std::cout << "Computed Pressure Values:" << std::endl;
    for (int j = 0; j < res_field->y; ++j) {
        for (int i = 0; i < res_field->x; ++i) {
            std::cout << "[";
            for (int z = 0; z < res_field->z; ++z) {
                std::cout << res_field->data[asIdx(i, j, z, res_field->y)] << ", ";
            }
            std::cout << "]\t";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    // Clean up memory
    delete[] test_field.data;
    delete[] res_field->data;
    delete res_field;
}