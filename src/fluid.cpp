/*
 * 
 */

// System Includes
#include <iostream>
#include <cassert>

// Local Includes
#include "../includes/fluid.hpp"

inline float jacobi(float xl, float xr, float xt, float xb, float alpha, float beta, float b){
    return (xl+xr+xt+xb+alpha*b)/beta;
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
    int width = vp->x, height = vp->y, depth = vp->z;
    float *data = vp->data;
    
    float alpha = -1.0f;
    float beta = 4.0f;
    
    // Temporary buffer to store new pressure values
    float *newPressure = new float[vp->x * vp->y * 4];

    // Compute divergence and store it in D
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            int index = (j * width + i) * 4;
            
            // Neighbor indices for velocity values
            int leftIdx = (j * width + i - 1) * 4;
            int rightIdx = (j * width + i + 1) * 4;
            int topIdx = ((j - 1) * width + i) * 4 + 1;
            int bottomIdx = ((j + 1) * width + i) * 4 + 1;
            int divergenceIdx = index + 3; 

            float uR = ((rightIdx < width * 4) ? data[rightIdx] : 0.0f) - ((leftIdx > -1) ? data[leftIdx] : 0.0f);
            float vT = ((bottomIdx < height * 4) ? data[bottomIdx] : 0.0f) - ((topIdx > -1) ? data[topIdx] : 0.0f);
            
            data[divergenceIdx] = 0.5f * (uR + vT);
        }
    }
    
    for (int iter = 0; iter < NUM_JACOBI_ITERS; iter++) {
        for (int j = 0; j < height; j++) {
            for (int i = 0; i < width; i++) {
                int index = (j * width + i) * 4;
                
                // Neighbor indices for pressure values
                int leftIdx = (j * width + i - 1) * 4 + 2;
                int rightIdx = (j * width + i + 1) * 4 + 2;
                int topIdx = ((j - 1) * width + i) * 4 + 2;
                int bottomIdx = ((j + 1) * width + i) * 4 + 2;
                int divergenceIdx = index + 3;

                float pL =  (leftIdx > -1 ? data[leftIdx] : 0.0f);
                float pR = (rightIdx < width * 4 ? data[rightIdx] : 0.0f);
                float pT = (topIdx > -1 ? data[topIdx] : 0.0f);
                float pB = (rightIdx < height * 4 ? data[bottomIdx] : 0.0f);
                float b = data[divergenceIdx];

                newPressure[index + 2] = jacobi(pL, pR, pT, pB, alpha, beta, b);
            }
        }  
    }

    // Copy back new pressure values to the original data array
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            int index = (j * width + i) * 4;
            data[index + 2] = newPressure[index + 2];
        }
    }
    
    delete[] newPressure;
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
    vp_field test_field;
    test_field.x = 3;
    test_field.y = 3;
    int size = test_field.x * test_field.y * 4;
    test_field.data = new float[size * 4];
    for (int i = 0; i < size; i++){
        test_field.data[i] = 1.0f;
    }

    computePressure(&test_field);

    // Print pressure values in matrix format
    std::cout << "Computed Pressure Values:" << std::endl;
    for (int y = 0; y < test_field.y; ++y) {
        for (int x = 0; x < test_field.x; ++x) {
            std::cout << test_field.data[((y * test_field.x) + x) * 4 + 2]<< " ";
        }
        std::cout << std::endl;
    }
    std::cout << std::endl;

    delete[] test_field.data;
}