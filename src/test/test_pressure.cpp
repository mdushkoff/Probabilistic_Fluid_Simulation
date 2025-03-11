#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib>

#include "../../includes/fluid.hpp"

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
            test_field.data[asIdx(i, j, 0, test_field.x, test_field.z)] = 1.0f;
            test_field.data[asIdx(i, j, 1, test_field.x, test_field.z)] = 1.0f;
            test_field.data[asIdx(i, j, 2, test_field.x, test_field.z)] = 1.0f;
        }
    }

    computePressure(&test_field, newPressure);

    std::cout << "Computed Pressure Values:" << std::endl;
    for (int j = 0; j < newPressure->y; ++j) {
        for (int i = 0; i < newPressure->x; ++i) {
            float pressure = newPressure->data[asIdx(i, j, 2, newPressure->x, newPressure->z)];
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
            test_field.data[asIdx(i, j, 0, test_field.x, test_field.z)] = 1.0f;
            test_field.data[asIdx(i, j, 1, test_field.x, test_field.z)] = 1.0f;
            test_field.data[asIdx(i, j, 2, test_field.x, test_field.z)] = 1.0f;
        }
    }

    subtractPressureGradient(&test_field, res_field);

    std::cout << "Computed Pressure Values:" << std::endl;
    for (int j = 0; j < res_field->y; ++j) {
        for (int i = 0; i < res_field->x; ++i) {
            std::cout << "[";
            for (int z = 0; z < res_field->z; ++z) {
                std::cout << res_field->data[asIdx(i, j, z, res_field->x, res_field->z)] << ", ";
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