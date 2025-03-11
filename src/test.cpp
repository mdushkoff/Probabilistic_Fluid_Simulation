// System Includes
#include <chrono>
#include <iostream>
#include <math.h>
#include <sstream>
#include <string>
#include <utility>

// Local Includes
#include "../includes/fluid.hpp"

extern void test_computePressure();
extern void test_subtractPressureGradient();

int main() {
    std::cout << "Running test for computePressure..." << std::endl;
    // Test for Compute Pressure
    test_computePressure();

    // Test for Subtract Pressure Gradient
    // test_subtractPressureGradient();
    std::cout << "Test completed successfully." << std::endl;
    return 0;
}