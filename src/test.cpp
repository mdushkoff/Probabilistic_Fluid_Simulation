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

int main() {
    std::cout << "Running test for computePressure..." << std::endl;
    test_computePressure();
    std::cout << "Test completed successfully." << std::endl;
    return 0;
}