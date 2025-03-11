#include <iostream>
#include <vector>
#include <cmath>
#include <iomanip>
#include <cstdlib>

#include "../../includes/fluid.hpp"

vp_field *init_vp_field(int width, int height, int depth, float init_value = 0.0f) {
    vp_field *vp = (vp_field *)malloc(sizeof(vp_field));
    if (!vp) {
        std::cerr << "Memory allocation failed for vp_field!" << std::endl;
        return nullptr;
    }

    vp->x = width;
    vp->y = height;
    vp->z = depth;
    vp->data = (float *)malloc(width * height * depth * sizeof(float));

    if (!vp->data) {
        std::cerr << "Memory allocation failed for data!" << std::endl;
        free(vp);
        return nullptr;
    }

    // Initialize data with a constant value
    for (int i = 0; i < width * height * depth; i++) {
        vp->data[i] = init_value;
    }

    return vp;
}

// Function to free vp_field memory
void free_vp_field(vp_field *vp) {
    if (vp) {
        free(vp->data);
        free(vp);
    }
}

void print_vp_field(const vp_field *vp) {
    int width = vp->x, height = vp->y;
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            std::cout << std::setw(6) << std::fixed << std::setprecision(2) << vp->data[i + j * width] << " ";
        }
        std::cout << std::endl;
    }
}

// Advect function declaration
void advect(vp_field *vp, vp_field *vp_out, double dt){
    // TODO: Perform advection
    int width = vp->x, height = vp->y, depth = vp->z;
    float *field = vp->data;
    float *new_field = vp_out->data;

    for (int j = 1; j < height - 1; j++) {
        for (int i = 1; i < width - 1; i++) { // Iterate through the inner part in the grid
            int cell_idx = (i + j * width); 

            double x_prev = i - dt * field[cell_idx] / width; // Compute the where the fluid came from at the previous time stamp // Divided by h so can convert velocity into grid units
            // Offset index to get velocity in y direction
            double y_prev = j - dt * field[cell_idx + (width * height * 1)] / height;

            x_prev = std::max(0.5, std::min(width - 1.5, x_prev)); // Get the center of cooridinates and prevent out-of-bounds errors
            y_prev = std::max(0.5, std::min(height - 1.5, y_prev));

            int i0 = int(x_prev), j0 = int(y_prev);
            int i1 = i0 + 1, j1 = j0 + 1; // Find the integer grid coordinates surrounding (x, y)
            double sx = x_prev - i0, sy = y_prev - j0; // Represent the fractional distances from (i0, j0)

                // Bilinear Interpolation
                // Try to understand why this is the weighted average value

            for (int vel_dir = 0; vel_dir < 2; vel_dir++) {
                int offset = (width * height * vel_dir);
                new_field[cell_idx + offset] = 
                    (1 - sx) * (1 - sy) * field[i0 + width * j0 + offset] +
                    sx * (1 - sy) * field[i1 + width * j0 + offset] +
                    (1 - sx) * sy * field[i0 + width * j1 + offset] +
                    sx * sy * field[i1 + width * j1 + offset];
            }
        }
    }
}


// Function to test the advect function
void test_advect() {
    int width = 10, height = 10, depth = 4;
    vp_field *vp = init_vp_field(width, height, depth);
    vp_field *vp_out = init_vp_field(width, height, depth);

    if (!vp || !vp_out) {
        std::cerr << "Failed to allocate memory for test" << std::endl;
        return;
    }

    // Set up a simple velocity field (e.g., constant rightward flow)
    for (int j = 0; j < height; j++) {
        for (int i = 0; i < width; i++) {
            int idx = i + j * width;
            vp->data[idx] = float(i) / float(width - 1);;                    // X velocity (constant rightward)
            vp->data[idx + width * height] = float(i) / float(height - 1);   // Y velocity (no vertical movement)
        }
    }

    double dt = 0.3; // Small time step
    advect(vp, vp_out, dt);
    print_vp_field(vp_out);

    // Free allocated memory
    free_vp_field(vp);
    free_vp_field(vp_out);
}

int main() {
    test_advect();
    return 0;
}
