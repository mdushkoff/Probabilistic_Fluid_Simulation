/*
 * Definitions for running fluid simulation.
 */


// Header guard
#ifndef FLUID_HPP_
#define FLUID_HPP_

// Definitions
#define NUM_JACOBI_ITERS (30)  // Total number of Jacobi iterations

/*
 * Definition of a velocity-pressure field where
 * data is stored in a 3D array of size X,Y,Z.
 */
typedef struct {
    int x;
    int y;
    int z;
    float *data; // The data array arranged in PNG order
} vp_field;

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *     vp_out - Result for modified vp
 */
#ifdef USE_CUDA
__global__
#endif // USE_CUDA
void advect(vp_field *vp, vp_field *vp_out);

/*
 * Inputs:
 *     image - The color image to advect
 *     out - The resulting output image
 *     vp - A 3D array representing the current velocity/pressure values
 *     dt - Change in time (time step)
 */
#ifdef USE_CUDA
__global__
#endif // USE_CUDA
void advect_color(vp_field *image, vp_field *out, vp_field *vp, float dt);

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *     vp_out - Result for modified vp
 *     viscosity - The viscosity of the fluid
 *     dt - Change in time (time step)
 */
#ifdef USE_CUDA
__global__
#endif // USE_CUDA
void diffuse(vp_field *vp, vp_field *vp_out, float viscosity, float dt);

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *         which is updated after running a single step of simulation
 *     forces - A 3D array of force values for each pixel in the grid
 */
#ifdef USE_CUDA
__global__
#endif // USE_CUDA
void addForces(vp_field *vp, float *forces);

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *     vp_out - Result for modified vp
 */
#ifdef USE_CUDA
__global__
#endif // USE_CUDA
void computePressure(vp_field *vp, vp_field *vp_out, float dt);

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *     vp_out - Result for modified vp
 */
#ifdef USE_CUDA
__global__
#endif // USE_CUDA
void subtractPressureGradient(vp_field *vp, vp_field *vp_out, float dt);

/*
 * Simulate the current fluid domain for a single
 * timestep and output the results.
 *
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *         which is updated after running a single step of simulation
 *     tmp - The temporary buffer for storing results
 *     dt - The simulation time step resolution
 *     viscosity - The viscosity of the fluid
 */
void simulate_fluid_step(vp_field *vp, vp_field *tmp, float dt, float viscosity);

/*
 *
 */
void advect_color_step(vp_field *image, vp_field *itmp, vp_field *vp, float dt);

#endif // FLUID_HPP_