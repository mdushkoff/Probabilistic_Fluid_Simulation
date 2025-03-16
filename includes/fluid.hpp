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
__global__ void advect(float *vp, float *vp_out, float dt, int vx, int vy, int vz);
#else
void advect(vp_field *vp, vp_field *vp_out, float dt);
#endif // USE_CUDA

/*
 * Inputs:
 *     image - The color image to advect
 *     out - The resulting output image
 *     vp - A 3D array representing the current velocity/pressure values
 *     dt - Change in time (time step)
 */
#ifdef USE_CUDA
__global__ void advect_color(float *image, float *out, float *vp, float dt, int ix, int iy, int iz, int vx, int vy, int vz);
#else
void advect_color(vp_field *image, vp_field *out, vp_field *vp, float dt);
#endif // USE_CUDA

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *     vp_out - Result for modified vp
 *     viscosity - The viscosity of the fluid
 *     dt - Change in time (time step)
 */
#ifdef USE_CUDA
__global__ void diffuse(float *vp, float *vp_out, float viscosity, float dt, int vx, int vy, int vz);
#else
void diffuse(vp_field *vp, vp_field *vp_out, float viscosity, float dt);
#endif // USE_CUDA

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *         which is updated after running a single step of simulation
 *     forces - A 3D array of force values for each pixel in the grid
 */
#ifdef USE_CUDA
__global__ void addForces(float *vp, float *forces, int vx, int vy, int vz);
#else
void addForces(vp_field *vp, float *forces);
#endif // USE_CUDA

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *     vp_out - Result for modified vp
 */
#ifdef USE_CUDA
__global__ void computePressure(float *vp, float *vp_out, float dt, int vx, int vy, int vz);
#else
void computePressure(vp_field *vp, vp_field *vp_out, float dt);
#endif // USE_CUDA

/*
 * Inputs:
 *     vp - A 3D array representing the current velocity/pressure values
 *     vp_out - Result for modified vp
 */
#ifdef USE_CUDA
__global__ void subtractPressureGradient(float *vp, float *vp_out, float dt, int vx, int vy, int vz);
#else
void subtractPressureGradient(vp_field *vp, vp_field *vp_out, float dt);
#endif // USE_CUDA

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
#ifdef USE_CUDA
void simulate_fluid_step(float **vp, float **tmp, float dt, float viscosity, int vx, int vy, int vz);
#else
void simulate_fluid_step(vp_field *vp, vp_field *tmp, float dt, float viscosity);
#endif // USE_CUDA

/*
 *
 */
#ifdef USE_CUDA
void advect_color_step(float **image, float **itmp, float **vp, float dt, int ix, int iy, int iz, int vx, int vy, int vz);
#else
void advect_color_step(vp_field *image, vp_field *itmp, vp_field *vp, float dt);
#endif // USE_CUDA

#endif // FLUID_HPP_