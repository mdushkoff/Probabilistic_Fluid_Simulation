/*
 * Various utility functions for doing operations on
 * PNG files.
 */

#ifndef UTILS_HPP_
#define UTILS_HPP_

// Includes
#include <png.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

// Conditional Includes
#ifdef USE_CUDA
#include <cuda_runtime.h>
#endif // USE_CUDA

/*
 * Read a png from a file into a floating point array.
 * 
 * Referenced from: https://github.com/pnggroup/libpng/blob/master/contrib/examples/pngtopng.c
 *
 * Inputs:
 *     image - The reference to the png_image structure
 *     fn - The filename to read from
 *     x - The output floating point array
 * Outputs:
 *     status - Whether or not the read was successful (0) or not (1)
 */
int read_png_to_array(png_image *image, char *fn, float **x){
    // Validate that the image a valid png file
    FILE *fp = fopen(fn, "rb");
    if (!fp){
        // Unsuccessful read
        return 1;
    }
    unsigned char sig[8];
    fread(sig, 1, 8, fp);  // Try to read signature bytes
    if (!png_check_sig(sig, 8)){
        // Unsuccessful signature check
        fclose(fp);
        return 1;
    }
    fclose(fp);

    /* Initialize the png_image structure */
    memset(image, 0, (sizeof image));
    image->version = PNG_IMAGE_VERSION;

    // Check if the file could be read from
    if (png_image_begin_read_from_file(image, fn)){
        // Declare a byte buffer for storing information
        png_bytep buffer;

        // Set the format
        image->format = PNG_FORMAT_RGBA;

        // Allocate memory to the buffer
        buffer = (png_bytep)malloc(PNG_IMAGE_SIZE(*image));

        // Check if the buffer allocated properly
        if (buffer != NULL){
            // Read PNG data to the buffer
            if (png_image_finish_read(image, NULL, buffer, 0, NULL)){
                // Allocate memory in the output pointer
                int n_elements = (image->height * image->width)*4;
#ifdef USE_CUDA
                cudaError_t cuda_status = cudaMallocHost((void**)x, sizeof(float)*n_elements, cudaHostAllocDefault); // Pinned host memory
                if (cuda_status != cudaSuccess){
                    // Failed to allocate pinned host memory
                    free(buffer);
                    return 1;
                }
                /*(*x) = (float*)malloc(sizeof(float)*n_elements);*/
#else
                (*x) = (float*)malloc(sizeof(float)*n_elements);
#endif // USE_CUDA

                // Convert data to floating point
                for (int i=0; i<n_elements; i++){
                    (*x)[i] = (float)buffer[i] / 255.0;
                }
            }
            else {
                // Failed to finish reading
                free(buffer);
                return 1;
            }
        }
        else {
            // Unsuccessful buffer allocation
            png_image_free(image);  // The only point where this has to be cleaned up due to lack of memory
            return 1;
        }

        // Clean up memory
        free(buffer);
    }
    else {
        // Unsuccessful png read from file
        return 1;
    }

    // Successful read
    return 0;
}

/*
 * Write a png to a file from a floating point array.
 *
 * Inputs:
 *     image - The reference to the png_image structure
 *     fn - The filename to write to
 *     x - The input floating point array
 * Outputs:
 *     status - Whether or not the write was successful (0) or not (1)
 */
int write_png_from_array(png_image *image, const char *fn, float **x){
    // Allocate space in the buffer
    png_bytep buffer;
    buffer = (png_bytep)malloc(PNG_IMAGE_SIZE(*image));

    // Check if the buffer allocated properly
    if (buffer != NULL){
        // Convert floating point data to byte
        int n_elements = (image->height * image->width)*4;
        for (int i=0; i<n_elements; i++){
            buffer[i] = (png_byte)((*x)[i] * 255.0);
        }

        // Write the buffer to a file
        if (png_image_write_to_file(image, fn, 0, buffer, 0, NULL)){
            // Successful write
            free(buffer);
            return 0;
        }
        else{
            // Unsuccessful write
            free(buffer);
            return 1;
        }
    }
    else {
        // Unsuccessful buffer allocation
        png_image_free(image);  // The only point where this has to be cleaned up due to lack of memory
        return 1;
    }
}

#endif // UTILS_HPP_