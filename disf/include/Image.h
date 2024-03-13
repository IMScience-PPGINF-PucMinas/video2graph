/**
* Image IO operations
* 
* @date September, 2019
*/
#ifndef IMAGE_H
#define IMAGE_H

#ifdef __cplusplus
extern "C" {
#endif

//=============================================================================
// Includes
//=============================================================================
#include "Utils.h"

//=============================================================================
// Structures
//=============================================================================
/**
* Simple image structure
*/
typedef struct
{
    int num_cols, num_rows, num_channels, num_pixels;
    // Each pixel whose index i < num_pixels contains num_channels of features,
    // which can be obtained through val[i]
    int **val;
} Image;

//=============================================================================
// Int Functions
//=============================================================================
/**
* Gets the maximum value within the image, considering the given channel. If 
* channel is -1, it will check all the channels.
*/
int getMaximumValue(Image *img, int channel);

/**
* Gets the minimum value within the image, considering the given channel. If 
* channel is -1, it will check all the channels.
*/
int getMinimumValue(Image *img, int channel);

/**
* Gets the normalization value (a.k.a. maximum value possible for a given bit depth)
*/
int getNormValue(Image *img);

//=============================================================================
// Image* Functions
//=============================================================================
/**
* Creates a black image with the number of rows, number of columns and number of
* channels given.
*/
Image **createVideo(int num_rows, int num_cols, int num_channels, int num_frames);

Image *createImage(int num_rows, int num_cols, int num_channels);

/**
* Loads an image pointed by the filepath given
*/
Image *loadImage(char* filepath);

/**
* Overlay the borders of the second image, into a copy of the first one. The RGB
* values are within the interval [0,1].
*/
Image *overlayBorders(Image *img, Image *border_img, float r, float g, float b);

//=============================================================================
// Void Functions
//=============================================================================
/**
* Deallocates the memory reserved for the image given in parameter
*/
 void freeVideo(Image **video, int num_frames);

void freeImage(Image **img);

/**
* Saves the image given as a PPM P6 file at the given filepath
*/
void writeImagePPM(Image *img, char* filepath);

/**
* Saves the image given as a PGM P5 file at the given filepath
*/
void writeImagePGM(Image *img, char* filepath);

#ifdef __cplusplus
}
#endif

#endif