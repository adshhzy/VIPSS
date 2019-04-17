//ppm.H
//
//This declares a few functions that will help save images
//as PPM files.
//
//Created: Thu Aug  5 12:23:04 PDT 2004
//Creator: Leon Barrett

#ifndef PPM_BY_LEON_RIGHT_HERE_DECLARE_H
#define PPM_BY_LEON_RIGHT_HERE_DECLARE_H

#include <iostream>

//Data should be arranged in row-major format (i.e. all of a row is together)
//with each pixel arranged as adjacent RGB.
//i.e. it should look like RGBRGBRGB...
void WritePPM( std::ostream & os, int width, int height, unsigned short * data );
void WritePPM( char * fileName, int width, int height, unsigned short * data );

#endif
