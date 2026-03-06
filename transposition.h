#ifndef TRANSPOSITION_H
#define TRANSPOSITION_H
#include "render.h"

void transpose_central_diff_borders(TD * t);
void wall_central_diff_borders(TD * t);
void zero_central_diff_borders(TD * t);
void transpose_sobel_borders(TD * t);
void wall_sobel_borders(TD * t);
void zero_sobel_borders(TD * t);
void transpose_borders(TD * t);
void transpose_mesh_borders(TD * t);


#endif