/* cfd.h - v0.1 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) computational fluid dynamics library(CFD).

This Test class defines cases to verify that we don't break the excepted behaviours in the future upon changes.

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#include "../cfd.h"                /* Vector graphics generator                        */
#include "../cfd_platform_write.h" /* Optional: OS-Specific write file implementations */

#include "test.h" /* Simple Testing framework */

#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

CFD_API CFD_INLINE double cfd_sin(double x)
{
  double term = x;
  double result = x;
  double x2 = x * x;
  int i;

  for (i = 1; i <= 10; ++i)
  {
    term *= -x2 / ((2 * i) * (2 * i + 1));
    result += term;
  }

  return result;
}

CFD_API CFD_INLINE double cfd_cos(double x)
{
  double term = 1.0;
  double result = 1.0;
  double x2 = x * x;
  int i;

  for (i = 1; i <= 10; ++i)
  {
    term *= -x2 / ((2 * i - 1) * (2 * i));
    result += term;
  }

  return result;
}

CFD_API CFD_INLINE double cfd_atan(double z)
{
  int sign = 1;
  double result;

  if (z < 0.0)
  {
    sign = -1;
    z = -z;
  }

  if (z > 1.0)
  {
    result = 1.57079632679 - cfd_atan(1.0 / z); /* pi/2 - atan(1/z) */
  }
  else
  {
    double z2 = z * z;
    result = z * (1.0 - z2 * (1.0 / 3.0 - z2 * (1.0 / 5.0 - z2 * (1.0 / 7.0))));
  }

  return sign * result;
}

CFD_API CFD_INLINE double cfd_atan2(double y, double x)
{
  const double PI = 3.14159265359;

  if (x > 0)
    return cfd_atan(y / x);
  if (x < 0 && y >= 0)
    return cfd_atan(y / x) + PI;
  if (x < 0 && y < 0)
    return cfd_atan(y / x) - PI;
  if (x == 0 && y > 0)
    return PI / 2;
  if (x == 0 && y < 0)
    return -PI / 2;
  return 0.0;
}

CFD_API CFD_INLINE double cfd_pow(double base, double exp)
{
  /* Only handles positive base for now */
  int i;
  double result = 1.0;

  if (exp == 0.0)
    return 1.0;
  if (exp == 1.0)
    return base;

  if ((int)exp == exp)
  {
    /* Integer exponent */
    int e = (int)exp;
    for (i = 0; i < e; ++i)
      result *= base;
    return result;
  }

  return -1.0;
}

CFD_API CFD_INLINE int cfd_abs(int x)
{
  return (x < 0) ? -x : x;
}

/* Structure to hold a single cfd_pixel_color's color data */
typedef struct cfd_cfd_pixel_color_color
{
  unsigned char r;
  unsigned char g;
  unsigned char b;

} cfd_pixel_color;

void cfd_lbm_draw_line(cfd_pixel_color *buffer, int width, int height, int x0, int y0, int x1, int y1, cfd_pixel_color color)
{
  /* Simple Bresenham's line algorithm */
  int dx = cfd_abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = -cfd_abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = dx + dy, e2;

  for (;;)
  {
    if (x0 >= 0 && x0 < width && y0 >= 0 && y0 < height)
    {
      buffer[x0 + y0 * width] = color;
    }

    if (x0 == x1 && y0 == y1)
    {
      break;
    }

    e2 = 2 * err;

    if (e2 >= dy)
    {
      err += dy;
      x0 += sx;
    }

    if (e2 <= dx)
    {
      err += dx;
      y0 += sy;
    }
  }
}

void cfd_lbm_draw_tracers(cfd_pixel_color *buffer, cfd_lbm_grid *grid, int pxPerSquare)
{
  cfd_pixel_color color = {150, 150, 150};
  int width = grid->xdim * pxPerSquare;
  int height = grid->ydim * pxPerSquare;

  int t;

  for (t = 0; t < CFD_LBM_NUMBER_TRACERS; t++)
  {
    int canvasX = (int)((grid->tracerX[t] + 0.5) * pxPerSquare);
    int canvasY = height - 1 - (int)((grid->tracerY[t] + 0.5) * pxPerSquare);

    int i;

    for (i = -1; i <= 1; i++)
    {
      int j;

      for (j = -1; j <= 1; j++)
      {
        int px = canvasX + i;
        int py = canvasY + j;
        if (px >= 0 && px < width && py >= 0 && py < height)
        {
          buffer[px + py * width] = color;
        }
      }
    }
  }
}

void cfd_lbm_draw_flowlines(cfd_pixel_color *buffer, cfd_lbm_grid *grid, int pxPerSquare)
{
  int width = grid->xdim * pxPerSquare;
  int height = grid->ydim * pxPerSquare;
  double sitesPerFlowline = 10.0 / pxPerSquare;
  double y;

  for (y = sitesPerFlowline / 2; y < grid->ydim; y += sitesPerFlowline)
  {
    double x;

    for (x = sitesPerFlowline / 2; x < grid->xdim; x += sitesPerFlowline)
    {
      int ix = (int)x;
      int iy = (int)y;
      double thisUx = grid->ux[ix + iy * grid->xdim];
      double thisUy = grid->uy[ix + iy * grid->xdim];
      double speed = cfd_sqrt(thisUx * thisUx + thisUy * thisUy);

      if (speed > 0.0001)
      {
        int px = (int)((x + 0.5) * pxPerSquare);
        int py = height - 1 - (int)((y + 0.5) * pxPerSquare);
        /* The scaling factor determines the length of the lines. */
        double scale = 0.25 * pxPerSquare * 10.0 / speed;
        int x1 = (int)(px - thisUx * scale);
        int y1 = (int)(py + thisUy * scale);
        int x2 = (int)(px + thisUx * scale);
        int y2 = (int)(py - thisUy * scale);

        cfd_pixel_color color = {80, 80, 80};
        cfd_lbm_draw_line(buffer, width, height, x1, y1, x2, y2, color);
      }
    }
  }
}

void cfd_lbm_draw_force_arrow(cfd_pixel_color *buffer, cfd_lbm_grid *grid, int pxPerSquare)
{
  if (grid->barrierCount == 0)
  {
    return;
  }

  {
    double x = grid->barrierxSum / grid->barrierCount;
    double y = grid->barrierySum / grid->barrierCount;
    double Fx = grid->barrierFx;
    double Fy = grid->barrierFy;

    int width = grid->xdim * pxPerSquare;
    int height = grid->ydim * pxPerSquare;

    int canvasX = (int)((x + 0.5) * pxPerSquare);
    int canvasY = height - 1 - (int)((y + 0.5) * pxPerSquare);

    double magF = cfd_sqrt(Fx * Fx + Fy * Fy);

    double scale = 4.0 * magF * 100.0;
    int x1 = canvasX;
    int y1 = canvasY;
    int x2 = (int)(canvasX + Fx / magF * scale);
    int y2 = (int)(canvasY - Fy / magF * scale);

    cfd_pixel_color color = {0, 0, 0};

    /* Draw arrowhead */
    double angle = cfd_atan2(-Fy, Fx);
    double arrowAngle = 25.0 * 3.14159 / 180.0;
    double arrowLength = 0.2 * scale;
    int xA1 = (int)(x2 - arrowLength * cfd_cos(angle - arrowAngle));
    int yA1 = (int)(y2 - arrowLength * cfd_sin(angle - arrowAngle));
    int xA2 = (int)(x2 - arrowLength * cfd_cos(angle + arrowAngle));
    int yA2 = (int)(y2 - arrowLength * cfd_sin(angle + arrowAngle));

    if (magF < 1e-6)
    {
      return;
    }

    cfd_lbm_draw_line(buffer, width, height, x1, y1, x2, y2, color);
    cfd_lbm_draw_line(buffer, width, height, x2, y2, xA1, yA1, color);
    cfd_lbm_draw_line(buffer, width, height, x2, y2, xA2, yA2, color);
  }
}

void cfd_lbm_draw_canvas(cfd_lbm_grid *grid, const char *filename, int plotType, double contrast, int pxPerSquare, int tracerCheck, int flowlineCheck, int forceCheck)
{
  int width = grid->xdim * pxPerSquare;
  int height = grid->ydim * pxPerSquare;
  cfd_pixel_color *buffer = (cfd_pixel_color *)malloc((size_t)(width * height) * sizeof(cfd_pixel_color));

  double contrastFactor = cfd_pow(1.2, contrast);
  int y;

  /* Step 1: Draw the main fluid plot */
  if (plotType == 4)
  {
    cfd_lbm_compute_curl(grid);
  }

  for (y = 0; y < grid->ydim; y++)
  {
    int x;

    for (x = 0; x < grid->xdim; x++)
    {
      double value = 0;
      int cIndex;
      cfd_pixel_color color;

      if (grid->barrier[x + y * grid->xdim])
      {
        value = -1; /* Special value for barrier */
      }
      else
      {
        switch (plotType)
        {
        case 0:
          value = (grid->rho[x + y * grid->xdim] - 1.0) * 6.0;
          break;
        case 1:
          value = grid->ux[x + y * grid->xdim] * 2.0;
          break;
        case 2:
          value = grid->uy[x + y * grid->xdim] * 2.0;
          break;
        case 3:
        {
          double speed = cfd_sqrt(grid->ux[x + y * grid->xdim] * grid->ux[x + y * grid->xdim] + grid->uy[x + y * grid->xdim] * grid->uy[x + y * grid->xdim]);
          value = speed * 4.0;
          break;
        }
        case 4:
          value = grid->curl[x + y * grid->xdim] * 5.0;
          break;
        }
        value = value * contrastFactor + 0.5;
      }

      cIndex = (int)(400 * value);

      if (value == -1)
      {
        cIndex = -1;
      }
      if (cIndex < 0 && cIndex != -1)
      {
        cIndex = 0;
      }
      if (cIndex > 400)
      {
        cIndex = 400;
      }

      if (cIndex == -1)
      {
        color.r = 0;
        color.g = 0;
        color.b = 0;
      }
      else
      {
        if (cIndex < 50)
        {
          color.r = 0;
          color.g = 0;
          color.b = (unsigned char)(255 * (cIndex + 50) / 100.0);
        }
        else if (cIndex < 150)
        {
          color.r = 0;
          color.g = (unsigned char)(255 * (cIndex - 50) / 100.0);
          color.b = 255;
        }
        else if (cIndex < 250)
        {
          color.r = (unsigned char)(255 * (cIndex - 150) / 100.0);
          color.g = 255;
          color.b = 255 - color.r;
        }
        else if (cIndex < 350)
        {
          color.r = 255;
          color.g = (unsigned char)(255 * (350 - cIndex) / 100.0);
          color.b = 0;
        }
        else
        {
          color.r = (unsigned char)(255 * (450 - cIndex) / 100.0);
          color.g = 0;
          color.b = 0;
        }
      }

      /* Color the square in the buffer */
      {
        int flippedy = grid->ydim - y - 1;
        int py;

        for (py = flippedy * pxPerSquare; py < (flippedy + 1) * pxPerSquare; py++)
        {
          int px;

          for (px = x * pxPerSquare; px < (x + 1) * pxPerSquare; px++)
          {
            buffer[px + py * width] = color;
          }
        }
      }
    }
  }

  /* Step 2: Draw overlays on top of the buffer */
  if (flowlineCheck)
  {
    cfd_lbm_draw_flowlines(buffer, grid, pxPerSquare);
  }
  if (tracerCheck)
  {
    cfd_lbm_draw_tracers(buffer, grid, pxPerSquare);
  }
  if (forceCheck)
  {
    cfd_lbm_draw_force_arrow(buffer, grid, pxPerSquare);
  }

  /* Step 3: Write the final buffer to a PPM file */
  {
    FILE *fp = fopen(filename, "wb");

    if (!fp)
    {
      perror("Failed to open file for writing");
      free(buffer);
      return;
    }

    (void)fprintf(fp, "P6\n%d %d\n255\n", width, height);
    (void)fwrite(buffer, sizeof(cfd_pixel_color), (size_t)(width * height), fp);
    (void)fclose(fp);

    free(buffer);
  }
}

int main(void)
{
  /* --- CONTROLS --- */
  int pxPerSquare = 2;       /* cfd_pixel_colors per grid site */
  double speedSlider = 0.1;  /* initial fluid speed */
  double contrastSlider = 0; /* contrastSlider */
  int stepsSlider = 20;      /* simulation steps per frame */
  double viscSlider = 0.02;  /* fluid viscosity */
  int plotSelect = 3;        /* plotSelect (0:density, 1:x-vel, 2:y-vel, 3:speed, 4:curl) */
  int tracerCheck = 1;       /* tracerCheck (0=off, 1=on) */
  int flowlineCheck = 1;     /* flowlineCheck (0=off, 1=on) */
  int forceCheck = 1;        /* forceCheck (0=off, 1=on) */
  int frameCount = 1000;     /* Number of frames to generate */

  /* --- SETUP --- */
  int xdim = 600 / pxPerSquare;
  int ydim = 240 / pxPerSquare;

  int frame;

  cfd_lbm_grid grid;
  size_t gridSize = (size_t)(xdim * ydim) * sizeof(double);

  grid.xdim = xdim;
  grid.ydim = ydim;

  grid.n0 = (double *)malloc(gridSize);
  grid.nN = (double *)malloc(gridSize);
  grid.nS = (double *)malloc(gridSize);
  grid.nE = (double *)malloc(gridSize);
  grid.nW = (double *)malloc(gridSize);
  grid.nNE = (double *)malloc(gridSize);
  grid.nSE = (double *)malloc(gridSize);
  grid.nNW = (double *)malloc(gridSize);
  grid.nSW = (double *)malloc(gridSize);
  grid.rho = (double *)malloc(gridSize);
  grid.ux = (double *)malloc(gridSize);
  grid.uy = (double *)malloc(gridSize);
  grid.curl = (double *)malloc(gridSize);
  grid.barrier = (int *)calloc((size_t)(xdim * ydim), sizeof(int));

  cfd_lbm_init_fluid(&grid, speedSlider);
  cfd_lbm_init_barriers(&grid);
  if (tracerCheck)
  {
    cfd_lbm_init_tracers(&grid);
  }

  /* --- SIMULATION LOOP --- */
  printf("Starting simulation...\n");
  for (frame = 0; frame < frameCount; frame++)
  {
    int step;
    char filename[32];

    for (step = 0; step < stepsSlider; step++)
    {
      cfd_lbm_collide(&grid, viscSlider);
      cfd_lbm_stream(&grid);
      if (tracerCheck)
      {
        cfd_lbm_move_tracers(&grid);
      }
    }

    sprintf(filename, "frame_%05d.ppm", frame);
    cfd_lbm_draw_canvas(&grid, filename, plotSelect, contrastSlider, pxPerSquare, tracerCheck, flowlineCheck, forceCheck);
    printf("Generated %s\n", filename);
  }

  printf("Simulation finished.\n");

  /* --- CLEANUP --- */
  free(grid.n0);
  free(grid.nN);
  free(grid.nS);
  free(grid.nE);
  free(grid.nW);
  free(grid.nNE);
  free(grid.nSE);
  free(grid.nNW);
  free(grid.nSW);
  free(grid.rho);
  free(grid.ux);
  free(grid.uy);
  free(grid.curl);
  free(grid.barrier);

  return 0;
}

/*
   ------------------------------------------------------------------------------
   This software is available under 2 licenses -- choose whichever you prefer.
   ------------------------------------------------------------------------------
   ALTERNATIVE A - MIT License
   Copyright (c) 2025 nickscha
   Permission is hereby granted, free of charge, to any person obtaining a copy of
   this software and associated documentation files (the "Software"), to deal in
   the Software without restriction, including without limitation the rights to
   use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
   of the Software, and to permit persons to whom the Software is furnished to do
   so, subject to the following conditions:
   The above copyright notice and this permission notice shall be included in all
   copies or substantial portions of the Software.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
   LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
   OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
   SOFTWARE.
   ------------------------------------------------------------------------------
   ALTERNATIVE B - Public Domain (www.unlicense.org)
   This is free and unencumbered software released into the public domain.
   Anyone is free to copy, modify, publish, use, compile, sell, or distribute this
   software, either in source code form or as a compiled binary, for any purpose,
   commercial or non-commercial, and by any means.
   In jurisdictions that recognize copyright laws, the author or authors of this
   software dedicate any and all copyright interest in the software to the public
   domain. We make this dedication for the benefit of the public at large and to
   the detriment of our heirs and successors. We intend this dedication to be an
   overt act of relinquishment in perpetuity of all present and future rights to
   this software under copyright law.
   THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
   IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
   FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
   AUTHORS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
   ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION
   WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
   ------------------------------------------------------------------------------
*/
