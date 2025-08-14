/* cfd.h - v0.1 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) computational fluid dynamics library (CFD).

This Test class defines cases to verify that we don't break the excepted behaviours in the future upon changes.

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#include "../cfd.h"   /* Computational Fluid Dynamics API */
#include "cfd_math.h" /* Math functions to replace math.h */
#include <stdlib.h>   /* malloc/free                      */
#include <stdio.h>    /* File IO */
#include "perf.h"     /* Simple performance profiler      */

/* Structure to hold a single cfd_pixel_color's color data */
typedef struct cfd_pixel_color
{
  unsigned char r;
  unsigned char g;
  unsigned char b;

} cfd_pixel_color;

/* Build once at startup */
static cfd_pixel_color cfd_color_map[401]; /* 0..400; use index -1 for barrier separately */

CFD_API CFD_INLINE void cfd_build_colormap(void)
{
  int i;
  for (i = 0; i <= 400; ++i)
  {
    cfd_pixel_color color;
    if (i < 50)
    {
      color.r = 0;
      color.g = 0;
      color.b = (unsigned char)(((float)(255 * (i + 50))) / 100.0f);
    }
    else if (i < 150)
    {
      color.r = 0;
      color.g = (unsigned char)(((float)(255 * (i - 50))) / 100.0f);
      color.b = 255;
    }
    else if (i < 250)
    {
      color.r = (unsigned char)(((float)(255 * (i - 150))) / 100.0f);
      color.g = 255;
      color.b = 255 - color.r;
    }
    else if (i < 350)
    {
      color.r = 255;
      color.g = (unsigned char)(((float)(255 * (350 - i))) / 100.0f);
      color.b = 0;
    }
    else
    {
      color.r = (unsigned char)(((float)(255 * (450 - i))) / 100.0f);
      color.g = 0;
      color.b = 0;
    }
    cfd_color_map[i] = color;
  }
}

CFD_API CFD_INLINE void cfd_lbm_draw_line(cfd_pixel_color *buffer, int full_width, int y_offset, int plot_height, int x0, int y0, int x1, int y1, cfd_pixel_color color)
{
  /* Simple Bresenham's line algorithm */
  int dx = cfd_abs(x1 - x0), sx = x0 < x1 ? 1 : -1;
  int dy = -cfd_abs(y1 - y0), sy = y0 < y1 ? 1 : -1;
  int err = dx + dy, e2;

  int y_start = y_offset;
  int y_end = y_offset + plot_height;

  for (;;)
  {
    if (x0 >= 0 && x0 < full_width && y0 >= y_start && y0 < y_end)
    {
      buffer[x0 + y0 * full_width] = color;
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

CFD_API CFD_INLINE void cfd_lbm_draw_tracers(cfd_pixel_color *buffer, cfd_lbm_grid *grid, int full_width, int y_offset, int pxPerSquare)
{
  cfd_pixel_color color = {150, 150, 150};
  int plot_height = grid->ydim * pxPerSquare;

  int t;
  for (t = 0; t < CFD_LBM_NUMBER_TRACERS; ++t)
  {
    int canvasX = (int)((grid->tracerX[t] + 0.5f) * (float)pxPerSquare);
    int canvasY = y_offset + (plot_height - 1 - (int)((grid->tracerY[t] + 0.5f) * (float)pxPerSquare));

    int i;
    for (i = -1; i <= 1; ++i)
    {
      int j;
      for (j = -1; j <= 1; ++j)
      {
        int px = canvasX + i;
        int py = canvasY + j;
        if (px >= 0 && px < full_width && py >= y_offset && py < y_offset + plot_height)
        {
          buffer[px + py * full_width] = color;
        }
      }
    }
  }
}

CFD_API CFD_INLINE void cfd_lbm_draw_flowlines(cfd_pixel_color *buffer, cfd_lbm_grid *grid, int full_width, int y_offset, int pxPerSquare)
{
  int plot_height = grid->ydim * pxPerSquare;
  float sitesPerFlowline = 10.0f / (float)pxPerSquare;
  float y;

  for (y = sitesPerFlowline / 2.0f; y < (float)grid->ydim; y += sitesPerFlowline)
  {
    float x;
    for (x = sitesPerFlowline / 2.0f; x < (float)grid->xdim; x += sitesPerFlowline)
    {
      int ix = (int)x;
      int iy = (int)y;
      float thisUx = grid->ux[ix + iy * grid->xdim];
      float thisUy = grid->uy[ix + iy * grid->xdim];
      float speed = cfd_sqrtf(thisUx * thisUx + thisUy * thisUy);

      if (speed > 0.0001f)
      {
        int px = (int)((x + 0.5f) * (float)pxPerSquare);
        int py = y_offset + (plot_height - 1 - (int)((y + 0.5f) * (float)pxPerSquare));
        float scale = 0.25f * (float)pxPerSquare * 10.0f / speed;
        int x1 = (int)((float)px - thisUx * scale);
        int y1 = (int)((float)py + thisUy * scale);
        int x2 = (int)((float)px + thisUx * scale);
        int y2 = (int)((float)py - thisUy * scale);

        cfd_pixel_color color = {80, 80, 80};
        cfd_lbm_draw_line(buffer, full_width, y_offset, plot_height, x1, y1, x2, y2, color);
      }
    }
  }
}

CFD_API CFD_INLINE void cfd_lbm_draw_force_arrow(cfd_pixel_color *buffer, cfd_lbm_grid *grid, int full_width, int y_offset, int pxPerSquare)
{

  float x = grid->barrierxSum / (float)grid->barrierCount;
  float y = grid->barrierySum / (float)grid->barrierCount;
  float Fx = grid->barrierFx;
  float Fy = grid->barrierFy;
  int plot_height = grid->ydim * pxPerSquare;

  int canvasX = (int)((x + 0.5f) * (float)pxPerSquare);
  int canvasY = y_offset + (plot_height - 1 - (int)((y + 0.5f) * (float)pxPerSquare));

  float magF = cfd_sqrtf(Fx * Fx + Fy * Fy);

  float scale = 4.0f * magF * 100.0f;
  int x1 = canvasX;
  int y1 = canvasY;
  int x2 = (int)((float)canvasX + Fx / magF * scale);
  int y2 = (int)((float)canvasY - Fy / magF * scale);

  cfd_pixel_color color = {0, 0, 0};

  /* Draw arrowhead */
  float angle = cfd_atan2f(-Fy, Fx);
  float arrowAngle = 25.0f * 3.14159f / 180.0f;
  float arrowLength = 0.2f * scale;
  int xA1 = (int)((float)x2 - arrowLength * cfd_cosf(angle - arrowAngle));
  int yA1 = (int)((float)y2 - arrowLength * cfd_sinf(angle - arrowAngle));
  int xA2 = (int)((float)x2 - arrowLength * cfd_cosf(angle + arrowAngle));
  int yA2 = (int)((float)y2 - arrowLength * cfd_sinf(angle + arrowAngle));

  if (grid->barrierCount == 0)
  {
    return;
  }

  if (magF < 1e-6f)
  {
    return;
  }

  cfd_lbm_draw_line(buffer, full_width, y_offset, plot_height, x1, y1, x2, y2, color);
  cfd_lbm_draw_line(buffer, full_width, y_offset, plot_height, x2, y2, xA1, yA1, color);
  cfd_lbm_draw_line(buffer, full_width, y_offset, plot_height, x2, y2, xA2, yA2, color);
}

CFD_API CFD_INLINE void cfd_lbm_draw_single_plot(cfd_pixel_color *buffer, cfd_lbm_grid *grid, int width, int y_offset, int plotType, float contrast, int pxPerSquare, int tracerCheck, int flowlineCheck, int forceCheck)
{
  float contrastFactor = cfd_powf(1.2f, contrast);
  int y;

  /* Step 1: Draw the main fluid plot */
  for (y = 0; y < grid->ydim; ++y)
  {
    int x;
    for (x = 0; x < grid->xdim; ++x)
    {
      float value = 0.0f;
      int cIndex;
      cfd_pixel_color color;

      if (grid->barrier[x + y * grid->xdim])
      {
        value = -1.0f; /* Special value for barrier */
      }
      else
      {
        switch (plotType)
        {
        case 0:
          value = cfd_lbm_calculate_density(grid, x, y);
          break;
        case 1:
          value = cfd_lbm_calculate_velocity_x(grid, x, y);
          break;
        case 2:
          value = cfd_lbm_calculate_velocity_y(grid, x, y);
          break;
        case 3:
          value = cfd_lbm_calculate_speed(grid, x, y);
          break;
        case 4:
          value = cfd_lbm_calculate_curl(grid, x, y);
          break;
        case 5:
          value = cfd_lbm_calculate_pressure(grid, x, y);
          break;
        case 6:
          value = cfd_lbm_calculate_wall_shear_stress(grid, x, y);
          break;
        }
        value = value * contrastFactor + 0.5f;
      }

      cIndex = value == -1.0f ? -1 : (int)(400 * value);

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
        color = cfd_color_map[cIndex];
      }

      /* Color the square in the buffer */
      {
        int flippedy = grid->ydim - y - 1;
        int py;
        for (py = flippedy * pxPerSquare; py < (flippedy + 1) * pxPerSquare; ++py)
        {
          int px;
          for (px = x * pxPerSquare; px < (x + 1) * pxPerSquare; ++px)
          {
            buffer[px + (py + y_offset) * width] = color;
          }
        }
      }
    }
  }

  /* Step 2: Draw overlays on top of the buffer */
  if (flowlineCheck)
  {
    cfd_lbm_draw_flowlines(buffer, grid, width, y_offset, pxPerSquare);
  }
  if (tracerCheck)
  {
    cfd_lbm_draw_tracers(buffer, grid, width, y_offset, pxPerSquare);
  }
  if (forceCheck)
  {
    cfd_lbm_draw_force_arrow(buffer, grid, width, y_offset, pxPerSquare);
  }
}

CFD_API CFD_INLINE void cfd_write_combined_ppm(cfd_pixel_color *buffer, int width, int total_height, int frame)
{
  char filename[32];
  FILE *fp;

  sprintf(filename, "frame_%05d.ppm", frame);

  fp = fopen(filename, "wb");
  if (!fp)
  {
    free(buffer);
    return;
  }
  (void)fprintf(fp, "P6\n%d %d\n255\n", width, total_height);
  (void)fwrite(buffer, sizeof(cfd_pixel_color), (size_t)(width * total_height), fp);
  (void)fclose(fp);
}

int main(void)
{
  /* --- CONTROLS --- */
  int pxPerSquare = 2;         /* cfd_pixel_colors per grid site */
  float speedSlider = 0.1f;    /* initial fluid speed */
  float contrastSlider = 0.0f; /* contrastSlider */
  int stepsSlider = 20;        /* simulation steps per frame */
  float viscSlider = 0.02f;    /* fluid viscosity */
  /*int plotSelect = 4; */     /* plotSelect (0:density, 1:x-vel, 2:y-vel, 3:speed, 4:curl, 5: pressure, 6: wall shear stress) */
  int num_plots = 7;           /* Number of plot types from 0 to 6 */
  int tracerCheck = 1;         /* tracerCheck (0=off, 1=on) */
  int flowlineCheck = 1;       /* flowlineCheck (0=off, 1=on) */
  int forceCheck = 1;          /* forceCheck (0=off, 1=on) */
  int frameCount = 250;        /* Number of frames to generate */
  int frame = 0;

  /* --- SETUP --- */
  int xdim = 600 / pxPerSquare;
  int ydim = 240 / pxPerSquare;

  void *memory = malloc(cfd_lbm_grid_memory_size(xdim, ydim));

  cfd_lbm_grid grid = {0};

  /* --- SETUP FOR COMBINED PLOTS --- */
  int width_per_plot = xdim * pxPerSquare;
  int height_per_plot = ydim * pxPerSquare;
  int total_height = height_per_plot * num_plots;
  cfd_pixel_color *combined_buffer = (cfd_pixel_color *)malloc((size_t)(width_per_plot * total_height) * sizeof(cfd_pixel_color));

  cfd_build_colormap();

  cfd_lbm_init_grid(&grid, memory, xdim, ydim);
  cfd_lbm_init_fluid(&grid, speedSlider);
  cfd_lbm_init_barriers(&grid);

  if (tracerCheck)
  {
    cfd_lbm_init_tracers(&grid);
  }

  /* --- SIMULATION LOOP --- */
  printf("Starting simulation...\n");
  for (frame = 0; frame < frameCount; ++frame)
  {
    int step;
    int plot_type;

    PERF_PROFILE_WITH_NAME({
    for (step = 0; step < stepsSlider; ++step)
    {
      cfd_lbm_collide(&grid, viscSlider);
      cfd_lbm_stream(&grid);
      if (tracerCheck)
      {
        cfd_lbm_move_tracers(&grid);
      }
    } }, "lbm_frame_step");

    /* Loop through all plot types and draw them into the combined buffer */
    PERF_PROFILE_WITH_NAME(
        for (plot_type = 0; plot_type < num_plots; ++plot_type) {
          int y_offset = plot_type * height_per_plot;
          cfd_lbm_draw_single_plot(combined_buffer, &grid, width_per_plot, y_offset, plot_type, contrastSlider, pxPerSquare, tracerCheck, flowlineCheck, forceCheck);
        },
        "lbm_draw_plots");

    PERF_PROFILE_WITH_NAME(
        cfd_write_combined_ppm(combined_buffer, width_per_plot, total_height, frame),
        "lbm_write_ppm");
  }

  printf("Simulation finished.\n");

  /* --- CLEANUP --- */
  free(memory);
  free(combined_buffer);

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
