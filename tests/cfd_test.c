/* cfd.h - v0.2 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) computational fluid dynamics library (CFD).

This Test class defines cases to verify that we don't break the excepted behaviours in the future upon changes.

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#include "../cfd.h"       /* Computational Fluid Dynamics API             */
#include "cfd_renderer.h" /* Utility to draw the LBM Grid to pixel buffer */
#include <stdlib.h>       /* malloc/free                                  */
#include <stdio.h>        /* File IO                                      */
#include "perf.h"         /* Simple performance profiler                  */

void cfd_test_write_combined_ppm(cfd_pixel_color *buffer, int width, int total_height, int frame)
{
  char filename[32];
  FILE *fp;

  sprintf(filename, "frame_%05d.ppm", frame);

  fp = fopen(filename, "wb");

  if (!fp)
  {
    return;
  }

  (void)fprintf(fp, "P6\n%d %d\n255\n", width, total_height);
  (void)fwrite(buffer, sizeof(cfd_pixel_color), (size_t)(width * total_height), fp);
  (void)fclose(fp);
}

void cfd_test_2d_simple_example(void)
{
  /* Setup the lattice boltzmann parameters */
  int xdim = 300;
  int ydim = 120;
  float viscosity = 0.02f;
  float omega = 1.0f / (3.0f * viscosity + 0.5f);
  float speed = 0.1f;

  int frame_count = 250;
  int frame = 0;
  int frame_steps = 20;

  unsigned long memory_grid_size = cfd_lbm_2d_grid_memory_size(xdim, ydim);
  unsigned long memory_size = memory_grid_size;

  /* Or nostdlib VirtualAlloc/mmap/... */
  void *memory = malloc(memory_size);

  cfd_lbm_2d_grid grid = {0};

  cfd_lbm_2d_init_grid(&grid, memory, xdim, ydim, omega);
  cfd_lbm_2d_init_fluid(&grid, speed);
  cfd_lbm_2d_init_barriers(&grid);
  cfd_lbm_2d_init_tracers(&grid);

  /* Run the simulation at 20 steps per frame */
  for (frame = 0; frame < frame_count; ++frame)
  {
    int step;
    int x, y;

    for (step = 0; step < frame_steps; ++step)
    {
      cfd_lbm_2d_collide_and_stream(&grid);
      cfd_lbm_2d_move_tracers(&grid);
    }

    /* The grid now contains the updated data for all cells */
    /* You can access them as follows                       */
    for (y = 0; y < grid.ydim; ++y)
    {
      for (x = 0; x < grid.xdim; ++x)
      {
        int idx = x + y * grid.xdim;

        if (!grid.barrier[idx])
        {
          float curl_value = cfd_lbm_2d_calculate_curl(&grid, x, y);
          (void)curl_value;
        }
      }
    }
  }

  free(memory);
}

void cfd_test_2d_full_visualization(void)
{
  /* --- CONTROLS --- */
  int pixel_per_square = 2; /* cfd_pixel_colors per grid site */
  float speed = 0.1f;       /* initial fluid speed            */
  float contrast = 0.0f;    /* contrast                       */
  float viscosity = 0.02f;  /* fluid viscosity                */
  float omega = 1.0f / (3.0f * viscosity + 0.5f);
  int frame_count = 250; /* Number of frames to generate */
  int frame_steps = 20;  /* simulation steps per frame   */
  int frame = 0;

  /* --- VISUALIZATIONS --- */
  int enable_tracer = 1;      /* enable_tracer      (0=off, 1=on) */
  int enable_flowline = 1;    /* enable_flowline    (0=off, 1=on) */
  int enable_force_arrow = 1; /* enable_force_arrow (0=off, 1=on) */
  int num_plots = 1;          /* Number of plot types from 0 to 6 */
  /*int plotSelect = 4; */    /* plotSelect (0:density, 1:x-vel, 2:y-vel, 3:speed, 4:curl, 5: pressure, 6: wall shear stress) */

  /* --- SETUP --- */
  int xdim = 600 / pixel_per_square;
  int ydim = 240 / pixel_per_square;

  /* --- SETUP FOR COMBINED PLOTS --- */
  int width_per_plot = xdim * pixel_per_square;
  int height_per_plot = ydim * pixel_per_square;
  int total_height = height_per_plot * num_plots;

  unsigned long memory_grid_size = cfd_lbm_2d_grid_memory_size(xdim, ydim);
  unsigned long memory_io_size = (unsigned long)(width_per_plot * total_height) * (unsigned long)sizeof(cfd_pixel_color);
  unsigned long memory_size = memory_grid_size + memory_io_size;
  void *memory = malloc(memory_size);
  cfd_pixel_color *memory_ppm = (cfd_pixel_color *)((char *)memory + memory_grid_size);

  cfd_lbm_2d_grid grid = {0};

  printf("[cfd][lbm][2d] mem. grid (MB): %10.4f\n", (double)memory_grid_size / 1024.0 / 1024.0);
  printf("[cfd][lbm][2d] mem. io   (MB): %10.4f\n", (double)memory_io_size / 1024.0 / 1024.0);
  printf("[cfd][lbm][2d]      viscosity: %10.4f\n", (double)viscosity);
  printf("[cfd][lbm][2d]          omega: %10.4f\n", (double)omega);
  printf("[cfd][lbm][2d]          speed: %10.4f\n", (double)speed);
  printf("[cfd][lbm][2d]              x: %10i\n", xdim);
  printf("[cfd][lbm][2d]              y: %10i\n", ydim);

  cfd_build_colormap();

  cfd_lbm_2d_init_grid(&grid, memory, xdim, ydim, omega);
  cfd_lbm_2d_init_fluid(&grid, speed);
  cfd_lbm_2d_init_barriers(&grid);

  if (enable_tracer)
  {
    cfd_lbm_2d_init_tracers(&grid);
  }

  /* --- SIMULATION LOOP --- */
  printf("[cfd][lbm][2d] Starting simulation...\n");
  for (frame = 0; frame < frame_count; ++frame)
  {
    int step;
    int plot_type;

    for (step = 0; step < frame_steps; ++step)
    {

      cfd_lbm_2d_collide_and_stream(&grid);

      if (enable_tracer)
      {
        cfd_lbm_2d_move_tracers(&grid);
      }
    }

    /* Loop through all plot types and draw them into the combined buffer
    PERF_PROFILE_WITH_NAME(
        for (plot_type = 0; plot_type < num_plots; ++plot_type) {
          int y_offset = plot_type * height_per_plot;
          cfd_lbm_2d_draw_single_plot(combined_buffer, &grid, width_per_plot, y_offset, plot_type, contrast, pixel_per_square, enable_tracer, enable_flowline, enable_force_arrow);
        },
        "lbm_2d_draw_plots");
        */

    (void)plot_type;

    cfd_lbm_2d_draw_single_plot(memory_ppm, &grid, width_per_plot, 0, 4, contrast, pixel_per_square, enable_tracer, enable_flowline, enable_force_arrow);
    cfd_test_write_combined_ppm(memory_ppm, width_per_plot, total_height, frame);
  }

  printf("[cfd][lbm][2d] Simulation finished.\n");

  free(memory);
}

void cfd_test_3d_simple_example(void)
{
  /* Setup the lattice boltzmann parameters */
  int xdim = 100;
  int ydim = 50;
  int zdim = 50;

  float viscosity = 0.02f;
  float omega = 1.0f / (3.0f * viscosity + 0.5f);
  float speed = 0.1f;

  int frame_count = 5;
  int frame = 0;
  int frame_steps = 20;

  unsigned long memory_grid_size = cfd_lbm_3d_grid_memory_size(xdim, ydim, zdim);
  unsigned long memory_size = memory_grid_size;

  /* Or nostdlib VirtualAlloc/mmap/... */
  void *memory = malloc(memory_size);

  cfd_lbm_3d_grid grid = {0};

  {
    printf("[cfd][lbm][3d] mem. grid (MB): %10.4f\n", (double)memory_grid_size / 1024.0 / 1024.0);
    printf("[cfd][lbm][3d]      viscosity: %10.4f\n", (double)viscosity);
    printf("[cfd][lbm][3d]          omega: %10.4f\n", (double)omega);
    printf("[cfd][lbm][3d]          speed: %10.4f\n", (double)speed);
    printf("[cfd][lbm][3d]              x: %10i\n", xdim);
    printf("[cfd][lbm][3d]              y: %10i\n", ydim);
    printf("[cfd][lbm][3d]              z: %10i\n", zdim);
  }

  cfd_lbm_3d_init_grid(&grid, memory, xdim, ydim, zdim, omega);
  cfd_lbm_3d_init_fluid(&grid, speed);
  cfd_lbm_3d_init_barriers(&grid);

  /* Run the simulation at 20 steps per frame */
  printf("[cfd][lbm][3d] Starting simulation...\n");

  for (frame = 0; frame < frame_count; ++frame)
  {
    int step;
    int x, y, z;

    PERF_PROFILE_WITH_NAME({
    for (step = 0; step < frame_steps; ++step)
    {
      cfd_lbm_3d_collide_and_stream(&grid);
    } }, "lbm_3d_step");

    /* The grid now contains the updated data for all cells */
    /* You can access them as follows                       */
    for (z = 0; z < grid.zdim; ++z)
    {
      for (y = 0; y < grid.ydim; ++y)
      {
        for (x = 0; x < grid.xdim; ++x)
        {
          int xdim_ydim = grid.xdim * grid.ydim;
          int idx = x + y * grid.xdim + z * xdim_ydim;

          if (!grid.barrier[idx])
          {
            float curl_value = cfd_lbm_3d_calculate_curl(&grid, x, y, z);
            (void)curl_value;
          }
        }
      }
    }
  }

  printf("[cfd][lbm][3d] Simulation finished.\n");

  free(memory);
}

int main(void)
{
  /* LBM D2Q9 */
  cfd_test_2d_simple_example();
  cfd_test_2d_full_visualization();

  /* LBM D3Q19 */
  cfd_test_3d_simple_example();

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
