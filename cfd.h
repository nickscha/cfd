/* cfd.h - v0.1 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) computational fluid dynamics library (CFD).

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#ifndef CFD_H
#define CFD_H

/* #############################################################################
 * # COMPILER SETTINGS
 * #############################################################################
 */
/* Check if using C99 or later (inline is supported) */
#if __STDC_VERSION__ >= 199901L
#define CFD_INLINE inline
#define CFD_API static
#elif defined(__GNUC__) || defined(__clang__)
#define CFD_INLINE __inline__
#define CFD_API static
#elif defined(_MSC_VER)
#define CFD_INLINE __inline
#define CFD_API static
#else
#define CFD_INLINE
#define CFD_API static
#endif

/* #############################################################################
 * # MATH Functions
 * #############################################################################
 */
#define CFD_RAND_MAX 32767

static unsigned long cfd_rand_next = 1234;

CFD_API CFD_INLINE int cfd_rand(void)
{
  cfd_rand_next = cfd_rand_next * 1103515245 + 12345;
  return (int)((cfd_rand_next >> 16) & CFD_RAND_MAX);
}

CFD_API CFD_INLINE float cfd_sqrtf(float x)
{
  float guess, prev;
  int i;

  if (x < 0.0f)
  {
    return -1.0f; /* Or handle as needed */
  }

  if (x == 0.0f)
  {
    return 0.0f;
  }

  guess = x > 1.0f ? x : 1.0f; /* Initial guess */
  for (i = 0; i < 20; ++i)
  {
    prev = guess;
    guess = 0.5f * (guess + x / guess);

    /* Optional convergence check */
    if (guess == prev)
    {
      break;
    }
  }

  return guess;
}

/* #############################################################################
 * # LBM D2Q9 Model
 * #############################################################################
 */
#define CFD_LBM_2D_FOUR_NINTHS (4.0f / 9.0f)
#define CFD_LBM_2D_ONE_NINTH (1.0f / 9.0f)
#define CFD_LBM_2D_ONE_36TH (1.0f / 36.0f)
#define CFD_LBM_2D_NUMBER_TRACERS 144

/* Structure to hold all the grid data for the LBM simulation */
typedef struct cfd_lbm_2d_grid
{
  /* Grid size */
  int xdim;
  int ydim;

  /* Force data */
  float barrierFx, barrierFy;
  int barrierCount;
  float barrierxSum, barrierySum;

  /* Tracer data */
  float tracerX[CFD_LBM_2D_NUMBER_TRACERS];
  float tracerY[CFD_LBM_2D_NUMBER_TRACERS];

  /* Particle distributions */
  float *nC, *nN, *nS, *nE, *nW, *nNE, *nSE, *nNW, *nSW;

  float *rho, *ux, *uy;

  unsigned char *barrier;

} cfd_lbm_2d_grid;

CFD_API CFD_INLINE unsigned long cfd_lbm_2d_grid_memory_size(int xdim, int ydim)
{
  unsigned long nSites = (unsigned long)(xdim * ydim);

  return (
      nSites * sizeof(float) * 12      /* 12 float arrays        */
      + nSites * sizeof(unsigned char) /* 1 byte array (barrier) */
  );
}

CFD_API CFD_INLINE void cfd_lbm_2d_init_grid(cfd_lbm_2d_grid *grid, void *memory, int xdim, int ydim)
{
  char *ptr = (char *)memory;

  unsigned long nSites = (unsigned long)(xdim * ydim);
  unsigned long sizeFloat = sizeof(float);
  unsigned long dist_size = nSites * sizeFloat;

  grid->xdim = xdim;
  grid->ydim = ydim;
  grid->nC = (float *)ptr;
  ptr += dist_size;
  grid->nN = (float *)ptr;
  ptr += dist_size;
  grid->nS = (float *)ptr;
  ptr += dist_size;
  grid->nE = (float *)ptr;
  ptr += dist_size;
  grid->nW = (float *)ptr;
  ptr += dist_size;
  grid->nNE = (float *)ptr;
  ptr += dist_size;
  grid->nSE = (float *)ptr;
  ptr += dist_size;
  grid->nNW = (float *)ptr;
  ptr += dist_size;
  grid->nSW = (float *)ptr;
  ptr += dist_size;
  grid->rho = (float *)ptr;
  ptr += dist_size;
  grid->ux = (float *)ptr;
  ptr += dist_size;
  grid->uy = (float *)ptr;
  ptr += dist_size;
  grid->barrier = (unsigned char *)ptr;
  ptr += nSites * sizeof(unsigned char);
}

/* This function places barriers on the grid. */
CFD_API CFD_INLINE void cfd_lbm_2d_init_barriers(cfd_lbm_2d_grid *grid)
{
  int barrierSize = 8;
  int x = (int)(grid->xdim / 3);
  int y;

  for (y = (grid->ydim / 2) - barrierSize; y <= (grid->ydim / 2) + barrierSize; ++y)
  {
    if (x >= 0 && x < grid->xdim && y >= 0 && y < grid->ydim)
    {
      grid->barrier[x + y * grid->xdim] = 1;
    }
  }
}

/* Set all densities in a cell to their equilibrium values for a given velocity and density. */
CFD_API CFD_INLINE void cfd_lbm_2d_init_equilibrium(cfd_lbm_2d_grid *grid, int x, int y, float newux, float newuy, float newrho)
{
  int i = x + y * grid->xdim;
  float ux3 = 3.0f * newux;
  float uy3 = 3.0f * newuy;
  float ux2 = newux * newux;
  float uy2 = newuy * newuy;
  float uxuy2 = 2.0f * newux * newuy;
  float u2 = ux2 + uy2;
  float u215 = 1.5f * u2;
  float one9th_newrho = CFD_LBM_2D_ONE_NINTH * newrho;
  float one36th_newrho = CFD_LBM_2D_ONE_36TH * newrho;

  grid->nC[i] = CFD_LBM_2D_FOUR_NINTHS * newrho * (1.0f - u215);
  grid->nE[i] = one9th_newrho * (1.0f + ux3 + 4.5f * ux2 - u215);
  grid->nW[i] = one9th_newrho * (1.0f - ux3 + 4.5f * ux2 - u215);
  grid->nN[i] = one9th_newrho * (1.0f + uy3 + 4.5f * uy2 - u215);
  grid->nS[i] = one9th_newrho * (1.0f - uy3 + 4.5f * uy2 - u215);
  grid->nNE[i] = one36th_newrho * (1.0f + ux3 + uy3 + 4.5f * (u2 + uxuy2) - u215);
  grid->nSE[i] = one36th_newrho * (1.0f + ux3 - uy3 + 4.5f * (u2 - uxuy2) - u215);
  grid->nNW[i] = one36th_newrho * (1.0f - ux3 + uy3 + 4.5f * (u2 - uxuy2) - u215);
  grid->nSW[i] = one36th_newrho * (1.0f - ux3 - uy3 + 4.5f * (u2 + uxuy2) - u215);
  grid->rho[i] = newrho;
  grid->ux[i] = newux;
  grid->uy[i] = newuy;
}

/* Initialize the fluid to a steady rightward flow. */
CFD_API CFD_INLINE void cfd_lbm_2d_init_fluid(cfd_lbm_2d_grid *grid, float u0)
{
  int y;
  for (y = 0; y < grid->ydim; ++y)
  {
    int x;
    for (x = 0; x < grid->xdim; ++x)
    {
      cfd_lbm_2d_init_equilibrium(grid, x, y, u0, 0.0f, 1.0f);
    }
  }
}

/* Place tracers in a grid formation. */
CFD_API CFD_INLINE void cfd_lbm_2d_init_tracers(cfd_lbm_2d_grid *grid)
{
  int nRows = (int)(cfd_sqrtf((float)CFD_LBM_2D_NUMBER_TRACERS) + 0.5f);

  float dx;
  float dy;
  float nextX;
  float nextY;

  int t;

  if (nRows <= 0)
  {
    nRows = 1;
  }

  dx = (float)grid->xdim / (float)nRows;
  dy = (float)grid->ydim / (float)nRows;

  nextX = dx * 0.5f;
  nextY = dy * 0.5f;

  for (t = 0; t < CFD_LBM_2D_NUMBER_TRACERS; ++t)
  {
    grid->tracerX[t] = nextX;
    grid->tracerY[t] = nextY;

    nextX += dx;

    if (nextX >= (float)grid->xdim)
    {
      nextX = dx * 0.5f;
      nextY += dy;
    }
  }
}

/* Collide particles within each cell. */
CFD_API CFD_INLINE void cfd_lbm_2d_collide(cfd_lbm_2d_grid *grid, float omega)
{
  int y;
  int x;

  for (y = 1; y < grid->ydim - 1; ++y)
  {
    for (x = 1; x < grid->xdim - 1; ++x)
    {
      int i = x + y * grid->xdim;
      float thisrho = grid->nC[i] + grid->nN[i] + grid->nS[i] + grid->nE[i] + grid->nW[i] + grid->nNW[i] + grid->nNE[i] + grid->nSW[i] + grid->nSE[i];
      float invRho = 1.0f / thisrho;
      float thisux = (grid->nE[i] + grid->nNE[i] + grid->nSE[i] - grid->nW[i] - grid->nNW[i] - grid->nSW[i]) * invRho;
      float thisuy = (grid->nN[i] + grid->nNE[i] + grid->nNW[i] - grid->nS[i] - grid->nSE[i] - grid->nSW[i]) * invRho;
      float one9thrho = CFD_LBM_2D_ONE_NINTH * thisrho;
      float one36thrho = CFD_LBM_2D_ONE_36TH * thisrho;
      float ux3 = 3.0f * thisux;
      float uy3 = 3.0f * thisuy;
      float ux2 = thisux * thisux;
      float uy2 = thisuy * thisuy;
      float uxuy2 = 2.0f * thisux * thisuy;
      float u2 = ux2 + uy2;
      float u215 = 1.5f * u2;

      grid->rho[i] = thisrho;
      grid->ux[i] = thisux;
      grid->uy[i] = thisuy;

      grid->nC[i] += omega * (CFD_LBM_2D_FOUR_NINTHS * thisrho * (1.0f - u215) - grid->nC[i]);
      grid->nE[i] += omega * (one9thrho * (1.0f + ux3 + 4.5f * ux2 - u215) - grid->nE[i]);
      grid->nW[i] += omega * (one9thrho * (1.0f - ux3 + 4.5f * ux2 - u215) - grid->nW[i]);
      grid->nN[i] += omega * (one9thrho * (1.0f + uy3 + 4.5f * uy2 - u215) - grid->nN[i]);
      grid->nS[i] += omega * (one9thrho * (1.0f - uy3 + 4.5f * uy2 - u215) - grid->nS[i]);
      grid->nNE[i] += omega * (one36thrho * (1.0f + ux3 + uy3 + 4.5f * (u2 + uxuy2) - u215) - grid->nNE[i]);
      grid->nSE[i] += omega * (one36thrho * (1.0f + ux3 - uy3 + 4.5f * (u2 - uxuy2) - u215) - grid->nSE[i]);
      grid->nNW[i] += omega * (one36thrho * (1.0f - ux3 + uy3 + 4.5f * (u2 - uxuy2) - u215) - grid->nNW[i]);
      grid->nSW[i] += omega * (one36thrho * (1.0f - ux3 - uy3 + 4.5f * (u2 + uxuy2) - u215) - grid->nSW[i]);
    }
  }
}

/* Move particles along their directions of motion. */
CFD_API CFD_INLINE void cfd_lbm_2d_stream(cfd_lbm_2d_grid *grid)
{
  int y;
  int x;

  grid->barrierCount = 0;
  grid->barrierxSum = 0.0f;
  grid->barrierySum = 0.0f;
  grid->barrierFx = 0.0f;
  grid->barrierFy = 0.0f;

  for (y = grid->ydim - 2; y > 0; --y)
  {
    for (x = 1; x < grid->xdim - 1; ++x)
    {
      grid->nN[x + y * grid->xdim] = grid->nN[x + (y - 1) * grid->xdim];
      grid->nNW[x + y * grid->xdim] = grid->nNW[x + 1 + (y - 1) * grid->xdim];
    }
  }

  for (y = grid->ydim - 2; y > 0; --y)
  {
    for (x = grid->xdim - 2; x > 0; --x)
    {
      grid->nE[x + y * grid->xdim] = grid->nE[x - 1 + y * grid->xdim];
      grid->nNE[x + y * grid->xdim] = grid->nNE[x - 1 + (y - 1) * grid->xdim];
    }
  }

  for (y = 1; y < grid->ydim - 1; ++y)
  {
    for (x = grid->xdim - 2; x > 0; --x)
    {
      grid->nS[x + y * grid->xdim] = grid->nS[x + (y + 1) * grid->xdim];
      grid->nSE[x + y * grid->xdim] = grid->nSE[x - 1 + (y + 1) * grid->xdim];
    }
  }

  for (y = 1; y < grid->ydim - 1; ++y)
  {
    for (x = 1; x < grid->xdim - 1; ++x)
    {
      grid->nW[x + y * grid->xdim] = grid->nW[x + 1 + y * grid->xdim];
      grid->nSW[x + y * grid->xdim] = grid->nSW[x + 1 + (y + 1) * grid->xdim];
    }
  }

  /* Handle bounce-back from barriers and calculate force. */
  for (y = 1; y < grid->ydim - 1; ++y)
  {
    for (x = 1; x < grid->xdim - 1; ++x)
    {
      int index = x + y * grid->xdim;

      if (grid->barrier[index])
      {
        grid->nE[x + 1 + y * grid->xdim] = grid->nW[index];
        grid->nW[x - 1 + y * grid->xdim] = grid->nE[index];
        grid->nN[x + (y + 1) * grid->xdim] = grid->nS[index];
        grid->nS[x + (y - 1) * grid->xdim] = grid->nN[index];
        grid->nNE[x + 1 + (y + 1) * grid->xdim] = grid->nSW[index];
        grid->nNW[x - 1 + (y + 1) * grid->xdim] = grid->nSE[index];
        grid->nSE[x + 1 + (y - 1) * grid->xdim] = grid->nNW[index];
        grid->nSW[x - 1 + (y - 1) * grid->xdim] = grid->nNE[index];

        /* Sum forces on barrier sites */
        grid->barrierCount++;
        grid->barrierxSum += (float)x;
        grid->barrierySum += (float)y;
        grid->barrierFx += grid->nE[index] + grid->nNE[index] + grid->nSE[index] - grid->nW[index] - grid->nNW[index] - grid->nSW[index];
        grid->barrierFy += grid->nN[index] + grid->nNE[index] + grid->nNW[index] - grid->nS[index] - grid->nSE[index] - grid->nSW[index];
      }
    }
  }
}

/* Move tracer particles according to the fluid velocity. */
CFD_API CFD_INLINE void cfd_lbm_2d_move_tracers(cfd_lbm_2d_grid *grid)
{
  int t;
  for (t = 0; t < CFD_LBM_2D_NUMBER_TRACERS; ++t)
  {
    int roundedX = (int)(grid->tracerX[t] + 0.5f);
    int roundedY = (int)(grid->tracerY[t] + 0.5f);

    int index;

    if (roundedX < 0 || roundedX >= grid->xdim || roundedY < 0 || roundedY >= grid->ydim)
    {
      continue;
    }

    index = roundedX + roundedY * grid->xdim;

    grid->tracerX[t] += grid->ux[index];
    grid->tracerY[t] += grid->uy[index];

    if (grid->tracerX[t] > (float)(grid->xdim - 1))
    {
      grid->tracerX[t] = 0.0f;
      grid->tracerY[t] = (float)cfd_rand() / (float)(CFD_RAND_MAX + 1U) * (float)grid->ydim;
    }
  }
}

/* #############################################################################
 * # LBM D2Q9 plot value calculations
 * #############################################################################
 */
CFD_API CFD_INLINE float cfd_lbm_2d_calculate_density(cfd_lbm_2d_grid *grid, int x, int y)
{
  return (grid->rho[x + y * grid->xdim] - 1.0f) * 6.0f;
}

CFD_API CFD_INLINE float cfd_lbm_2d_calculate_velocity_x(cfd_lbm_2d_grid *grid, int x, int y)
{
  return grid->ux[x + y * grid->xdim] * 2.0f;
}

CFD_API CFD_INLINE float cfd_lbm_2d_calculate_velocity_y(cfd_lbm_2d_grid *grid, int x, int y)
{
  return grid->uy[x + y * grid->xdim] * 2.0f;
}

CFD_API CFD_INLINE float cfd_lbm_2d_calculate_speed(cfd_lbm_2d_grid *grid, int x, int y)
{
  return cfd_sqrtf(grid->ux[x + y * grid->xdim] * grid->ux[x + y * grid->xdim] + grid->uy[x + y * grid->xdim] * grid->uy[x + y * grid->xdim]) * 4.0f;
}

CFD_API CFD_INLINE float cfd_lbm_2d_calculate_curl(cfd_lbm_2d_grid *grid, int x, int y)
{
  return (grid->uy[x + 1 + y * grid->xdim] - grid->uy[x - 1 + y * grid->xdim] - grid->ux[x + (y + 1) * grid->xdim] + grid->ux[x + (y - 1) * grid->xdim]) * 5.0f;
}

CFD_API CFD_INLINE float cfd_lbm_2d_calculate_pressure(cfd_lbm_2d_grid *grid, int x, int y)
{
  return (grid->rho[x + y * grid->xdim] - 1.0f) * 20.0f;
}

CFD_API CFD_INLINE float cfd_lbm_2d_calculate_wall_shear_stress(cfd_lbm_2d_grid *grid, int x, int y)
{
  float shear = 0.0f;

  if (grid->barrier[x - 1 + y * grid->xdim] ||
      grid->barrier[x + 1 + y * grid->xdim] ||
      grid->barrier[x + (y - 1) * grid->xdim] ||
      grid->barrier[x + (y + 1) * grid->xdim])
  {
    shear = grid->nE[x + y * grid->xdim] + grid->nNE[x + y * grid->xdim] + grid->nSE[x + y * grid->xdim] -
            grid->nW[x + y * grid->xdim] - grid->nNW[x + y * grid->xdim] - grid->nSW[x + y * grid->xdim];
  }

  return shear * 10.0f;
}

#endif /* CFD_H */

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
