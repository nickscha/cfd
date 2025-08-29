/* cfd.h - v0.2 - public domain data structures - nickscha 2025

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

CFD_API CFD_INLINE float cfd_cbrtf(float x)
{
  float guess;
  int i;

  if (x == 0.0f)
  {
    return 0.0f;
  }

  guess = x;

  /* Use Newton's method to find the cube root */
  for (i = 0; i < 20; ++i)
  {
    guess = (2.0f * guess + x / (guess * guess)) / 3.0f;
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

  float omega;

  /* Barrier Force data */
  int barrierCount;
  float barrierFx, barrierFy;
  float barrierxSum, barrierySum;

  /* Tracer data */
  float tracerX[CFD_LBM_2D_NUMBER_TRACERS];
  float tracerY[CFD_LBM_2D_NUMBER_TRACERS];

  /* Particle distributions */
  float *nC, *nN, *nS, *nE, *nW, *nNE, *nSE, *nNW, *nSW;

  float *rho; /* Density    */
  float *ux;  /* Velocity X */
  float *uy;  /* Velocity Y */

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

CFD_API CFD_INLINE void cfd_lbm_2d_init_grid(cfd_lbm_2d_grid *grid, void *memory, int xdim, int ydim, float omega)
{
  char *ptr = (char *)memory;

  unsigned long nSites = (unsigned long)(xdim * ydim);
  unsigned long sizeFloat = sizeof(float);
  unsigned long dist_size = nSites * sizeFloat;

  grid->xdim = xdim;
  grid->ydim = ydim;
  grid->omega = omega;

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
CFD_API CFD_INLINE void cfd_lbm_2d_collide(cfd_lbm_2d_grid *grid)
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
      float omega = grid->omega;

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

CFD_API CFD_INLINE void cfd_lbm_2d_collide_and_stream(cfd_lbm_2d_grid *grid)
{
  cfd_lbm_2d_collide(grid);
  cfd_lbm_2d_stream(grid);
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

/* #############################################################################
 * # LBM D3Q19 Model
 * #############################################################################
 */
#define CFD_LBM_3D_W0 (1.0f / 3.0f)   /* Weight for center particle */
#define CFD_LBM_3D_W1 (1.0f / 18.0f)  /* Weight for axial particles */
#define CFD_LBM_3D_W2 (1.0f / 36.0f)  /* Weight for diagonal particles */
#define CFD_LBM_3D_NUMBER_TRACERS 144 /* Must be <= 256 for this implementation */

typedef struct cfd_lbm_3d_grid
{
  /* Grid size */
  int xdim, ydim, zdim;
  float omega;

  /* Barrier Force data */
  int barrierCount;
  float barrierFx, barrierFy, barrierFz;
  float barrierxSum, barrierySum, barrierzSum;

  /* Tracer data */
  float tracerX[CFD_LBM_3D_NUMBER_TRACERS];
  float tracerY[CFD_LBM_3D_NUMBER_TRACERS];
  float tracerZ[CFD_LBM_3D_NUMBER_TRACERS];

  /* Particle distributions (19 directions) */
  float *f0, *fE, *fW, *fN, *fS, *fT, *fB;
  float *fNE, *fNW, *fSE, *fSW;
  float *fET, *fWT, *fEB, *fWB;
  float *fNT, *fST, *fNB, *fSB;

  /* Macroscopic properties */
  float *rho; /* Density    */
  float *ux;  /* Velocity X */
  float *uy;  /* Velocity Y */
  float *uz;  /* Velocity Z */

  /* Barrier grid */
  unsigned char *barrier;

} cfd_lbm_3d_grid;

CFD_API CFD_INLINE unsigned long cfd_lbm_3d_grid_memory_size(int xdim, int ydim, int zdim)
{
  unsigned long nSites = (unsigned long)(xdim * ydim * zdim);
  return (nSites * sizeof(float) * 23 + nSites * sizeof(unsigned char));
}

CFD_API CFD_INLINE void cfd_lbm_3d_init_grid(cfd_lbm_3d_grid *grid, void *memory, int xdim, int ydim, int zdim, float omega)
{
  char *ptr = (char *)memory;
  unsigned long nSites = (unsigned long)(xdim * ydim * zdim);
  unsigned long dist_size = nSites * sizeof(float);

  grid->xdim = xdim;
  grid->ydim = ydim;
  grid->zdim = zdim;
  grid->omega = omega;

  /* Assign pointers for 19 distribution arrays */
  grid->f0 = (float *)ptr;
  ptr += dist_size;
  grid->fE = (float *)ptr;
  ptr += dist_size;
  grid->fW = (float *)ptr;
  ptr += dist_size;
  grid->fN = (float *)ptr;
  ptr += dist_size;
  grid->fS = (float *)ptr;
  ptr += dist_size;
  grid->fT = (float *)ptr;
  ptr += dist_size;
  grid->fB = (float *)ptr;
  ptr += dist_size;
  grid->fNE = (float *)ptr;
  ptr += dist_size;
  grid->fNW = (float *)ptr;
  ptr += dist_size;
  grid->fSE = (float *)ptr;
  ptr += dist_size;
  grid->fSW = (float *)ptr;
  ptr += dist_size;
  grid->fET = (float *)ptr;
  ptr += dist_size;
  grid->fWT = (float *)ptr;
  ptr += dist_size;
  grid->fEB = (float *)ptr;
  ptr += dist_size;
  grid->fWB = (float *)ptr;
  ptr += dist_size;
  grid->fNT = (float *)ptr;
  ptr += dist_size;
  grid->fST = (float *)ptr;
  ptr += dist_size;
  grid->fNB = (float *)ptr;
  ptr += dist_size;
  grid->fSB = (float *)ptr;
  ptr += dist_size;

  /* Assign pointers for macroscopic properties */
  grid->rho = (float *)ptr;
  ptr += dist_size;
  grid->ux = (float *)ptr;
  ptr += dist_size;
  grid->uy = (float *)ptr;
  ptr += dist_size;
  grid->uz = (float *)ptr;
  ptr += dist_size;

  /* Assign pointer for barrier grid */
  grid->barrier = (unsigned char *)ptr;
}

CFD_API CFD_INLINE void cfd_lbm_3d_init_barriers(cfd_lbm_3d_grid *grid)
{
  int barrierSize = 8;
  int x = grid->xdim / 3;
  int y, z;
  int xdim = grid->xdim;
  int ydim = grid->ydim;

  /* Create a flat plate barrier */
  for (z = (grid->zdim / 2) - barrierSize; z <= (grid->zdim / 2) + barrierSize; ++z)
  {
    for (y = (grid->ydim / 2) - barrierSize; y <= (grid->ydim / 2) + barrierSize; ++y)
    {
      if (x >= 0 && x < grid->xdim && y >= 0 && y < grid->ydim && z >= 0 && z < grid->zdim)
      {
        grid->barrier[x + y * xdim + z * xdim * ydim] = 1;
      }
    }
  }
}

CFD_API CFD_INLINE void cfd_lbm_3d_init_equilibrium(cfd_lbm_3d_grid *grid, int x, int y, int z, float newux, float newuy, float newuz, float newrho)
{
  int i = x + y * grid->xdim + z * grid->xdim * grid->ydim;
  float u_sq = newux * newux + newuy * newuy + newuz * newuz;
  float term_u_sq = 1.5f * u_sq;

  /* Center */
  grid->f0[i] = CFD_LBM_3D_W0 * newrho * (1.0f - term_u_sq);

  /* Axial */
  grid->fE[i] = CFD_LBM_3D_W1 * newrho * (1.0f + 3.0f * newux + 4.5f * newux * newux - term_u_sq);
  grid->fW[i] = CFD_LBM_3D_W1 * newrho * (1.0f - 3.0f * newux + 4.5f * newux * newux - term_u_sq);
  grid->fN[i] = CFD_LBM_3D_W1 * newrho * (1.0f + 3.0f * newuy + 4.5f * newuy * newuy - term_u_sq);
  grid->fS[i] = CFD_LBM_3D_W1 * newrho * (1.0f - 3.0f * newuy + 4.5f * newuy * newuy - term_u_sq);
  grid->fT[i] = CFD_LBM_3D_W1 * newrho * (1.0f + 3.0f * newuz + 4.5f * newuz * newuz - term_u_sq);
  grid->fB[i] = CFD_LBM_3D_W1 * newrho * (1.0f - 3.0f * newuz + 4.5f * newuz * newuz - term_u_sq);

  /* Diagonal */
  grid->fNE[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (newux + newuy) + 4.5f * (newux + newuy) * (newux + newuy) - term_u_sq);
  grid->fNW[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (-newux + newuy) + 4.5f * (-newux + newuy) * (-newux + newuy) - term_u_sq);
  grid->fSE[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (newux - newuy) + 4.5f * (newux - newuy) * (newux - newuy) - term_u_sq);
  grid->fSW[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (-newux - newuy) + 4.5f * (-newux - newuy) * (-newux - newuy) - term_u_sq);
  grid->fET[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (newux + newuz) + 4.5f * (newux + newuz) * (newux + newuz) - term_u_sq);
  grid->fWT[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (-newux + newuz) + 4.5f * (-newux + newuz) * (-newux + newuz) - term_u_sq);
  grid->fEB[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (newux - newuz) + 4.5f * (newux - newuz) * (newux - newuz) - term_u_sq);
  grid->fWB[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (-newux - newuz) + 4.5f * (-newux - newuz) * (-newux - newuz) - term_u_sq);
  grid->fNT[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (newuy + newuz) + 4.5f * (newuy + newuz) * (newuy + newuz) - term_u_sq);
  grid->fST[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (-newuy + newuz) + 4.5f * (-newuy + newuz) * (-newuy + newuz) - term_u_sq);
  grid->fNB[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (newuy - newuz) + 4.5f * (newuy - newuz) * (newuy - newuz) - term_u_sq);
  grid->fSB[i] = CFD_LBM_3D_W2 * newrho * (1.0f + 3.0f * (-newuy - newuz) + 4.5f * (-newuy - newuz) * (-newuy - newuz) - term_u_sq);

  grid->rho[i] = newrho;
  grid->ux[i] = newux;
  grid->uy[i] = newuy;
  grid->uz[i] = newuz;
}

CFD_API CFD_INLINE void cfd_lbm_3d_init_fluid(cfd_lbm_3d_grid *grid, float u0)
{
  int z;

  for (z = 0; z < grid->zdim; ++z)
  {
    int y;

    for (y = 0; y < grid->ydim; ++y)
    {
      int x;

      for (x = 0; x < grid->xdim; ++x)
      {
        cfd_lbm_3d_init_equilibrium(grid, x, y, z, u0, 0.0f, 0.0f, 1.0f);
      }
    }
  }
}

CFD_API CFD_INLINE void cfd_lbm_3d_init_tracers(cfd_lbm_3d_grid *grid)
{
  int nRows = (int)(cfd_cbrtf((float)CFD_LBM_3D_NUMBER_TRACERS));
  float dx, dy, dz;
  float nextX, nextY, nextZ;
  int t;
  if (nRows <= 0)
    nRows = 1;

  dx = (float)grid->xdim / (float)nRows;
  dy = (float)grid->ydim / (float)nRows;
  dz = (float)grid->zdim / (float)nRows;

  nextX = dx * 0.5f;
  nextY = dy * 0.5f;
  nextZ = dz * 0.5f;

  for (t = 0; t < CFD_LBM_3D_NUMBER_TRACERS; ++t)
  {
    grid->tracerX[t] = nextX;
    grid->tracerY[t] = nextY;
    grid->tracerZ[t] = nextZ;
    nextX += dx;
    if (nextX >= (float)grid->xdim)
    {
      nextX = dx * 0.5f;
      nextY += dy;
      if (nextY >= (float)grid->ydim)
      {
        nextY = dy * 0.5f;
        nextZ += dz;
      }
    }
  }
}

CFD_API CFD_INLINE void cfd_lbm_3d_collide(cfd_lbm_3d_grid *grid)
{
  int x, y, z;
  int xdim = grid->xdim;
  int ydim = grid->ydim;
  int zdim = grid->zdim;
  float omega = grid->omega;

  for (z = 1; z < zdim - 1; ++z)
  {
    for (y = 1; y < ydim - 1; ++y)
    {
      for (x = 1; x < xdim - 1; ++x)
      {
        int i = x + y * xdim + z * xdim * ydim;
        float thisrho, thisux, thisuy, thisuz;
        float fEq[19];
        float u_sq, term_u_sq;

        /* Calculate macroscopic properties */
        thisrho = grid->f0[i] + grid->fE[i] + grid->fW[i] + grid->fN[i] + grid->fS[i] + grid->fT[i] + grid->fB[i] +
                  grid->fNE[i] + grid->fNW[i] + grid->fSE[i] + grid->fSW[i] +
                  grid->fET[i] + grid->fWT[i] + grid->fEB[i] + grid->fWB[i] +
                  grid->fNT[i] + grid->fST[i] + grid->fNB[i] + grid->fSB[i];

        thisux = (grid->fE[i] - grid->fW[i] + grid->fNE[i] - grid->fNW[i] + grid->fSE[i] - grid->fSW[i] +
                  grid->fET[i] - grid->fWT[i] + grid->fEB[i] - grid->fWB[i]) /
                 thisrho;

        thisuy = (grid->fN[i] - grid->fS[i] + grid->fNE[i] + grid->fNW[i] - grid->fSE[i] - grid->fSW[i] +
                  grid->fNT[i] - grid->fST[i] + grid->fNB[i] - grid->fSB[i]) /
                 thisrho;

        thisuz = (grid->fT[i] - grid->fB[i] + grid->fET[i] + grid->fWT[i] - grid->fEB[i] - grid->fWB[i] +
                  grid->fNT[i] + grid->fST[i] - grid->fNB[i] - grid->fSB[i]) /
                 thisrho;

        grid->rho[i] = thisrho;
        grid->ux[i] = thisux;
        grid->uy[i] = thisuy;
        grid->uz[i] = thisuz;

        /* Calculate equilibrium distribution */
        u_sq = thisux * thisux + thisuy * thisuy + thisuz * thisuz;
        term_u_sq = 1.5f * u_sq;

        fEq[0] = CFD_LBM_3D_W0 * thisrho * (1.0f - term_u_sq);
        fEq[1] = CFD_LBM_3D_W1 * thisrho * (1.0f + 3.0f * thisux + 4.5f * thisux * thisux - term_u_sq);
        fEq[2] = CFD_LBM_3D_W1 * thisrho * (1.0f - 3.0f * thisux + 4.5f * thisux * thisux - term_u_sq);
        fEq[3] = CFD_LBM_3D_W1 * thisrho * (1.0f + 3.0f * thisuy + 4.5f * thisuy * thisuy - term_u_sq);
        fEq[4] = CFD_LBM_3D_W1 * thisrho * (1.0f - 3.0f * thisuy + 4.5f * thisuy * thisuy - term_u_sq);
        fEq[5] = CFD_LBM_3D_W1 * thisrho * (1.0f + 3.0f * thisuz + 4.5f * thisuz * thisuz - term_u_sq);
        fEq[6] = CFD_LBM_3D_W1 * thisrho * (1.0f - 3.0f * thisuz + 4.5f * thisuz * thisuz - term_u_sq);
        fEq[7] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (thisux + thisuy) + 4.5f * (thisux + thisuy) * (thisux + thisuy) - term_u_sq);
        fEq[8] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (-thisux + thisuy) + 4.5f * (-thisux + thisuy) * (-thisux + thisuy) - term_u_sq);
        fEq[9] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (thisux - thisuy) + 4.5f * (thisux - thisuy) * (thisux - thisuy) - term_u_sq);
        fEq[10] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (-thisux - thisuy) + 4.5f * (-thisux - thisuy) * (-thisux - thisuy) - term_u_sq);
        fEq[11] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (thisux + thisuz) + 4.5f * (thisux + thisuz) * (thisux + thisuz) - term_u_sq);
        fEq[12] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (-thisux + thisuz) + 4.5f * (-thisux + thisuz) * (-thisux + thisuz) - term_u_sq);
        fEq[13] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (thisux - thisuz) + 4.5f * (thisux - thisuz) * (thisux - thisuz) - term_u_sq);
        fEq[14] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (-thisux - thisuz) + 4.5f * (-thisux - thisuz) * (-thisux - thisuz) - term_u_sq);
        fEq[15] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (thisuy + thisuz) + 4.5f * (thisuy + thisuz) * (thisuy + thisuz) - term_u_sq);
        fEq[16] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (-thisuy + thisuz) + 4.5f * (-thisuy + thisuz) * (-thisuy + thisuz) - term_u_sq);
        fEq[17] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (thisuy - thisuz) + 4.5f * (thisuy - thisuz) * (thisuy - thisuz) - term_u_sq);
        fEq[18] = CFD_LBM_3D_W2 * thisrho * (1.0f + 3.0f * (-thisuy - thisuz) + 4.5f * (-thisuy - thisuz) * (-thisuy - thisuz) - term_u_sq);

        /* Perform collision */
        grid->f0[i] += omega * (fEq[0] - grid->f0[i]);
        grid->fE[i] += omega * (fEq[1] - grid->fE[i]);
        grid->fW[i] += omega * (fEq[2] - grid->fW[i]);
        grid->fN[i] += omega * (fEq[3] - grid->fN[i]);
        grid->fS[i] += omega * (fEq[4] - grid->fS[i]);
        grid->fT[i] += omega * (fEq[5] - grid->fT[i]);
        grid->fB[i] += omega * (fEq[6] - grid->fB[i]);
        grid->fNE[i] += omega * (fEq[7] - grid->fNE[i]);
        grid->fNW[i] += omega * (fEq[8] - grid->fNW[i]);
        grid->fSE[i] += omega * (fEq[9] - grid->fSE[i]);
        grid->fSW[i] += omega * (fEq[10] - grid->fSW[i]);
        grid->fET[i] += omega * (fEq[11] - grid->fET[i]);
        grid->fWT[i] += omega * (fEq[12] - grid->fWT[i]);
        grid->fEB[i] += omega * (fEq[13] - grid->fEB[i]);
        grid->fWB[i] += omega * (fEq[14] - grid->fWB[i]);
        grid->fNT[i] += omega * (fEq[15] - grid->fNT[i]);
        grid->fST[i] += omega * (fEq[16] - grid->fST[i]);
        grid->fNB[i] += omega * (fEq[17] - grid->fNB[i]);
        grid->fSB[i] += omega * (fEq[18] - grid->fSB[i]);
      }
    }
  }
}

CFD_API CFD_INLINE void cfd_lbm_3d_stream(cfd_lbm_3d_grid *grid)
{
  int x, y, z;
  int xdim = grid->xdim;
  int ydim = grid->ydim;
  int zdim = grid->zdim;
  int xdim_ydim = xdim * ydim;

  grid->barrierCount = 0;
  grid->barrierxSum = 0.0f;
  grid->barrierySum = 0.0f;
  grid->barrierzSum = 0.0f;
  grid->barrierFx = 0.0f;
  grid->barrierFy = 0.0f;
  grid->barrierFz = 0.0f;

  /* This streaming is done in-place with careful loop ordering */
  for (z = zdim - 2; z > 0; --z)
  {
    for (y = ydim - 2; y > 0; --y)
    {
      for (x = xdim - 2; x > 0; --x)
      {
        int i = x + y * xdim + z * xdim_ydim;
        grid->fE[i] = grid->fE[i - 1];
        grid->fN[i] = grid->fN[i - xdim];
        grid->fT[i] = grid->fT[i - xdim_ydim];
        grid->fNE[i] = grid->fNE[i - 1 - xdim];
        grid->fET[i] = grid->fET[i - 1 - xdim_ydim];
        grid->fNT[i] = grid->fNT[i - xdim - xdim_ydim];
      }
    }
  }

  for (z = zdim - 2; z > 0; --z)
  {
    for (y = ydim - 2; y > 0; --y)
    {
      for (x = 1; x < xdim - 1; ++x)
      {
        int i = x + y * xdim + z * xdim_ydim;
        grid->fW[i] = grid->fW[i + 1];
        grid->fNW[i] = grid->fNW[i + 1 - xdim];
        grid->fWT[i] = grid->fWT[i + 1 - xdim_ydim];
      }
    }
  }

  for (z = zdim - 2; z > 0; --z)
  {
    for (y = 1; y < ydim - 1; ++y)
    {
      for (x = xdim - 2; x > 0; --x)
      {
        int i = x + y * xdim + z * xdim_ydim;
        grid->fS[i] = grid->fS[i + xdim];
        grid->fSE[i] = grid->fSE[i - 1 + xdim];
        grid->fST[i] = grid->fST[i + xdim - xdim_ydim];
      }
    }
  }

  for (z = zdim - 2; z > 0; --z)
  {
    for (y = 1; y < ydim - 1; ++y)
    {
      for (x = 1; x < xdim - 1; ++x)
      {
        int i = x + y * xdim + z * xdim_ydim;
        grid->fSW[i] = grid->fSW[i + 1 + xdim];
      }
    }
  }

  for (z = 1; z < zdim - 1; ++z)
  {
    for (y = ydim - 2; y > 0; --y)
    {
      for (x = xdim - 2; x > 0; --x)
      {
        int i = x + y * xdim + z * xdim_ydim;
        grid->fB[i] = grid->fB[i + xdim_ydim];
        grid->fEB[i] = grid->fEB[i - 1 + xdim_ydim];
        grid->fNB[i] = grid->fNB[i - xdim + xdim_ydim];
      }
    }
  }

  for (z = 1; z < zdim - 1; ++z)
  {
    for (y = ydim - 2; y > 0; --y)
    {
      for (x = 1; x < xdim - 1; ++x)
      {
        int i = x + y * xdim + z * xdim_ydim;
        grid->fWB[i] = grid->fWB[i + 1 + xdim_ydim];
      }
    }
  }
  for (z = 1; z < zdim - 1; ++z)
  {
    for (y = 1; y < ydim - 1; ++y)
    {
      for (x = xdim - 2; x > 0; --x)
      {
        int i = x + y * xdim + z * xdim_ydim;
        grid->fSB[i] = grid->fSB[i + xdim + xdim_ydim];
      }
    }
  }

  /* Handle bounce-back from barriers */
  for (z = 1; z < zdim - 1; ++z)
  {
    for (y = 1; y < ydim - 1; ++y)
    {
      for (x = 1; x < xdim - 1; ++x)
      {
        int i = x + y * xdim + z * xdim_ydim;

        if (grid->barrier[i])
        {
          float t;

          t = grid->fE[i];

          grid->fE[i] = grid->fW[i];
          grid->fW[i] = t;

          t = grid->fN[i];

          grid->fN[i] = grid->fS[i];
          grid->fS[i] = t;

          t = grid->fT[i];

          grid->fT[i] = grid->fB[i];
          grid->fB[i] = t;

          t = grid->fNE[i];

          grid->fNE[i] = grid->fSW[i];
          grid->fSW[i] = t;

          t = grid->fNW[i];

          grid->fNW[i] = grid->fSE[i];
          grid->fSE[i] = t;

          t = grid->fET[i];

          grid->fET[i] = grid->fWB[i];
          grid->fWB[i] = t;

          t = grid->fWT[i];

          grid->fWT[i] = grid->fEB[i];
          grid->fEB[i] = t;

          t = grid->fNT[i];

          grid->fNT[i] = grid->fSB[i];
          grid->fSB[i] = t;

          t = grid->fST[i];

          grid->fST[i] = grid->fNB[i];
          grid->fNB[i] = t;

          /* Sum forces on barrier sites */
          grid->barrierCount++;
          grid->barrierxSum += (float)x;
          grid->barrierySum += (float)y;
          grid->barrierzSum += (float)z;
          grid->barrierFx += 2.0f * (grid->fE[i] + grid->fNE[i] + grid->fSE[i] + grid->fET[i] + grid->fEB[i] -
                                     (grid->fW[i] + grid->fNW[i] + grid->fSW[i] + grid->fWT[i] + grid->fWB[i]));
          grid->barrierFy += 2.0f * (grid->fN[i] + grid->fNE[i] + grid->fNW[i] + grid->fNT[i] + grid->fNB[i] -
                                     (grid->fS[i] + grid->fSE[i] + grid->fSW[i] + grid->fST[i] + grid->fSB[i]));
          grid->barrierFz += 2.0f * (grid->fT[i] + grid->fET[i] + grid->fWT[i] + grid->fNT[i] + grid->fST[i] -
                                     (grid->fB[i] + grid->fEB[i] + grid->fWB[i] + grid->fNB[i] + grid->fSB[i]));
        }
      }
    }
  }
}

CFD_API CFD_INLINE void cfd_lbm_3d_collide_and_stream(cfd_lbm_3d_grid *grid)
{
  cfd_lbm_3d_collide(grid);
  cfd_lbm_3d_stream(grid);
}

CFD_API CFD_INLINE void cfd_lbm_3d_move_tracers(cfd_lbm_3d_grid *grid)
{
  int t;
  for (t = 0; t < CFD_LBM_3D_NUMBER_TRACERS; ++t)
  {
    int roundedX = (int)(grid->tracerX[t] + 0.5f);
    int roundedY = (int)(grid->tracerY[t] + 0.5f);
    int roundedZ = (int)(grid->tracerZ[t] + 0.5f);
    int index;

    if (roundedX < 0 || roundedX >= grid->xdim || roundedY < 0 || roundedY >= grid->ydim || roundedZ < 0 || roundedZ >= grid->zdim)
    {
      continue;
    }

    index = roundedX + roundedY * grid->xdim + roundedZ * grid->xdim * grid->ydim;
    grid->tracerX[t] += grid->ux[index];
    grid->tracerY[t] += grid->uy[index];
    grid->tracerZ[t] += grid->uz[index];

    /* Wrap tracers around boundaries */
    if (grid->tracerX[t] > (float)(grid->xdim - 1))
    {
      grid->tracerX[t] = 0.0f;
      grid->tracerY[t] = (float)cfd_rand() / (float)(CFD_RAND_MAX + 1U) * (float)grid->ydim;
      grid->tracerZ[t] = (float)cfd_rand() / (float)(CFD_RAND_MAX + 1U) * (float)grid->zdim;
    }
  }
}

/* #############################################################################
 * # LBM D3Q19 plot value calculations
 * #############################################################################
 */
CFD_API CFD_INLINE float cfd_lbm_3d_calculate_speed(cfd_lbm_3d_grid *grid, int x, int y, int z)
{
  int i = x + y * grid->xdim + z * grid->xdim * grid->ydim;
  float u_sq = grid->ux[i] * grid->ux[i] + grid->uy[i] * grid->uy[i] + grid->uz[i] * grid->uz[i];
  return cfd_sqrtf(u_sq) * 4.0f;
}

CFD_API CFD_INLINE float cfd_lbm_3d_calculate_curl(cfd_lbm_3d_grid *grid, int x, int y, int z)
{
  int xdim = grid->xdim;
  int ydim = grid->ydim;
  int xdim_ydim = xdim * ydim;
  float curl_x, curl_y, curl_z;

  curl_x = grid->uz[x + (y + 1) * xdim + z * xdim_ydim] - grid->uz[x + (y - 1) * xdim + z * xdim_ydim] -
           (grid->uy[x + y * xdim + (z + 1) * xdim_ydim] - grid->uy[x + y * xdim + (z - 1) * xdim_ydim]);
  curl_y = grid->ux[x + y * xdim + (z + 1) * xdim_ydim] - grid->ux[x + y * xdim + (z - 1) * xdim_ydim] -
           (grid->uz[(x + 1) + y * xdim + z * xdim_ydim] - grid->uz[(x - 1) + y * xdim + z * xdim_ydim]);
  curl_z = grid->uy[(x + 1) + y * xdim + z * xdim_ydim] - grid->uy[(x - 1) + y * xdim + z * xdim_ydim] -
           (grid->ux[x + (y + 1) * xdim + z * xdim_ydim] - grid->ux[x + (y - 1) * xdim + z * xdim_ydim]);

  return cfd_sqrtf(curl_x * curl_x + curl_y * curl_y + curl_z * curl_z) * 5.0f;
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
