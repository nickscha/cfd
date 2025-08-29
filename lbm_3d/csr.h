/* csr.h - v0.2 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) software renderer (CSR).

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#ifndef CSR_H
#define CSR_H

/* #############################################################################
 * # COMPILER SETTINGS
 * #############################################################################
 */
/* Check if using C99 or later (inline is supported) */
#if __STDC_VERSION__ >= 199901L
#define CSR_INLINE inline
#define CSR_API static
#elif defined(__GNUC__) || defined(__clang__)
#define CSR_INLINE __inline__
#define CSR_API static
#elif defined(_MSC_VER)
#define CSR_INLINE __inline
#define CSR_API static
#else
#define CSR_INLINE
#define CSR_API static
#endif

/* If we are on a platform that does not use SSE we undefine CSR_USE_SSE if accidently enabled by the user */
#if defined(CSR_USE_SSE) && !(defined(__x86_64__) || defined(__i386__))
#undef CSR_USE_SSE
#endif

#ifdef CSR_USE_SSE
#include <xmmintrin.h>
#endif

/* #############################################################################
 * # MATRIX LAYOUT
 * #############################################################################
 */
/* Define CSR_MATRIX_ROW_MAJOR before including this file to use row-major order.
 * By default, the library uses column-major order (OpenGL standard).
 */
#ifdef CSR_MATRIX_ROW_MAJOR
#define CSR_M4X4_AT(row, col) ((row) * 4 + (col))
#else /* Column Major layout */
#define CSR_M4X4_AT(row, col) ((col) * 4 + (row))
#endif

/* #############################################################################
 * # MATH Functions
 * #############################################################################
 */
CSR_API CSR_INLINE float csr_minf(float a, float b)
{
  return (a < b) ? a : b;
}

CSR_API CSR_INLINE float csr_maxf(float a, float b)
{
  return (a > b) ? a : b;
}

CSR_API CSR_INLINE int csr_absi(int x)
{
  return (x < 0 ? -x : x);
}

CSR_API CSR_INLINE int csr_mini(int a, int b)
{
  return (a < b) ? a : b;
}

CSR_API CSR_INLINE int csr_maxi(int a, int b)
{
  return (a > b) ? a : b;
}

CSR_API CSR_INLINE void csr_pos_init(float *pos, float x, float y, float z, float w)
{
  pos[0] = x;
  pos[1] = y;
  pos[2] = z;
  pos[3] = w;
}

CSR_API CSR_INLINE void csr_v4_divf(float result[4], float v[4], float f)
{
#ifdef CSR_USE_SSE
  __m128 vec = _mm_loadu_ps(v);
  __m128 scalar = _mm_set1_ps(f);
  __m128 inv_scalar = _mm_rcp_ps(scalar);
  __m128 res_vec = _mm_mul_ps(vec, inv_scalar);

  _mm_storeu_ps(result, res_vec);
#else
  float inv_f = 1.0f / f;

  result[0] = v[0] * inv_f;
  result[1] = v[1] * inv_f;
  result[2] = v[2] * inv_f;
  result[3] = v[3] * inv_f;
#endif
}

/* Multiplies a 4x4 matrix by a 4D vector, respecting the defined matrix layout. */
CSR_API CSR_INLINE void csr_m4x4_mul_v4(float result[4], float m[16], float v[4])
{
#ifdef CSR_USE_SSE
  __m128 vec_x, vec_y, vec_z, vec_w;
  __m128 res;
  __m128 col0, col1, col2, col3;

  vec_x = _mm_set1_ps(v[0]);
  vec_y = _mm_set1_ps(v[1]);
  vec_z = _mm_set1_ps(v[2]);
  vec_w = _mm_set1_ps(v[3]);

  col0 = _mm_loadu_ps(&m[CSR_M4X4_AT(0, 0)]);
  col1 = _mm_loadu_ps(&m[CSR_M4X4_AT(0, 1)]);
  col2 = _mm_loadu_ps(&m[CSR_M4X4_AT(0, 2)]);
  col3 = _mm_loadu_ps(&m[CSR_M4X4_AT(0, 3)]);

  res = _mm_mul_ps(col0, vec_x);
  res = _mm_add_ps(res, _mm_mul_ps(col1, vec_y));
  res = _mm_add_ps(res, _mm_mul_ps(col2, vec_z));
  res = _mm_add_ps(res, _mm_mul_ps(col3, vec_w));

  _mm_storeu_ps(result, res);
#else
  result[0] = m[CSR_M4X4_AT(0, 0)] * v[0] + m[CSR_M4X4_AT(0, 1)] * v[1] + m[CSR_M4X4_AT(0, 2)] * v[2] + m[CSR_M4X4_AT(0, 3)] * v[3];
  result[1] = m[CSR_M4X4_AT(1, 0)] * v[0] + m[CSR_M4X4_AT(1, 1)] * v[1] + m[CSR_M4X4_AT(1, 2)] * v[2] + m[CSR_M4X4_AT(1, 3)] * v[3];
  result[2] = m[CSR_M4X4_AT(2, 0)] * v[0] + m[CSR_M4X4_AT(2, 1)] * v[1] + m[CSR_M4X4_AT(2, 2)] * v[2] + m[CSR_M4X4_AT(2, 3)] * v[3];
  result[3] = m[CSR_M4X4_AT(3, 0)] * v[0] + m[CSR_M4X4_AT(3, 1)] * v[1] + m[CSR_M4X4_AT(3, 2)] * v[2] + m[CSR_M4X4_AT(3, 3)] * v[3];
#endif
}

/* #############################################################################
 * # RENDERING Functions
 * #############################################################################
 */
typedef struct csr_color
{
  unsigned char r;
  unsigned char g;
  unsigned char b;

} csr_color;

typedef enum csr_render_mode
{
  CSR_RENDER_SOLID = 0,
  CSR_RENDER_WIREFRAME = 1

} csr_render_mode;

typedef enum csr_culling_mode
{

  CSR_CULLING_DISABLED = 0,      /* Default is no culling                                   */
  CSR_CULLING_CCW_BACKFACE = 1,  /* Counter-clockwise winding order, back faces are culled  */
  CSR_CULLING_CCW_FRONTFACE = 2, /* Counter-clockwise winding order, front faces are culled */
  CSR_CULLING_CW_BACKFACE = 3,   /* Clockwise winding order, back faces are culled          */
  CSR_CULLING_CW_FRONTFACE = 4   /* Clockwise winding order, front faces are culled         */

} csr_culling_mode;

typedef struct csr_context
{

  int width;              /* render area width in pixels            */
  int height;             /* render area height in pixels           */
  csr_color *framebuffer; /* memory pointer for framebuffer         */
  float *zbuffer;         /* memory pointer for zbuffer             */

} csr_context;

CSR_API CSR_INLINE unsigned long csr_memory_size(int width, int height)
{
  unsigned long area = (unsigned long)(width * height);

  return (unsigned long)(area * (unsigned long)sizeof(csr_color) + /* framebuffer size */
                         area * (unsigned long)sizeof(float)       /* zbuffer size     */
  );
}

CSR_API CSR_INLINE int csr_init_model(csr_context *context, void *memory, unsigned long memory_size, int width, int height)
{
  unsigned long memory_framebuffer_size = (unsigned long)(width * height) * (unsigned long)sizeof(csr_color);

  if (memory_size < csr_memory_size(width, height))
  {
    return 0;
  }

  context->width = width;
  context->height = height;
  context->framebuffer = (csr_color *)memory;
  context->zbuffer = (float *)((char *)memory + memory_framebuffer_size);

  return 1;
}

CSR_API CSR_INLINE csr_color csr_init_color(unsigned char r, unsigned char g, unsigned char b)
{
  csr_color result;
  result.r = r;
  result.g = g;
  result.b = b;

  return result;
}

/* Converts a point from normalized device coordinates(NDC) to screen space. */
CSR_API CSR_INLINE void csr_ndc_to_screen(csr_context *context, float result[3], float ndc_pos[4])
{
  result[0] = (ndc_pos[0] + 1.0f) * 0.5f * (float)context->width;
  result[1] = (1.0f - ndc_pos[1]) * 0.5f * (float)context->height;
  result[2] = ndc_pos[2];
}

CSR_API CSR_INLINE void csr_render_clear_screen(csr_context *context, csr_color clear_color)
{
  int size = context->width * context->height;

  int i = 0;

  for (; i + 4 <= size; i += 4)
  {
    context->framebuffer[i] = clear_color;
    context->framebuffer[i + 1] = clear_color;
    context->framebuffer[i + 2] = clear_color;
    context->framebuffer[i + 3] = clear_color;
    context->zbuffer[i] = 1.0f;
    context->zbuffer[i + 1] = 1.0f;
    context->zbuffer[i + 2] = 1.0f;
    context->zbuffer[i + 3] = 1.0f;
  }

  for (; i < size; ++i)
  {
    context->framebuffer[i] = clear_color;
    context->zbuffer[i] = 1.0f;
  }
}

/* Draws a line with depth testing using Bresenham's algorithm. */
CSR_API CSR_INLINE void csr_draw_line(csr_context *context, float p0[3], float p1[3], csr_color color)
{
  int x0 = (int)p0[0], y0 = (int)p0[1];
  int x1 = (int)p1[0], y1 = (int)p1[1];
  float z0 = p0[2], z1 = p1[2];

  int dx = csr_absi(x1 - x0);
  int sx = x0 < x1 ? 1 : -1;
  int dy = -csr_absi(y1 - y0);
  int sy = y0 < y1 ? 1 : -1;

  int err = dx + dy, e2;
  float z = z0;
  float dz = (float)(dx > -dy ? dx : -dy);

  dz = (dz == 0) ? 0.0f : (z1 - z0) / dz;

  while (1)
  {
    if (x0 >= 0 && x0 < context->width && y0 >= 0 && y0 < context->height)
    {
      int index = y0 * context->width + x0;

      if (z < context->zbuffer[index])
      {
        context->framebuffer[index] = color;
        context->zbuffer[index] = z;
      }
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
      z += dz;
    }

    if (e2 <= dx)
    {
      err += dx;
      y0 += sy;
      z += dz;
    }
  }
}

/* Fills a triangle using the barycentric coordinate method with color interpolation. */
CSR_API CSR_INLINE void csr_draw_triangle(csr_context *context, float p0[3], float p1[3], float p2[3], csr_color c0, csr_color c1, csr_color c2)
{
  /* Bounding box for the triangle */
  int min_x = (int)csr_minf(p0[0], csr_minf(p1[0], p2[0]));
  int min_y = (int)csr_minf(p0[1], csr_minf(p1[1], p2[1]));
  int max_x = (int)csr_maxf(p0[0], csr_maxf(p1[0], p2[0]));
  int max_y = (int)csr_maxf(p0[1], csr_maxf(p1[1], p2[1]));

  /* Pre-calculate constants for barycentric coordinates */
  float area = (p1[1] - p2[1]) * (p0[0] - p2[0]) + (p2[0] - p1[0]) * (p0[1] - p2[1]);

  if (area == 0.0f)
  {
    return;
  }

  /* Clamp bounding box to screen dimensions */
  min_x = csr_maxi(0, min_x);
  min_y = csr_maxi(0, min_y);
  max_x = csr_mini(context->width - 1, max_x);
  max_y = csr_mini(context->height - 1, max_y);

  {
    float inv_area = 1.0f / area;

    /* Calculate barycentric coordinate derivatives with respect to x and y */
    float w0_dx = (p1[1] - p2[1]) * inv_area;
    float w1_dx = (p2[1] - p0[1]) * inv_area;
    float w2_dx = -w0_dx - w1_dx;

    float w0_dy = (p2[0] - p1[0]) * inv_area;
    float w1_dy = (p0[0] - p2[0]) * inv_area;
    float w2_dy = -w0_dy - w1_dy;

    /* Initialize barycentric coordinates at the top-left of the bounding box */
    float w0_start = ((p1[1] - p2[1]) * ((float)min_x - p2[0]) + (p2[0] - p1[0]) * ((float)min_y - p2[1])) * inv_area;
    float w1_start = ((p2[1] - p0[1]) * ((float)min_x - p0[0]) + (p0[0] - p2[0]) * ((float)min_y - p0[1])) * inv_area;
    float w2_start = 1.0f - w0_start - w1_start;

    /* Pre-calculate color channel differences for interpolation */
    float dr_dx = (c1.r - c0.r) * w1_dx + (c2.r - c0.r) * w2_dx;
    float dg_dx = (c1.g - c0.g) * w1_dx + (c2.g - c0.g) * w2_dx;
    float db_dx = (c1.b - c0.b) * w1_dx + (c2.b - c0.b) * w2_dx;

    float dr_dy = (c1.r - c0.r) * w1_dy + (c2.r - c0.r) * w2_dy;
    float dg_dy = (c1.g - c0.g) * w1_dy + (c2.g - c0.g) * w2_dy;
    float db_dy = (c1.b - c0.b) * w1_dy + (c2.b - c0.b) * w2_dy;

    float r_start = c0.r + (c1.r - c0.r) * w1_start + (c2.r - c0.r) * w2_start;
    float g_start = c0.g + (c1.g - c0.g) * w1_start + (c2.g - c0.g) * w2_start;
    float b_start = c0.b + (c1.b - c0.b) * w1_start + (c2.b - c0.b) * w2_start;

    int x, y;

    for (y = min_y; y <= max_y; ++y)
    {
      float w0 = w0_start;
      float w1 = w1_start;
      float w2 = w2_start;

      float current_r = r_start;
      float current_g = g_start;
      float current_b = b_start;

      int index_row_start = y * context->width + min_x;

      for (x = min_x; x <= max_x; ++x)
      {
        if (w0 >= 0.0f && w1 >= 0.0f && w2 >= 0.0f)
        {
          /* Interpolate Z-depth using w values */
          float z = p0[2] * w0 + p1[2] * w1 + p2[2] * w2;

          int index = index_row_start + (x - min_x);

          /* Depth testing: only draw if the new pixel is closer than the existing one */
          if (z < context->zbuffer[index])
          {
            csr_color pixel_color;
            pixel_color.r = (unsigned char)current_r;
            pixel_color.g = (unsigned char)current_g;
            pixel_color.b = (unsigned char)current_b;

            context->framebuffer[index] = pixel_color;
            context->zbuffer[index] = z;
          }
        }

        /* Increment barycentric coordinates and colors with pre-calculated deltas */
        w0 += w0_dx;
        w1 += w1_dx;
        w2 += w2_dx;
        current_r += dr_dx;
        current_g += dg_dx;
        current_b += db_dx;
      }

      /* Reset w values and colors for the start of the next row */
      w0_start += w0_dy;
      w1_start += w1_dy;
      w2_start += w2_dy;
      r_start += dr_dy;
      g_start += dg_dy;
      b_start += db_dy;
    }
  }
}

CSR_API CSR_INLINE void csr_render(csr_context *context, csr_render_mode render_mode, csr_culling_mode culling_mode, int stride, float *vertices, unsigned long num_vertices, int *indices, unsigned long num_indices, float projection_view_model_matrix[16], csr_color c0, csr_color c1, csr_color c2)
{
  unsigned long i;

  (void)num_vertices;

  for (i = 0; i < num_indices; i += 3)
  {
    /* Get vertex indices for the current triangle */
    int i0 = indices[i];
    int i1 = indices[i + 1];
    int i2 = indices[i + 2];

    /* Get the vertex data from the main vertex array using the indices, and convert to homogeneous coordinates */
    float pos0[4];
    float pos1[4];
    float pos2[4];

    /* 1. Vertex Processing (Model, View, Projection) */
    float v0_transformed[4];
    float v1_transformed[4];
    float v2_transformed[4];

    float v0_ndc[4];
    float v1_ndc[4];
    float v2_ndc[4];

    float v0_screen[3];
    float v1_screen[3];
    float v2_screen[3];

    csr_pos_init(pos0, vertices[i0 * stride + 0], vertices[i0 * stride + 1], vertices[i0 * stride + 2], 1.0f);
    csr_pos_init(pos1, vertices[i1 * stride + 0], vertices[i1 * stride + 1], vertices[i1 * stride + 2], 1.0f);
    csr_pos_init(pos2, vertices[i2 * stride + 0], vertices[i2 * stride + 1], vertices[i2 * stride + 2], 1.0f);

    csr_m4x4_mul_v4(v0_transformed, projection_view_model_matrix, pos0);
    csr_m4x4_mul_v4(v1_transformed, projection_view_model_matrix, pos1);
    csr_m4x4_mul_v4(v2_transformed, projection_view_model_matrix, pos2);

    /* Check if the triangle is behind the camera (clipping) */
    if (v0_transformed[3] <= 0.0f || v1_transformed[3] <= 0.0f || v2_transformed[3] <= 0.0f)
    {
      continue;
    }

    /* 2. Perspective Divide (Clip Space to NDC) */
    csr_v4_divf(v0_ndc, v0_transformed, v0_transformed[3]);
    csr_v4_divf(v1_ndc, v1_transformed, v1_transformed[3]);
    csr_v4_divf(v2_ndc, v2_transformed, v2_transformed[3]);

    /* 3. Viewport Transform (NDC to Screen Space) */
    csr_ndc_to_screen(context, v0_screen, v0_ndc);
    csr_ndc_to_screen(context, v1_screen, v1_ndc);
    csr_ndc_to_screen(context, v2_screen, v2_ndc);

    /* 4. Culling based on winding order */
    if (culling_mode != CSR_CULLING_DISABLED)
    {
      float ax = v1_screen[0] - v0_screen[0];
      float ay = v1_screen[1] - v0_screen[1];
      float bx = v2_screen[0] - v0_screen[0];
      float by = v2_screen[1] - v0_screen[1];
      float face = ax * by - ay * bx;

      int is_ccw_face = (face >= 0.0f);
      int is_cw_face = (face <= 0.0f);

      int should_cull = 0;

      should_cull |= (culling_mode == CSR_CULLING_CCW_BACKFACE) & is_cw_face;
      should_cull |= (culling_mode == CSR_CULLING_CCW_FRONTFACE) & is_ccw_face;
      should_cull |= (culling_mode == CSR_CULLING_CW_BACKFACE) & is_ccw_face;
      should_cull |= (culling_mode == CSR_CULLING_CW_FRONTFACE) & is_cw_face;

      if (should_cull)
      {
        continue;
      }
    }

    /* 5. Rasterization & Depth Testing */
    if (render_mode == CSR_RENDER_SOLID)
    {
      csr_draw_triangle(context, v0_screen, v1_screen, v2_screen, c0, c1, c2);
    }
    else
    {
      csr_draw_line(context, v0_screen, v1_screen, c0);
      csr_draw_line(context, v1_screen, v2_screen, c0);
      csr_draw_line(context, v2_screen, v0_screen, c0);
    }
  }
}

#endif /* CSR_H */

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
