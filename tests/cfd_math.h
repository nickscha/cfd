/* cfd.h - v0.1 - public domain data structures - nickscha 2025

A C89 standard compliant, single header, nostdlib (no C Standard Library) math helper functions.

LICENSE

  Placed in the public domain and also MIT licensed.
  See end of file for detailed license information.

*/
#ifndef CFD_MATH_H
#define CFD_MATH_H

/* #############################################################################
 * # COMPILER SETTINGS
 * #############################################################################
 */
/* Check if using C99 or later (inline is supported) */
#if __STDC_VERSION__ >= 199901L
#define CFD_MATH_INLINE inline
#define CFD_MATH_API extern
#elif defined(__GNUC__) || defined(__clang__)
#define CFD_MATH_INLINE __inline__
#define CFD_MATH_API static
#elif defined(_MSC_VER)
#define CFD_MATH_INLINE __inline
#define CFD_MATH_API static
#else
#define CFD_MATH_INLINE
#define CFD_MATH_API static
#endif

CFD_MATH_API CFD_MATH_INLINE int cfd_abs(int x)
{
  return (x < 0) ? -x : x;
}

CFD_MATH_API CFD_MATH_INLINE float cfd_sinf(float x)
{
  float term = x;
  float result = x;
  float x2 = x * x;
  int i;

  for (i = 1; i <= 10; ++i)
  {
    term *= -x2 / ((float)(2 * i) * (float)(2 * i + 1));
    result += term;
  }

  return result;
}

CFD_MATH_API CFD_MATH_INLINE float cfd_cosf(float x)
{
  float term = 1.0f;
  float result = 1.0f;
  float x2 = x * x;
  int i;

  for (i = 1; i <= 10; ++i)
  {
    term *= -x2 / ((float)(2 * i - 1) * (float)(2 * i));
    result += term;
  }

  return result;
}

CFD_MATH_API CFD_MATH_INLINE float cfd_atanf(float z)
{
  int sign = 1;
  float result;

  if (z < 0.0f)
  {
    sign = -1;
    z = -z;
  }

  if (z > 1.0f)
  {
    result = 1.57079632679f - cfd_atanf(1.0f / z); /* pi/2 - atan(1/z) */
  }
  else
  {
    float z2 = z * z;
    result = z * (1.0f - z2 * (1.0f / 3.0f - z2 * (1.0f / 5.0f - z2 * (1.0f / 7.0f))));
  }

  return (float)sign * result;
}

CFD_MATH_API CFD_MATH_INLINE float cfd_atan2f(float y, float x)
{
  const float PI = 3.14159265359f;

  if (x > 0.0f)
    return cfd_atanf(y / x);
  if (x < 0.0f && y >= 0.0f)
    return cfd_atanf(y / x) + PI;
  if (x < 0.0f && y < 0.0f)
    return cfd_atanf(y / x) - PI;
  if (x == 0.0f && y > 0.0f)
    return PI / 2.0f;
  if (x == 0.0f && y < 0.0f)
    return -PI / 2.0f;
  return 0.0f;
}

CFD_MATH_API CFD_MATH_INLINE float cfd_floorf(float x)
{
  int i = (int)x;
  if (x < (float)i)
  {
    return (float)(i - 1);
  }
  return (float)i;
}

CFD_MATH_API CFD_MATH_INLINE float cfd_powf(float base, float exp)
{
  /* Only handles positive base for now */
  int i;
  float result = 1.0f;

  if (exp == 0.0f)
  {
    return 1.0f;
  }
  if (exp == 1.0f)
  {
    return base;
  }

  /*
   * Check if exp is an integer.
   * This is a robust check, but the original logic also worked for this specific use case.
   */
  if (cfd_floorf(exp) == exp)
  {
    /* Integer exponent */
    int e = (int)exp;

    /* Handle positive exponents */
    if (e > 0)
    {
      for (i = 0; i < e; ++i)
      {
        result *= base;
      }
      return result;
    }
    /* Handle negative exponents */
    else if (e < 0)
    {
      /* Use the absolute value of the exponent for the loop */
      int positive_e = (e > 0) ? e : -e;
      for (i = 0; i < positive_e; ++i)
      {
        result *= base;
      }
      /* Check for division by zero */
      if (result == 0.0f)
      {
        /* Return an appropriate value for a math error */
        return 0.0f;
      }
      return 1.0f / result;
    }
    else /* e == 0 */
    {
      return 1.0f;
    }
  }

  /* Exponent is not an integer. Return -1.0f as originally implemented. */
  return -1.0f;
}

#endif /* CFD_MATH_H */

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
