#include "fractal_myflpt.h"
#include <assert.h>
#include <stdint.h>
#include <stdio.h>
#include <swap.h>

//! \brief  Mandelbrot fractal point calculation function
//! \param  cx    x-coordinate
//! \param  cy    y-coordinate
//! \param  n_max maximum number of iterations
//! \return       number of performed iterations at coordinate (cx, cy)
//
typedef uint32_t myfloat;
#define BIAS 12
#define HIDDEN 0x08000000u
#define MANT_BITS 27
#define MANT_MASK 0x07FFFFFFu
#define EXP_BITS 4
#define DIFF 0x39800000u
#define EXP_MASK 0x78000000u
#define SIGN_MASK 0x80000000u
#define SHIFT 13
#define ONE_EXP 0x08000000u
#define MAX_NEG 0xFFFFFFFFu
#define MAX_POS 0xFFFFFFFFu
#define TWO 0x68000000u
#define FOUR 0x70000000u
myfloat to_myfloat(float n) {
  union {
    float f;
    uint32_t u;
  } c = {.f = n};

  uint32_t sign = c.u & SIGN_MASK;
  uint32_t mag = c.u & 0x7FFFFFFFu;
  return sign + (((mag - DIFF) << 4) & 0x7FFFFFFFu); /* shift by 4 */
}

myfloat calc_add(myfloat a, myfloat b) {
  uint32_t signa = a & SIGN_MASK;
  uint32_t signb = b & SIGN_MASK;
  uint32_t expa = a & EXP_MASK;
  uint32_t expb = b & EXP_MASK;
  if (expa == 0)
    return b;
  if (expb == 0)
    return a;
  uint32_t manta = HIDDEN + (a & MANT_MASK);
  uint32_t mantb = HIDDEN + (b & MANT_MASK);

  if (expa > expb) {
    mantb >>= ((expa - expb) >> MANT_BITS);
  } else if (expb > expa) {
    manta >>= ((expb - expa) >> MANT_BITS);
    expa = expb;
  }
  uint32_t mant, sign;
  if (signa != signb) {
    if (manta > mantb) {
      sign = signa;
      mant = manta - mantb;
    } else {
      mant = mantb - manta;
      sign = signb;
    }
    if (mant == 0)
      return 0;
    while ((mant & HIDDEN) == 0) {
      if (expa <= ONE_EXP) /* would become code 0 */
        return 0;
      mant <<= 1;
      expa -= ONE_EXP;
    }
  } else {
    sign = signa;
    mant = manta + mantb;
    if (mant & (HIDDEN << 1)) {
      expa += ONE_EXP;
      mant >>= 1;
    }
  }
  mant -= HIDDEN;
  return sign + expa + mant;
}

uint16_t calc_mandelbrot_point_soft(myfloat cx, myfloat cy, uint16_t n_max) {
  myfloat x_sign = cx & SIGN_MASK;
  myfloat x_exp = cx & EXP_MASK;
  myfloat x_mant = (cx & MANT_MASK) | (x_exp ? HIDDEN : 0);
  myfloat y_sign = cy & SIGN_MASK;
  myfloat y_exp = cy & EXP_MASK;
  myfloat y_mant = (cy & MANT_MASK) | (y_exp ? HIDDEN : 0);
  uint16_t n = 0;
  myfloat xx_exp, yy_exp, xx_mant, yy_mant, xy_sign, xy_exp, xy_mant;
  while (n < n_max) {
    myfloat cx_sign = cx & SIGN_MASK;
    myfloat cx_exp = cx & EXP_MASK;
    myfloat cx_mant = (cx & MANT_MASK) | (cx_exp ? HIDDEN : 0);

    myfloat cy_sign = cy & SIGN_MASK;
    myfloat cy_exp = cy & EXP_MASK;
    myfloat cy_mant = (cy & MANT_MASK) | (cy_exp ? HIDDEN : 0);
    ++n;
    if (x_exp >= TWO || y_exp >= TWO)
      break;
    // Start doing moltiplication for x
    if (x_exp == 0) {
      x_exp = 0;
      x_mant = 0;
      xx_exp = 0; // GOTO SOMEWHERE
      xx_mant = 0;
      goto end_of_xx;
    }

    xx_exp = (x_exp << 1) >= 0x60000000u ? (x_exp << 1) - 0x60000000u : 0;

    uint32_t x_smant = x_mant >> 12;
    xx_mant = x_smant * x_smant;
    if (xx_mant & 0x80000000u) {
      xx_exp += ONE_EXP;
      xx_mant >>= 4;
    } else {
      xx_mant >>= 3;
    }
  end_of_xx:
    // Start doing multiplicatin for y
    if (y_exp == 0) {
      y_exp = 0;
      y_mant = 0; // GOTO SOMEWHERE
      yy_exp = 0;
      yy_mant = 0;
      goto end_of_yy;
    }

    yy_exp = (y_exp << 1) >= 0x60000000u ? (y_exp << 1) - 0x60000000u : 0;

    uint32_t y_smant = y_mant >> 12;
    yy_mant = y_smant * y_smant;
    if (yy_mant & 0x80000000u) {
      yy_exp += ONE_EXP;
      yy_mant >>= 4;
    } else {
      yy_mant >>= 3;
    }
  end_of_yy:

    if (x_exp == 0 || y_exp == 0)
      goto one_zero;

    if (xx_exp > yy_exp) {
      yy_mant >>= ((xx_exp - yy_exp) >> MANT_BITS);
      yy_exp = xx_exp;
    } else if (yy_exp > xx_exp) {
      xx_mant >>= ((yy_exp - xx_exp) >> MANT_BITS);
      xx_exp = yy_exp;
    }

    myfloat xx_yy_mant = xx_mant + yy_mant;
    myfloat sum_exp = xx_exp;
    if (xx_yy_mant & (HIDDEN << 1))
      sum_exp += ONE_EXP;
    if (sum_exp >= FOUR)
      break; // xx = calc_mult(x, x);
    // yy = calc_mult(y, y);

    // Working on xy
    xy_sign = x_sign ^ y_sign;

    xy_exp = (x_exp + y_exp) >= 0x60000000u ? (x_exp + y_exp) - 0x60000000u : 0;

    xy_mant = x_smant * y_smant;
    if (xy_mant & 0x80000000u) {
      xy_exp += ONE_EXP;
      xy_mant >>= 4;
    } else {
      xy_mant >>= 3;
    }
    xy_exp += ONE_EXP; // multiply by 2

    // ADDITION PART
    // calculation new x
    if (cx_exp > yy_exp) {
      yy_mant >>= ((cx_exp - yy_exp) >> MANT_BITS);
      xx_mant >>= ((cx_exp - xx_exp) >> MANT_BITS);
      yy_exp = cx_exp;
      xx_exp = cx_exp;
    } else if (yy_exp > cx_exp) {
      cx_mant >>= ((yy_exp - cx_exp) >> MANT_BITS);
      cx_exp = yy_exp;
    }

    x_exp = cx_exp;
    if (cx_sign == 0) {
      if (cx_mant + xx_mant < yy_mant) {
        x_sign = SIGN_MASK;
        x_mant = -cx_mant - xx_mant + yy_mant;
      } else if (cx_mant + xx_mant > yy_mant) {
        x_sign = 0;
        x_mant = xx_mant + cx_mant - yy_mant;
      } else {
        x_sign = 0;
        x_mant = 0;
        x_exp = 0;
      }

    } else {
      if (xx_mant > cx_mant + yy_mant) {
        x_sign = 0;
        x_mant = xx_mant - cx_mant - yy_mant;
      } else if (xx_mant < cx_mant + yy_mant) {
        x_sign = SIGN_MASK;
        x_mant = -xx_mant + cx_mant + yy_mant;
      } else {
        x_sign = 0;
        x_mant = 0;
        x_exp = 0;
      }
    };

    if (x_mant & (HIDDEN << 1)) {
      x_exp += ONE_EXP;
      x_mant >>= 1;
    }
    while ((x_mant & HIDDEN) == 0) {
      if (x_exp <= ONE_EXP) { /* would become code 0 */
        x_mant = 0;
        x_exp = 0;
        break;
      }
      x_mant <<= 1;
      x_exp -= ONE_EXP;
    }

    // calculation new y
    if (cy_exp > xy_exp) {
      xy_mant >>= ((cy_exp - xy_exp) >> MANT_BITS);
      xy_exp = cy_exp;
    } else if (xy_exp > cy_exp) {
      cy_mant >>= ((xy_exp - cy_exp) >> MANT_BITS);
      cy_exp = xy_exp;
    }

    y_exp = cy_exp;

    if (cy_sign ^ xy_sign) {
      if (cy_mant < xy_mant) {
        y_sign = xy_sign;
        y_mant = xy_mant - cy_mant;
      } else if (cy_mant > xy_mant) {
        y_sign = cy_sign;
        y_mant = cy_mant - xy_mant;
      } else {
        y_exp = 0;
        y_mant = 0;
      }
    } else {
      y_sign = cy_sign;
      y_mant = xy_mant + cy_mant;
    };

    if (y_mant & (HIDDEN << 1)) {
      y_exp += ONE_EXP;
      y_mant >>= 1;
    }
    while ((y_mant & HIDDEN) == 0) {
      if (y_exp <= ONE_EXP) { /* would become code 0 */
        y_mant = 0;
        y_exp = 0;
        break;
      }
      y_mant <<= 1;
      y_exp -= ONE_EXP;
    }

    continue;
  one_zero:
    xy_exp = 0;
    xy_mant = 0;

    if (cx_exp == 0) {
      x_mant = xx_mant + yy_mant;
      x_exp = xx_exp + yy_exp;
      x_sign = yy_exp ? SIGN_MASK : 0;
      goto skip_x;
    }

    if (xx_exp) {
      if (cx_exp > xx_exp) {
        xx_mant >>= ((cx_exp - xx_exp) >> MANT_BITS);
        xx_exp = cx_exp;
      } else if (xx_exp > cx_exp) {
        cx_mant >>= ((xx_exp - cx_exp) >> MANT_BITS);
        cx_exp = xx_exp;
      }
      x_exp = cx_exp;

      if (cx_sign == 0) {
        x_sign = 0;
        x_mant = xx_mant + cx_mant;
      } else {
        if (xx_mant > cx_mant) {
          x_sign = 0;
          x_mant = xx_mant - cx_mant;
        } else if (xx_mant < cx_mant) {
          x_sign = 1;
          x_mant = cx_mant - xx_mant;
        } else {
          x_exp = 0;
          x_mant = 0;
          x_sign = 0;
        }
      }
    } else if (yy_exp) {
      if (cx_exp > yy_exp) {
        yy_mant >>= ((cx_exp - yy_exp) >> MANT_BITS);
        yy_exp = cx_exp;
      } else if (yy_exp > cx_exp) {
        cx_mant >>= ((yy_exp - cx_exp) >> MANT_BITS);
        cx_exp = yy_exp;
      }
      x_exp = cx_exp;

      if (cx_sign != 0) {
        x_sign = 1;
        x_mant = yy_mant + cx_mant;
      } else {
        if (yy_mant > cx_mant) {
          x_sign = 1;
          x_mant = yy_mant - cx_mant;
        } else if (yy_mant < cx_mant) {
          x_sign = 0;
          x_mant = cx_mant - yy_mant;
        } else {
          x_exp = 0;
          x_mant = 0;
          x_sign = 0;
          goto skip_x;
        }
      }
    } else {
      x_sign = cx;
      x_exp = cx_exp;
      x_mant = cx_mant;
    }

    if (x_mant & (HIDDEN << 1)) {
      x_exp += ONE_EXP;
      x_mant >>= 1;
    }
    while ((x_mant & HIDDEN) == 0) {
      if (x_exp <= ONE_EXP) { /* would become code 0 */
        x_mant = 0;
        x_exp = 0;
        break;
      }
      x_mant <<= 1;
      x_exp -= ONE_EXP;
    }
  skip_x:
    y_sign = cy_sign;
    y_exp = cy_exp;
    y_mant = cy_mant;
  }
  return n;
} //! \brief  Map number of performed iterations to black and white
//! \param  iter  performed number of iterations
//! \param  n_max maximum number of iterations
//! \return       colour
rgb565 iter_to_bw(uint16_t iter, uint16_t n_max) {
  if (iter == n_max) {
    return 0x0000;
  }
  return 0xffff;
}

//! \brief  Map number of performed iterations to grayscale
//! \param  iter  performed number of iterations
//! \param  n_max maximum number of iterations
//! \return       colour
rgb565 iter_to_grayscale(uint16_t iter, uint16_t n_max) {
  if (iter == n_max) {
    return 0x0000;
  }
  uint16_t brightness = iter & 0xf;
  return swap_u16(((brightness << 12) | ((brightness << 7) | brightness << 1)));
}

//! \brief Calculate binary logarithm for unsigned integer argument x
//! \note  For x equal 0, the function returns -1.
int ilog2(unsigned x) {
  if (x == 0)
    return -1;
  int n = 1;
  if ((x >> 16) == 0) {
    n += 16;
    x <<= 16;
  }
  if ((x >> 24) == 0) {
    n += 8;
    x <<= 8;
  }
  if ((x >> 28) == 0) {
    n += 4;
    x <<= 4;
  }
  if ((x >> 30) == 0) {
    n += 2;
    x <<= 2;
  }
  n -= x >> 31;
  return 31 - n;
}

//! \brief  Map number of performed iterations to a colour
//! \param  iter  performed number of iterations
//! \param  n_max maximum number of iterations
//! \return colour in rgb565 format little Endian (big Endian for openrisc)
rgb565 iter_to_colour(uint16_t iter, uint16_t n_max) {
  if (iter == n_max) {
    return 0x0000;
  }
  uint16_t brightness = (iter & 1) << 4 | 0xF;
  uint16_t r = (iter & (1 << 3)) ? brightness : 0x0;
  uint16_t g = (iter & (1 << 2)) ? brightness : 0x0;
  uint16_t b = (iter & (1 << 1)) ? brightness : 0x0;
  return swap_u16(((r & 0x1f) << 11) | ((g & 0x1f) << 6) | ((b & 0x1f)));
}

rgb565 iter_to_colour1(uint16_t iter, uint16_t n_max) {
  if (iter == n_max) {
    return 0x0000;
  }
  uint16_t brightness = ((iter & 0x78) >> 2) ^ 0x1F;
  uint16_t r = (iter & (1 << 2)) ? brightness : 0x0;
  uint16_t g = (iter & (1 << 1)) ? brightness : 0x0;
  uint16_t b = (iter & (1 << 0)) ? brightness : 0x0;
  return swap_u16(((r & 0xf) << 12) | ((g & 0xf) << 7) | ((b & 0xf) << 1));
}

//! \brief  Draw fractal into frame buffer
//! \param  width  width of frame buffer
//! \param  height height of frame buffer
//! \param  cfp_p  pointer to fractal function
//! \param  i2c_p  pointer to function mapping number of iterations to
//! colour
//! \param  cx_0   start x-coordinate
//! \param  cy_0   start y-coordinate
//! \param  delta  increment for x- and y-coordinate
//! \param  n_max  maximum number of iterations
void draw_fractal(rgb565 *fbuf, int width, int height, calc_frac_point_p cfp_p,
                  iter_to_colour_p i2c_p, myfloat cx_0, myfloat cy_0,
                  myfloat delta, uint16_t n_max) {
  rgb565 *pixel = fbuf;
  myfloat cy = cy_0;
  for (int k = 0; k < height; ++k) {
    myfloat cx = cx_0;
    for (int i = 0; i < width; ++i) {
      uint16_t n_iter = (*cfp_p)(cx, cy, n_max);
      *(pixel++) = (*i2c_p)(n_iter, n_max);
      cx = calc_add(cx, delta);
    }
    cy = calc_add(cy, delta); /* was: cy += ... */
  }
}
