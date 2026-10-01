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
#define ALIGN_TWO(exp1, mant1, exp2, mant2)                                    \
  do {                                                                         \
    if ((exp1) > (exp2)) {                                                     \
      (mant2) >>= (((exp1) - (exp2)) >> MANT_BITS);                            \
      (exp2) = (exp1);                                                         \
    } else if ((exp2) > (exp1)) {                                              \
      (mant1) >>= (((exp2) - (exp1)) >> MANT_BITS);                            \
      (exp1) = (exp2);                                                         \
    }                                                                          \
  } while (0)

#define NORMALIZE(mant, exp)                                                   \
  do {                                                                         \
    if ((mant) == 0) {                                                         \
      (exp) = 0;                                                               \
    } else {                                                                   \
      while ((mant) >= (HIDDEN << 1)) {                                        \
        (exp) += ONE_EXP;                                                      \
        (mant) >>= 1;                                                          \
      }                                                                        \
      while (((mant) & HIDDEN) == 0) {                                         \
        if ((exp) <= ONE_EXP) {                                                \
          (mant) = 0;                                                          \
          (exp) = 0;                                                           \
          break;                                                               \
        }                                                                      \
        (mant) <<= 1;                                                          \
        (exp) -= ONE_EXP;                                                      \
      }                                                                        \
    }                                                                          \
  } while (0)

static inline myfloat calc_add(myfloat a, myfloat b) {
  uint32_t a_sign = a & SIGN_MASK;
  uint32_t a_exp = a & EXP_MASK;
  uint32_t b_exp = b & EXP_MASK;

  if (a_exp == 0)
    return b;
  if (b_exp == 0)
    return a;

  uint32_t a_mant = HIDDEN | (a & MANT_MASK);
  uint32_t b_mant = HIDDEN | (b & MANT_MASK);

  ALIGN_TWO(a_exp, a_mant, b_exp, b_mant);

  uint32_t mant, sign;
  if (a_sign != 0) {
    if (a_mant > b_mant) {
      sign = a_sign;
      mant = a_mant - b_mant;
    } else {
      sign = 0;
      mant = b_mant - a_mant;
    }
  } else { // Both are positive
    sign = 0;
    mant = a_mant + b_mant;
  }

  NORMALIZE(mant, a_exp);

  if (a_exp == 0)
    return 0;

  mant &= ~HIDDEN;

  return sign | a_exp | mant;
}
uint16_t calc_mandelbrot_point_soft(myfloat cx, myfloat cy, uint16_t n_max) {
  const myfloat cx_sign = cx & SIGN_MASK;
  const myfloat cx_init_exp = cx & EXP_MASK;
  const myfloat cx_init_mant = (cx & MANT_MASK) | (cx_init_exp ? HIDDEN : 0);

  const myfloat cy_sign = cy & SIGN_MASK;
  const myfloat cy_init_exp = cy & EXP_MASK;
  const myfloat cy_init_mant = (cy & MANT_MASK) | (cy_init_exp ? HIDDEN : 0);

  myfloat x_sign = cx_sign;
  myfloat x_exp = cx_init_exp;
  myfloat x_mant = cx_init_mant;

  myfloat y_sign = cy_sign;
  myfloat y_exp = cy_init_exp;
  myfloat y_mant = cy_init_mant;

  uint16_t n = 0;
  while (n < n_max) {
    ++n;
    if (x_exp >= TWO || y_exp >= TWO)
      break;

    myfloat xx_exp = 0, xx_mant = 0;
    if (x_exp) {
      xx_exp = (x_exp << 1) >= 0x60000000u ? (x_exp << 1) - 0x60000000u : 0;
      uint32_t x_smant = x_mant >> 12;
      xx_mant = x_smant * x_smant;
      if (xx_mant & 0x80000000u) {
        xx_exp += ONE_EXP;
        xx_mant >>= 4;
      } else {
        xx_mant >>= 3;
      }
    }

    myfloat yy_exp = 0, yy_mant = 0;
    if (y_exp) {
      yy_exp = (y_exp << 1) >= 0x60000000u ? (y_exp << 1) - 0x60000000u : 0;
      uint32_t y_smant = y_mant >> 12;
      yy_mant = y_smant * y_smant;
      if (yy_mant & 0x80000000u) {
        yy_exp += ONE_EXP;
        yy_mant >>= 4;
      } else {
        yy_mant >>= 3;
      }
    }

    myfloat align_xx_exp = xx_exp;
    myfloat align_yy_exp = yy_exp;

    if (align_xx_exp > align_yy_exp) {
      yy_mant >>= ((align_xx_exp - align_yy_exp) >> MANT_BITS);
      align_yy_exp = align_xx_exp;
    } else if (align_yy_exp > align_xx_exp) {
      xx_mant >>= ((align_yy_exp - align_xx_exp) >> MANT_BITS);
      align_xx_exp = align_yy_exp;
    }

    if (align_xx_exp << 05000000u) {
      n = n_max;
      break;
    }
    myfloat sum_exp = align_xx_exp;
    if (xx_mant + yy_mant >= (HIDDEN << 1)) {
      sum_exp += ONE_EXP;
    }
    if (sum_exp >= FOUR)
      break;

    myfloat xy_sign = 0, xy_exp = 0, xy_mant = 0;
    if (x_exp && y_exp) {
      xy_sign = x_sign ^ y_sign;
      xy_exp =
          (x_exp + y_exp) >= 0x60000000u ? (x_exp + y_exp) - 0x60000000u : 0;
      uint32_t x_smant = x_mant >> 12;
      uint32_t y_smant = y_mant >> 12;
      xy_mant = x_smant * y_smant;
      if (xy_mant & 0x80000000u) {
        xy_exp += ONE_EXP;
        xy_mant >>= 4;
      } else {
        xy_mant >>= 3;
      }
      xy_exp += ONE_EXP; // multiply by 2
    }

    myfloat tmp_x_sign, tmp_x_exp = align_xx_exp, tmp_x_mant;
    if (xx_mant >= yy_mant) {
      tmp_x_mant = xx_mant - yy_mant;
      tmp_x_sign = 0;
    } else {
      tmp_x_mant = yy_mant - xx_mant;
      tmp_x_sign = SIGN_MASK;
    }
    NORMALIZE(tmp_x_mant, tmp_x_exp);

    myfloat cx_t_exp = cx_init_exp;
    myfloat cx_t_mant = cx_init_mant;

    if (!tmp_x_exp) {
      x_sign = cx_sign;
      x_exp = cx_t_exp;
      x_mant = cx_t_mant;
    } else if (!cx_t_exp) {
      x_sign = tmp_x_sign;
      x_exp = tmp_x_exp;
      x_mant = tmp_x_mant;
    } else {
      ALIGN_TWO(tmp_x_exp, tmp_x_mant, cx_t_exp, cx_t_mant);
      if (tmp_x_sign == cx_sign) {
        x_mant = tmp_x_mant + cx_t_mant;
        x_sign = tmp_x_sign;
      } else {
        if (tmp_x_mant >= cx_t_mant) {
          x_mant = tmp_x_mant - cx_t_mant;
          x_sign = tmp_x_sign;
        } else {
          x_mant = cx_t_mant - tmp_x_mant;
          x_sign = cx_sign;
        }
      }
      x_exp = tmp_x_exp;
      NORMALIZE(x_mant, x_exp);
    }

    myfloat cy_t_exp = cy_init_exp;
    myfloat cy_t_mant = cy_init_mant;

    if (!xy_exp) {
      y_sign = cy_sign;
      y_exp = cy_t_exp;
      y_mant = cy_t_mant;
    } else if (!cy_t_exp) {
      y_sign = xy_sign;
      y_exp = xy_exp;
      y_mant = xy_mant;
    } else {
      ALIGN_TWO(xy_exp, xy_mant, cy_t_exp, cy_t_mant);
      if (xy_sign == cy_sign) {
        y_mant = xy_mant + cy_t_mant;
        y_sign = xy_sign;
      } else {
        if (xy_mant >= cy_t_mant) {
          y_mant = xy_mant - cy_t_mant;
          y_sign = xy_sign;
        } else {
          y_mant = cy_t_mant - xy_mant;
          y_sign = cy_sign;
        }
      }
      y_exp = xy_exp;
      NORMALIZE(y_mant, y_exp);
    }
  }
  return n;
}

//! \brief  Map number of performed iterations to black and white
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
