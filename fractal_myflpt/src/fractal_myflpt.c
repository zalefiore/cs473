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

#define LNORM(m, e)                                                            \
  do {                                                                         \
    if (!((m) & HIDDEN)) {                                                     \
      myfloat d_ = 0;                                                          \
      if ((m) < (HIDDEN >> 7)) {                                               \
        (m) <<= 8;                                                             \
        d_ = ONE_EXP << 3;                                                     \
      }                                                                        \
      if ((m) < (HIDDEN >> 3)) {                                               \
        (m) <<= 4;                                                             \
        d_ += ONE_EXP << 2;                                                    \
      }                                                                        \
      if ((m) < (HIDDEN >> 1)) {                                               \
        (m) <<= 2;                                                             \
        d_ += ONE_EXP << 1;                                                    \
      }                                                                        \
      if (!((m) & HIDDEN)) {                                                   \
        (m) <<= 1;                                                             \
        d_ += ONE_EXP;                                                         \
      }                                                                        \
      if ((e) > d_)                                                            \
        (e) -= d_;                                                             \
      else                                                                     \
        (m) = (e) = 0;                                                         \
    }                                                                          \
  } while (0)

static inline myfloat calc_add_delta(myfloat a, myfloat d) {
  const myfloat de = d & EXP_MASK;
  myfloat ae = a & EXP_MASK;
  if (!ae)
    return d;

  myfloat am = HIDDEN | (a & MANT_MASK), dm = HIDDEN | (d & MANT_MASK), m,
          s = 0;

  if (!(a & SIGN_MASK)) {
    if (ae > de)
      dm >>= (ae - de) >> MANT_BITS;
    else if (de > ae) {
      am >>= (de - ae) >> MANT_BITS;
      ae = de;
    }
    m = am + dm;
    if (m >= (HIDDEN << 1)) {
      m >>= 1;
      ae += ONE_EXP;
    }
    return ae | (m - HIDDEN);
  }

  if (ae > de) {
    m = am - (dm >> ((ae - de) >> MANT_BITS));
    s = SIGN_MASK;
  } else if (de > ae) {
    m = dm - (am >> ((de - ae) >> MANT_BITS));
    ae = de;
  } else if (am > dm) {
    m = am - dm;
    s = SIGN_MASK;
  } else
    m = dm - am;

  LNORM(m, ae);
  return ae ? s | ae | (m - HIDDEN) : 0;
}

uint16_t calc_mandelbrot_point_soft(myfloat cx, myfloat cy, uint16_t n_max) {
  const myfloat cxs = cx & SIGN_MASK, cxe = cx & EXP_MASK,
                cxm = (cx & MANT_MASK) | (cxe ? HIDDEN : 0);
  const myfloat cys = cy & SIGN_MASK, cye = cy & EXP_MASK,
                cym = (cy & MANT_MASK) | (cye ? HIDDEN : 0);
  myfloat xs = cxs, xe = cxe, xm = cxm, ys = cys, ye = cye, ym = cym;
  myfloat pxs, pxe, pxm, pys, pye, pym;
  uint_fast16_t n = 0;

  while (n < n_max) {
    if (n & (n - 1)) {
      if (xm == pxm && ym == pym && xe == pxe && ye == pye && xs == pxs &&
          ys == pys)
        return n_max;
    } else {
      pxs = xs;
      pxe = xe;
      pxm = xm;
      pys = ys;
      pye = ye;
      pym = ym;
    }
    ++n;

    const uint32_t a = xm >> 12, b = ym >> 12;
    myfloat e, tm, ts, xxe = 0, xxm = 0, yye = 0, yym = 0;

    if (xe) {
      e = xe << 1;
      xxe = e >= MIN_EXP ? e - MIN_EXP : 0;
      xxm = a * a;
      if (xxm & 0x80000000u) {
        xxe += ONE_EXP;
        xxm >>= 4;
      } else
        xxm >>= 3;
    }
    if (ye) {
      e = ye << 1;
      yye = e >= MIN_EXP ? e - MIN_EXP : 0;
      yym = b * b;
      if (yym & 0x80000000u) {
        yye += ONE_EXP;
        yym >>= 4;
      } else
        yym >>= 3;
    }

    if (xxe > yye) {
      yym >>= (xxe - yye) >> MANT_BITS;
      yye = xxe;
    } else if (yye > xxe) {
      xxm >>= (yye - xxe) >> MANT_BITS;
      xxe = yye;
    }

    if (xxe >= FOUR - ONE_EXP && (xxe >= FOUR || xxm + yym >= (HIDDEN << 1)))
      break;

    myfloat xys = 0, xye = 0, xym = 0;
    if (xe && ye) {
      e = xe + ye;
      xys = xs ^ ys;
      xye = (e >= MIN_EXP ? e - MIN_EXP : 0) + ONE_EXP;
      xym = a * b;
      if (xym & 0x80000000u) {
        xye += ONE_EXP;
        xym >>= 4;
      } else
        xym >>= 3;
    }

    if (xxm >= yym) {
      tm = xxm - yym;
      ts = 0;
    } else {
      tm = yym - xxm;
      ts = SIGN_MASK;
    }
    LNORM(tm, xxe);

    if (!xxe) {
      xs = cxs;
      xe = cxe;
      xm = cxm;
    } else if (!cxe) {
      xs = ts;
      xe = xxe;
      xm = tm;
    } else {
      if (xxe > cxe)
        xm = cxm >> ((xxe - cxe) >> MANT_BITS);
      else {
        xm = cxm;
        if (cxe > xxe) {
          tm >>= (cxe - xxe) >> MANT_BITS;
          xxe = cxe;
        }
      }
      xe = xxe;
      if (ts == cxs) {
        xs = ts;
        xm += tm;
        if (xm >= (HIDDEN << 1)) {
          xm >>= 1;
          xe += ONE_EXP;
        }
      } else {
        if (tm >= xm) {
          xm = tm - xm;
          xs = ts;
        } else {
          xm -= tm;
          xs = cxs;
        }
        LNORM(xm, xe);
      }
    }

    if (!xye) {
      ys = cys;
      ye = cye;
      ym = cym;
    } else if (!cye) {
      ys = xys;
      ye = xye;
      ym = xym;
    } else {
      if (xye > cye)
        ym = cym >> ((xye - cye) >> MANT_BITS);
      else {
        ym = cym;
        if (cye > xye) {
          xym >>= (cye - xye) >> MANT_BITS;
          xye = cye;
        }
      }
      ye = xye;
      if (xys == cys) {
        ys = xys;
        ym += xym;
        if (ym >= (HIDDEN << 1)) {
          ym >>= 1;
          ye += ONE_EXP;
        }
      } else {
        if (xym >= ym) {
          ym = xym - ym;
          ys = xys;
        } else {
          ym -= xym;
          ys = cys;
        }
        LNORM(ym, ye);
      }
    }
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
      cx = calc_add_delta(cx, delta);
    }
    cy = calc_add_delta(cy, delta); /* was: cy += ... */
  }
}
