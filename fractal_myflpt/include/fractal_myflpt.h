#ifndef FRACTAL_MYFLPT_H
#define FRACTAL_MYFLPT_H

#include <stdint.h>

//! Colour type (5-bit red, 6-bit green, 5-bit blue)
typedef uint16_t rgb565;
typedef uint32_t myfloat;

myfloat to_myfloat(float n);
//! \brief Pointer to fractal point calculation function
typedef uint16_t (*calc_frac_point_p)(myfloat cx, myfloat cy, uint16_t n_max);

uint16_t calc_mandelbrot_point_soft(myfloat cx, myfloat cy, uint16_t n_max);

//! Pointer to function mapping iteration to colour value
typedef rgb565 (*iter_to_colour_p)(uint16_t iter, uint16_t n_max);

rgb565 iter_to_bw(uint16_t iter, uint16_t n_max);
rgb565 iter_to_grayscale(uint16_t iter, uint16_t n_max);
rgb565 iter_to_colour(uint16_t iter, uint16_t n_max);

void draw_fractal(rgb565 *fbuf, int width, int height, calc_frac_point_p cfp_p,
                  iter_to_colour_p i2c_p, myfloat cx_0, myfloat cy_0,
                  myfloat delta, uint16_t n_max);

#endif // FRACTAL_MYFLPT_H
