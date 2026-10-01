#include "cache.h"
#include "fractal_fxpt.h"
#include "swap.h"
#include "vga.h"
#include <stddef.h>
#include <stdio.h>

#include <stdint.h>
// Constants describing the output device
const int SCREEN_WIDTH = 512;  //!< screen width
const int SCREEN_HEIGHT = 512; //!< screen height

// Constants describing the initial view port on the fractal function
const int FRAC_WIDTH = 3;  //!< default fractal width (3.0 in Q4.28)
const uint16_t N_MAX = 64; //!< maximum number of iterations

static rgb565 frameBuffer[512 * 512]; // not on the stack

int main() {
  volatile unsigned int *vga = (unsigned int *)0x50000020;
  volatile unsigned int reg, hi;

  fxpt_3_29 CX_0 = -2 * (1 << FXPT_FRAC_BITS);       // -2.0
  fxpt_3_29 CY_0 = -3 * (1 << (FXPT_FRAC_BITS - 1)); // -1.5
  fxpt_3_29 delta =
      (fxpt_3_29)(FRAC_WIDTH << (FXPT_FRAC_BITS - 9)); // 3/512 = 3<<20
  int i;
  vga_clear();
#ifdef __OR1300__
  /* enable the caches */
  icache_write_cfg(CACHE_DIRECT_MAPPED | CACHE_SIZE_8K | CACHE_REPLACE_FIFO);
  dcache_write_cfg(CACHE_FOUR_WAY | CACHE_SIZE_8K | CACHE_REPLACE_LRU |
                   CACHE_WRITE_BACK);
  icache_enable(1);
  dcache_enable(1);
#endif
  /* Enable the vga-controller's graphic mode */
  vga[0] = swap_u32(SCREEN_WIDTH);
  vga[1] = swap_u32(SCREEN_HEIGHT);
  vga[2] = swap_u32(1);
  vga[3] = swap_u32((unsigned int)&frameBuffer[0]);
  /* Clear screen */
  for (i = 0; i < SCREEN_WIDTH * SCREEN_HEIGHT; i++)
    frameBuffer[i] = 0;
  printf("Starting drawing a fractal\n");
  iter_to_colour_p c = 1 ? &iter_to_colour : *&iter_to_grayscale;
  draw_fractal(frameBuffer, SCREEN_WIDTH, SCREEN_HEIGHT,
               &calc_mandelbrot_point_soft, c, CX_0, CY_0, delta, N_MAX);
#ifdef __OR1300__
  dcache_flush();
#endif
  printf("Done\n");
} // Constants describing the output device
