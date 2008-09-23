/** 
 * @file   rta_onepole.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Fri Aug 29 12:38:46 2008
 * 
 * @brief  One-pole one-zero filters
 * 
 * Simple low-pass and high-pass filters.
 * \see rta_biquad.h
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */

#include "rta_onepole.h"

inline rta_real_t rta_onepole_lowpass(const rta_real_t x, const rta_real_t f0,
                                      rta_real_t * state)
{
  *state = x * f0 + *state * (1. - f0);
  return *state;
}

inline rta_real_t rta_onepole_highpass(const rta_real_t x, const rta_real_t f0,
                                       rta_real_t * state)
{
  /* highpass = x - lowpass */

  rta_real_t y = f0 * x + *state;
  *state = (1. - f0) * y;
  return (x - y);
}

void rta_onepole_lowpass_vector(
  rta_real_t * y,
  const rta_real_t * x, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state)
{
  unsigned int i;

  for(i=0; i<x_size; i++)
  {
    y[i] = rta_onepole_lowpass(x[i], f0, state);
  }

  return;
}

void rta_onepole_lowpass_vector_stride(
  rta_real_t * y, const int y_stride,
  const rta_real_t * x, const int x_stride, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state)
{
  int ix, iy;

  for(ix = 0, iy = 0;
      ix < x_size*x_stride;
      ix += x_stride, iy += y_stride)
  {
    y[iy] = rta_onepole_lowpass(x[ix], f0, state);
  }

  return;
}

void rta_onepole_highpass_vector(
  rta_real_t * y,
  const rta_real_t * x, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state)
{
  unsigned int i;

  for(i=0; i<x_size; i++)
  {
    y[i] = rta_onepole_highpass(x[i], f0, state);
  }

  return;
}

void rta_onepole_highpass_vector_stride(
  rta_real_t * y, const int y_stride,
  const rta_real_t * x, const int x_stride, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state)
{
  int ix, iy;

  for(ix = 0, iy = 0;
      ix < x_size*x_stride;
      ix += x_stride, iy += y_stride)
  {
    y[iy] = rta_onepole_highpass(x[ix], f0, state);
  }

  return;
}
