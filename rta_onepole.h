/** 
 * @file   rta_onepole.h
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

#ifndef _RTA_ONEPOLE_H_
#define _RTA_ONEPOLE_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * One-pole low-pass filter computed as:
 * y(n) = f0 * x(n) - (f0 - 1) * y(n-1)
 * \see rta_onepole_highpass
 * 
 * @param x is an input sample
 * @param f0 is the cutoff frequency, normalised by the nyquist frequency.
 * @param state is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 * 
 * @return the output sample y
 */
inline rta_real_t rta_onepole_lowpass(rta_real_t x, const rta_real_t f0,
                                      rta_real_t * state);

/** 
 * One-pole high-pass filter computed as the difference between the
 * input and a low-pass filtered input:
 * y(n) = x(n) - ( f0 * x(n) - (f0 - 1) * y(n-1) )
 * \see rta_onepole_lowpass
 * 
 * @param x is an input sample
 * @param f0 is the cutoff frequency, normalised by the nyquist frequency.
 * @param state is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 * 
 * @return the output sample y
 */
inline rta_real_t rta_onepole_highpass(rta_real_t x, const rta_real_t f0,
                                       rta_real_t * state);
/** 
 * One-pole low-pass computation on a vector of samples.
 * \see rta_onepole_lowpass
 * 
 * @param y is a vector of output samples. Its size is 'x_size'
 * @param x is a vector of input samples. Its size is 'x_size'
 * @param x_size is the size of 'y' and 'x'
 * @param f0 is the cutoff frequency, normalised by the nyquist frequency.
 * @param state is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 */
void rta_onepole_lowpass_vector(
  rta_real_t * y,
  const rta_real_t * x, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state);

/** 
 * One-pole low-pass computation on a vector of samples.
 * \see rta_onepole_lowpass
 * 
 * @param y is a vector of output samples. Its size is 'x_size'
 * @param y_stride is 'y' stride
 * @param x is a vector of input samples. Its size is 'x_size'
 * @param x_stride is 'x' stride
 * @param x_size is the size of 'y' and 'x'
 * @param f0 is the cutoff frequency, normalised by the nyquist frequency.
 * @param state is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 */
void rta_onepole_lowpass_vector_stride(
  rta_real_t * y, const int y_stride,
  const rta_real_t * x, const int x_stride, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state);

/** 
 * One-pole high-pass computation on a vector of samples.
 * \see rta_onepole_highpass
 * 
 * @param y is a vector of output samples. Its size is 'x_size'
 * @param x is a vector of input samples. Its size is 'x_size'
 * @param x_size is the size of 'y' and 'x'
 * @param f0 is the cutoff frequency, normalised by the nyquist frequency.
 * @param state is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 */
void rta_onepole_highpass_vector(
  rta_real_t * y,
  const rta_real_t * x, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state);

/** 
 * One-pole high-pass computation on a vector of samples.
 * \see rta_onepole_highpass
 * 
 * @param y is a vector of output samples. Its size is 'x_size'
 * @param y_stride is 'y' stride
 * @param x is a vector of input samples. Its size is 'x_size'
 * @param x_stride is 'x' stride
 * @param x_size is the size of 'y' and 'x'
 * @param f0 is the cutoff frequency, normalised by the nyquist frequency.
 * @param state is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 */
void rta_onepole_highpass_vector_stride(
  rta_real_t * y, const int y_stride,
  const rta_real_t * x, const int x_stride, const unsigned int x_size, 
  const rta_real_t f0, rta_real_t * state);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_ONEPOLE_H_ */
