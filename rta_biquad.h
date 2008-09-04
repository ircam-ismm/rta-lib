/** 
 * @file   rta_biquad.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Fri Aug 29 12:38:46 2008
 *
 * @brief  Biquad filter and coefficients calculations.
 *
 * Based on the "Cookbook formulae for audio EQ biquad filter
 * coefficients" by Robert Bristow-Johnson.
 *
 * \htmlonly <pre> 
 * y(n) = b0 x(n) + b1 x(n-1) + b2 x(n-2) 
 *                - a1 x(n-1) - a2 x(n-2) 
 * </pre> \endhtmlonly
 * 
 * (This is Matlab convention, MaxMSP biquad~ swaps the names for a
 * and b.)
 *
 * a0 is always 1. as each coefficient is normalised by a0, including
 * a0.
 *
 * For every function, a1 is a[0] and a2 is a[1]. b0 is b[0], b1 is
 * b[1] and b2 is b[2].
 *
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 *
 */

#ifndef _RTA_BIQUAD_H_
#define _RTA_BIQUAD_H_ 1

#include "rta.h"
#include "rta_filter.h" /* filter types */

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * Biquad coefficients for a low-pass filter.
 * H(s) = 1 / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency: the filter is closed if f0 == 0. and open if f0 == 1.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_lowpass_coefs(rta_real_t * b, rta_real_t * a,
                              const rta_real_t f0, const rta_real_t q);

/**  
 * Biquad coefficients for a low-pass filter.
 * H(s) = 1 / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency: the filter is closed if f0 == 0. and open if f0 == 1.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_lowpass_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0, const rta_real_t q);

/** 
 * Biquad coefficients for a high-pass filter.
 * H(s) = s^2 / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency: the filter is closed if f0 == 1. and open if f0 == 0.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_highpass_coefs(rta_real_t * b, rta_real_t * a,
                               const rta_real_t f0, const rta_real_t q);

/** 
 * Biquad coefficients for a high-pass filter.
 * H(s) = s^2 / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency: the filter is closed if f0 == 1. and open if f0 == 0.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_highpass_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0, const rta_real_t q);

/** 
 * Biquad coefficients for a band-pass filter with constant skirt. The
 * peak gain is 'q'.
 * H(s) = s / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_bandpass_constant_skirt_coefs(rta_real_t * b, rta_real_t * a,
                                              const rta_real_t f0,
                                              const rta_real_t q);

/** 
 * Biquad coefficients for a band-pass filter with constant skirt. The
 * peak gain is 'q'.
 * H(s) = s / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_bandpass_constant_skirt_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0,
  const rta_real_t q);

/** 
 * Biquad coefficients for a band-pass filter with constant 0 dB peak.
 * H(s) = (s/Q) / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_bandpass_constant_peak_coefs(rta_real_t * b, rta_real_t * a,
                                             const rta_real_t f0,
                                             const rta_real_t q);

/** 
 * Biquad coefficients for a band-pass filter with constant 0 dB peak.
 * H(s) = (s/Q) / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_bandpass_constant_peak_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0,
  const rta_real_t q);

/** 
 * Biquad coefficients for a notch filter.
 * H(s) = (s^2 + 1) / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_notch_coefs(rta_real_t * b, rta_real_t * a,
                            const rta_real_t f0, const rta_real_t q);

/** 
 * Biquad coefficients for a notch filter.
 * H(s) = (s^2 + 1) / (s^2 + s/q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_notch_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0, const rta_real_t q);

/** 
 * Biquad coefficients for an all-pass filter. 
 * H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_allpass_coefs(rta_real_t * b, rta_real_t * a,
                              const rta_real_t f0, const rta_real_t q);

/** 
 * Biquad coefficients for an all-pass filter. 
 * H(s) = (s^2 - s/Q + 1) / (s^2 + s/Q + 1)
 *
 * @param b is a vector of feed-forward coefficients. To apply a
 * (linear) gain, simply multiply the b coefficients by the gain.
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 */
void rta_biquad_allpass_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0, const rta_real_t q);

/** 
 * Biquad coefficients for an peaking filter.
 * H(s) = (s^2 + s*(g/q) + 1) / (s^2 + s/(g*q) + 1),
 * g = sqrt('gain'),
 * 'gain' is linear.
 *
 * @param b is a vector of feed-forward coefficients
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_peaking_coefs(rta_real_t * b, rta_real_t * a,
                              const rta_real_t f0, const rta_real_t q,
                              const rta_real_t gain);

/** 
 * Biquad coefficients for an peaking filter.
 * H(s) = (s^2 + s*(g/q) + 1) / (s^2 + s/(g*q) + 1),
 * g = sqrt('gain'),
 * 'gain' is linear.
 *
 * @param b is a vector of feed-forward coefficients
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_peaking_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0, const rta_real_t q,
  const rta_real_t gain);

/** 
 * Biquad coefficients for an low-shelf filter.
 * H(s) = g * (s^2 + (sqrt(g)/q)*s + g)/(g*s^2 + (sqrt(g)/q)*s + 1)
 * g = sqrt('gain'),
 * 'gain' is linear.
 *
 * @param b is a vector of feed-forward coefficients
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_lowshelf_coefs(rta_real_t * b, rta_real_t * a,
                               const rta_real_t f0, const rta_real_t q,
                               const rta_real_t gain);

/** 
 * Biquad coefficients for an low-shelf filter.
 * H(s) = g * (s^2 + (sqrt(g)/q)*s + g)/(g*s^2 + (sqrt(g)/q)*s + 1)
 * g = sqrt('gain'),
 * 'gain' is linear.
 *
 * @param b is a vector of feed-forward coefficients
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_lowshelf_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0, const rta_real_t q,
  const rta_real_t gain);

/** 
 * Biquad coefficients for an high-shelf filter.
 * H(s) = g * (g*s^2 + (sqrt(g)/q)*s + 1)/(s^2 + (sqrt(g)/q)*s + g)
 * g = sqrt('gain'),
 * 'gain' is linear.
 *
 * @param b is a vector of feed-forward coefficients
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_highshelf_coefs(rta_real_t * b, rta_real_t * a,
                                const rta_real_t f0, const rta_real_t q,
                                const rta_real_t gain);

/** 
 * Biquad coefficients for an high-shelf filter.
 * H(s) = g * (g*s^2 + (sqrt(g)/q)*s + 1)/(s^2 + (sqrt(g)/q)*s + g)
 * g = sqrt('gain'),
 * 'gain' is linear.
 *
 * @param b is a vector of feed-forward coefficients
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_highshelf_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_real_t f0, const rta_real_t q,
  const rta_real_t gain);

/** 
 * Helper function calling the proper biquad coefficients calculation
 * function, depending on the filter type.
 *
 * @param b is a vector of feed-forward coefficients
 * @param a is a vector of feed-backward coefficients
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_coefs(rta_real_t * b, rta_real_t * a,
                      const rta_filter_t type,
                      const rta_real_t f0, const rta_real_t q,
                      const rta_real_t gain);


/** 
 * Helper function calling the proper biquad coefficients calculation
 * function, depending on the filter type.
 *
 * @param b is a vector of feed-forward coefficients
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients
 * @param a_stride is 'a' stride
 * @param f0 is the cutoff frequency, normalised by the nyquist
 * frequency.
 * @param q must be > 0. and is generally >= 0.5 for audio
 * filtering. q <= 1./sqrt(2.) is the limit for monotonic response.
 * @param gain is linear and must be > 0.
 */
void rta_biquad_coefs_stride(
  rta_real_t * b, const int b_stride,
  rta_real_t * a, const int a_stride,
  const rta_filter_t type,
  const rta_real_t f0, const rta_real_t q,
  const rta_real_t gain);


/** 
 * Biquad computation, using a transposed direct form II.
 *
 * <pre>
 * x
 * --------+---------------+---------------+
 *         |               |               |
 *         |b2             |b1             |b0
 *         V               V               V            y
 *        (+)--->[z-1]--->(+)--->[z-1]--->(+)----+------->
 *         ^               ^                     |
 *         |-a2            |-a1                  |
 *         |               |                     |
 *         +---------------+---------------------+
 * 
 * </pre>
 *
 * @param x is an input sample
 * @param b is a vector of feed-forward coefficients. b0 is b[0], b1
 * is b[1] and b2 is b[2].
 * @param a is a vector of feed-backward coefficients. Note that a1 is
 * a[0] and a2 is a[1] (and a0 is supposed to be 1.).
 * @param state_1 is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 * @param state_2 is the two samples delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 *
 * @return the output sample y
 */
inline rta_real_t rta_biquad(const rta_real_t x,
                             const rta_real_t * b, const rta_real_t * a,
                             rta_real_t * state_1, rta_real_t * state_2);

/** 
 * Biquad computation, using a transposed direct form II.
 *
 * <pre>
 * 
 * x
 * --------+---------------+---------------+
 *         |               |               |
 *         |b2             |b1             |b0
 *         V               V               V            y
 *        (+)--->[z-1]--->(+)--->[z-1]--->(+)----+------->
 *         ^               ^                     |
 *         |-a2            |-a1                  |
 *         |               |                     |
 *         +---------------+---------------------+
 * 
 * </pre> 
 *
 * @param x is an input sample
 * @param b is a vector of feed-forward coefficients. b0 is b[0], b1
 * is b[1] and b2 is b[2].
 * @param b_stride is 'b' stride
 * @param a is a vector of feed-backward coefficients. Note that a1 is
 * a[0] and a2 is a[1] (and a0 is supposed to be 1.).
 * @param a_stride is 'a' stride
 * @param state_1 is the one sample delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 * @param state_2 is the two samples delay state. It can be initialised
 * with 0. or the last computed value, which is updated by this
 * function.
 *
 * @return the output sample y
 */
inline rta_real_t rta_biquad_stride(
  const rta_real_t x,
  const rta_real_t * b, const int b_stride,
  const rta_real_t * a, const int a_stride,
  rta_real_t * state_1, rta_real_t * state_2);



#ifdef __cplusplus
}
#endif

#endif /* _RTA_BIQUAD_H_ */
