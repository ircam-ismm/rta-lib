/**
 * @file   rta_mean_variance.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Aug 25 16:13:42 2008
 * 
 * @brief Mean and variance from an input vector
 * 
 * Copyright (C) 2008 by IRCAM-Centre Georges Pompidou, Paris, France.
 */

#ifndef _RTA_MEAN_VARIANCE_H_
#define _RTA_MEAN_VARIANCE_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * Variance is computed as var(x) = E(x^2) - mean(x)^2 so as it is
 * done in one loop. Note that this may lead to inacuracies when
 * E(x^2) and mean(x)^2 are similar in magnitude. The mean and
 * variance are normalised by 'i_size', hence the bias.
 * \see rta_variance
 * \see rta_mean_variance_unbiased
 * 
 * @param mean is a pointer to the mean result
 * @param variance is a pointer to the variance result
 * @param input is the input vector of size 'i_size'
 * @param i_size is the size of 'input' and must be > 0
 */
void
rta_mean_variance(rta_real_t * mean, rta_real_t * variance,
                  rta_real_t * input, const unsigned int i_size);

/** 
 * Variance is computed as var(x) = E(x^2) - mean(x)^2 so as it is
 * done in one loop. Note that this may lead to inacuracies when
 * E(x^2) and mean(x)^2 are similar in magnitude. The mean and
 * variance are normalised by 'i_size', hence the bias.
 * \see rta_variance
 * \see rta_mean_variance_unbiased
 * 
 * @param mean is a pointer to the mean result
 * @param variance is a pointer to the variance result
 * @param input is the input vector of size 'i_size'
 * @param i_stride is the 'input' stride
 * @param i_size is the size of 'input' and must be > 0
 */
void
rta_mean_variance_stride(
  rta_real_t * mean, rta_real_t * variance,
  rta_real_t * input, const int i_stride, const unsigned int i_size);

/** 
 * Variance is computed as var(x) = E(x^2) - mean(x)^2 so as it is
 * done in one loop. Note that this may lead to inacuracies when
 * E(x^2) and mean(x)^2 are similar in magnitude. The mean and
 * variance are normalised by ('i_size' - 1).
 * \see rta_variance
 * \see rta_mean_variance
 * 
 * @param mean is a pointer to the mean result
 * @param variance is a pointer to the variance result
 * @param input is the input vector of size 'i_size'
 * @param i_size is the size of 'input' and must be > 0
 */
void
rta_mean_variance_unbiased(rta_real_t * mean, rta_real_t * variance,
                           rta_real_t * input, const unsigned int i_size);

/** 
 * Variance is computed as var(x) = E(x^2) - mean(x)^2 so as it is
 * done in one loop. Note that this may lead to inacuracies when
 * E(x^2) and mean(x)^2 are similar in magnitude. The mean and
 * variance are normalised by ('i_size' - 1).
 * \see rta_variance
 * \see rta_mean_variance
 * 
 * @param mean is a pointer to the mean result
 * @param variance is a pointer to the variance result
 * @param input is the input vector of size 'i_size'
 * @param i_stride is the 'input' stride
 * @param i_size is the size of 'input' and must be > 0
 */
void
rta_mean_variance_unbiased_stride(
  rta_real_t * mean, rta_real_t * variance,
  rta_real_t * input, const int i_stride, const unsigned int i_size);

/** 
 * Mean of 'input'
 * 
 * @param input is the input vector of size 'i_size'
 * @param i_size is the size of 'input' and must be > 0
 * 
 * @return mean of 'input'
 */
rta_real_t
rta_mean(rta_real_t * input, const unsigned int i_size);

/** 
 * Mean of 'input'
 * 
 * @param input is the input vector of size 'i_size'
 * @param i_stride is the 'input' stride
 * @param i_size is the size of 'input' and must be > 0
 * 
 * @return mean of 'input'
 */
rta_real_t
rta_mean_stride(
  rta_real_t * input, const int i_stride, const unsigned int i_size);

/** 
 * Variance is computed as var(x) = E( (x - mean(x))^2 ) and is
 * normalised by 'i_size', hence the bias.
 * \see rta_variance_unbiased
 * \see rta_mean_variance
 * \see rta_mean
 * 
 * @param input is the input vector of size 'i_size'
 * @param i_size is the size of 'input' and must be > 0
 * @param mean is the mean of 'input' 
 * 
 * @return variance of 'input'
 */
rta_real_t
rta_variance(rta_real_t * input, const unsigned int i_size, rta_real_t mean);

/** 
 * Variance is computed as var(x) = E( (x - mean(x))^2 ) and is
 * normalised by 'i_size', hence the bias.
 * \see rta_variance_unbiased
 * \see rta_mean_variance
 * \see rta_mean
 * 
 * @param input is the input vector of size 'i_size'
 * @param i_stride is the 'input' stride
 * @param i_size is the size of 'input' and must be > 0
 * @param mean is the mean of 'input' 
 * 
 * @return variance of 'input'
 */
rta_real_t
rta_variance_stride(
  rta_real_t * input, const int i_stride, const unsigned int i_size, 
  rta_real_t mean);

/** 
 * Variance is computed as var(x) = E( (x - mean(x))^2 ) and is
 * normalised by ('i_size' - 1).
 * \see rta_variance
 * \see rta_mean_variance
 * \see rta_mean
 * 
 * @param input is the input vector of size 'i_size'
 * @param i_size is the size of 'input' and must be > 0
 * @param mean is the mean of 'input' 
 * 
 * @return variance of 'input'
 */
rta_real_t
rta_variance_unbiased(rta_real_t * input, const unsigned int i_size, 
                      rta_real_t mean);

/** 
 * Variance is computed as var(x) = E( (x - mean(x))^2 ) and is
 * normalised by ('i_size' - 1).
 * \see rta_variance
 * \see rta_mean_variance
 * \see rta_mean
 * 
 * @param input is the input vector of size 'i_size'
 * @param i_stride is the 'input' stride
 * @param i_size is the size of 'input' and must be > 0
 * @param mean is the mean of 'input' 
 * 
 * @return variance of 'input'
 */
rta_real_t rta_variance_unbiased_stride(
  rta_real_t * input, const int i_stride, const unsigned int i_size, 
  rta_real_t mean);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_MEAN_VARIANCE_H_ */

