/**
 * @file   rta_moments.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Thu Dec 13 15:28:26 2007
 * 
 * @brief  Statistical moments functions
 * 
 * The moments are calculated over the indexes and weighted by the
 * input values (eg. the amplitudes of the spectrum regularly
 * sampled). Note that all moments (but the first) are centered. The
 * results unit is index (starting from 0).
 *
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */

#ifndef _RTA_MOMENTS_H_
#define _RTA_MOMENTS_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

/** 
 * First order moment over indexes weighted by 'input' values:
 * weighted mean, centroid.
 * m1 = centroid = sum(i, i*input(i)) / sum(i, input(i))
 * 
 * @param input_sum is calculated by the function and may be used for
 * higher-order moments. It is the sum of all 'input' values.
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * 
 * @return the first moment, 0.5 * ('input_size' - 1) if 'input_sum' == 0.
 */
rta_real_t
rta_weighted_moment_1_indexes(
  rta_real_t * input_sum,
  const rta_real_t * input, const unsigned int input_size);

/** 
 * First order moment over indexes weighted by 'input' values:
 * weighted mean, centroid.
 * m1 = centroid = sum(i, i*input(i)) / sum(i, input(i))
 * 
 * @param input_sum is calculated by the function and may be used for
 * higher-order moments. It is the sum of all 'input' values.
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * 
 * @return the first moment, 0.5 * ('input_size' - 1) if 'input_sum' == 0.
 */
rta_real_t
rta_weighted_moment_1_indexes_stride(
  rta_real_t * input_sum,
  const rta_real_t * input, const int i_stride,
  const unsigned int input_size);

/** 
 * Second order weighted central moment over indexes: spread, weighted
 * variance.
 * m2 = spread = sum(i,input(i) * (i-centroid)^2) / sum(i, input(i))
 * 
 * Note that standard deviation is std = sqrt(spread)
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * 
 * @return spread
 */
rta_real_t
rta_weighted_moment_2_indexes(
  const rta_real_t * input, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum);

/** 
 * Second order weighted central moment over indexes: spread, weighted
 * variance.
 * m2 = spread = sum(i,input(i) * (i-centroid)^2) / sum(i, input(i))
 * 
 * Note that standard deviation is std = sqrt(spread)
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes_stride
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * 
 * @return spread
 */
rta_real_t
rta_weighted_moment_2_indexes_stride(
  const rta_real_t * input, const int i_stride, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum);

/** 
 * Third order weighted central moment over indexes.
 * m3 = sum(i,input(i) * (i-centroid)^3) / sum(i, input(i))
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'. 
 * \see rta_weighted_moment_1_indexes
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * 
 * @return third central moment
 */
rta_real_t
rta_weighted_moment_3_indexes(
  const rta_real_t * input, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum);

/** 
 * Third order weighted central moment over indexes.
 * m3 = sum(i,input(i) * (i-centroid)^3) / sum(i, input(i))
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'. 
 * \see rta_weighted_moment_1_indexes_stride
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * 
 * @return third central moment
 */
rta_real_t
rta_weighted_moment_3_indexes_stride(
  const rta_real_t * input, const int i_stride, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum);

/** 
 * Third order standardised weighted central moment over indexes: skewness.
 * skewness = m3 / std^3
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param deviation must be != 0. It is the standard deviation.
 * \see rta_weighted_moment_2_indexes
 * 
 * @return skewness
 */
rta_real_t 
rta_std_weighted_moment_3_indexes(
  const rta_real_t * input, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t deviation);

/** 
 * Third order standardised weighted central moment over indexes: skewness.
 * skewness = m3 / std^3
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes_stride 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param deviation must be != 0. It is the standard deviation.
 * \see rta_weighted_moment_2_indexes_stride
 * 
 * @return skewness
 */
rta_real_t 
rta_std_weighted_moment_3_indexes_stride(
  const rta_real_t * input, const int i_stride, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t deviation);

/** 
 * Fourth order weighted central moment over indexes.
 * m4 = sum(i,input(i) * (i-centroid)^4) / sum(i, input(i))
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'. 
 * \see rta_weighted_moment_1_indexes 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * 
 * @return fourth central moment
 */
rta_real_t
rta_weighted_moment_4_indexes(
  const rta_real_t * input, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum);

/** 
 * Fourth order weighted central moment over indexes.
 * m4 = sum(i,input(i) * (i-centroid)^4) / sum(i, input(i))
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'. 
 * \see rta_weighted_moment_1_indexes_stride 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * 
 * @return fourth central moment
 */
rta_real_t
rta_weighted_moment_4_indexes_stride(
  const rta_real_t * input, const int i_stride, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum);

/** 
 * Fourth order standardised weighted central moment over indexes: kurtosis.
 * kurtosis = m4 / std^4
 * 
 * Note that the kurtosis is often defined as the fourth cumulant
 * divided by the square root of the variance, which gives
 * kurtosis = m4 / std^4 - 3. This function does not include the "- 3"
 * term.
 *
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'. 
 * \see rta_weighted_moment_1_indexes 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param deviation must be != 0. It is the standard deviation.
 * \see rta_weighted_moment_2_indexes
 * 
 * @return kurtosis
 */
rta_real_t
rta_std_weighted_moment_4_indexes(
  const rta_real_t * input, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t deviation);

/** 
 * Fourth order standardised weighted central moment over indexes: kurtosis.
 * kurtosis = m4 / std^4
 * 
 * Note that the kurtosis is often defined as the fourth cumulant
 * divided by the square root of the variance, which gives
 * kurtosis = m4 / std^4 - 3. This function does not include the "- 3"
 * term.
 *
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'. 
 * \see rta_weighted_moment_1_indexes_stride 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param deviation must be != 0. It is the standard deviation.
 * \see rta_weighted_moment_2_indexes_stride
 * 
 * @return kurtosis
 */
rta_real_t
rta_std_weighted_moment_4_indexes_stride(
  const rta_real_t * input, const int i_stride, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t deviation);

/** 
 * General order weighted central moment over indexes.
 * m_order = sum(i,input(i) * (i-centroid)^order) / sum(i, input(i))
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param order is the moment order.
 * 
 * @return moment
 */
rta_real_t
rta_weighted_moment_indexes(
  const rta_real_t * input, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t order);

/** 
 * General order weighted central moment over indexes.
 * m_order = sum(i,input(i) * (i-centroid)^order) / sum(i, input(i))
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes_stride 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param order is the moment order.
 * 
 * @return moment
 */
rta_real_t
rta_weighted_moment_indexes_stride(
  const rta_real_t * input, const int i_stride, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t order);

/** 
 * General order standardised weighted central moment over indexes.
 * m_order / std^order
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param deviation must be != 0. It is the standard deviation.
 * \see rta_weighted_moment_2_indexes
 * @param order is the moment order.
 * 
 * @return standardised moment
 */
rta_real_t
rta_std_weighted_moment_indexes(
  const rta_real_t * input, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t deviation,
  const rta_real_t order);

/** 
 * General order standardised weighted central moment over indexes.
 * m_order / std^order
 * 
 * @param input is usually amplitudes or weights. Each element of
 * 'input' must be >=0. 
 * @param i_stride is 'input' stride
 * @param input_size is 'input' size
 * @param centroid is the first moment of 'input'.
 * \see rta_weighted_moment_1_indexes 
 * @param input_sum must be != 0. It is the sum of all 'input' values.
 * @param deviation must be != 0. It is the standard deviation.
 * \see rta_weighted_moment_2_indexes
 * @param order is the moment order.
 * 
 * @return standardised moment
 */
rta_real_t
rta_std_weighted_moment_indexes_stride(
  const rta_real_t * input, const int i_stride, const unsigned int input_size,
  const rta_real_t centroid, const rta_real_t input_sum,
  const rta_real_t deviation,
  const rta_real_t order);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_MOMENTS_H_ */
