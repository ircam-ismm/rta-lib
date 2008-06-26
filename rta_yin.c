/**
 *
 * @file   rta_yin.c
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Mon Nov 12 18:21:06 2007
 * 
 * @brief  Yin algorithm for periodicity analysis
 * 
 * Based on Norbert Schnell's real-time implementation after the YIN
 * algorithm by Alain de Cheveigne' et Hideki Kawahara.
 * See "YIN, a fundamental frequency estimator for speech and music"
 * published in J. Acoust. Soc. Am. 111 (4), April 2002.
 *
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#include "rta_yin.h"

#include "rta_stdlib.h" /* rta_malloc, rta_free */
#include "rta_correlation.h" /* rta_correlation_fast */

/* private structure for yin minima search */
typedef struct rta_yin_mins rta_yin_mins_t;
struct rta_yin_mins
{
  rta_real_t min;
  rta_real_t lag;
} ;
  
/* public typedef, private structure */
struct rta_yin_setup
{
  unsigned int max_mins;
  rta_yin_mins_t * mins;
};

int rta_yin_setup_new(rta_yin_setup_t ** yin_setup, unsigned int max_mins)
{
  int ret = 0;
  *yin_setup = (rta_yin_setup_t *) rta_malloc(sizeof(rta_yin_setup_t));

  if(*yin_setup != NULL)
  {
    (*yin_setup)->mins = (rta_yin_mins_t *) rta_malloc(
      sizeof(rta_yin_mins_t) * max_mins);

    if((*yin_setup)->mins != NULL)
    {
      (*yin_setup)->max_mins = max_mins;
      ret = 1;
    }
    else
    {
      rta_free(*yin_setup);
      *yin_setup = NULL;
    }
  }

  return ret;
}

void rta_yin_setup_delete(rta_yin_setup_t * yin_setup)
{
  if(yin_setup != NULL)
  {
    if(yin_setup->mins != NULL)
    {
      rta_free(yin_setup->mins);
    }
    rta_free(yin_setup);
  }
  return;
}

/* Contract: input_size >= 2 * max_lag */
/*           threshold in [0., 1.] == (1. - confidence)^2 */
rta_real_t rta_yin(rta_real_t * abs_min, rta_real_t * autocorrelation, 
                   const rta_real_t * input, const unsigned int input_size,
                   const rta_yin_setup_t * yin_setup,
                   const unsigned int max_lag, const rta_real_t threshold)
{
  rta_real_t biased_threshold = threshold;
  rta_real_t abs_lag = (rta_real_t) max_lag;
  const unsigned int max_mins = yin_setup->max_mins;
  rta_yin_mins_t * mins = yin_setup->mins;

  unsigned int n_mins = 0;    /* number of minimums */
  rta_real_t x;               /* current input sample */
  rta_real_t xm;              /* (current + max_lag) sample */
  rta_real_t energy;
  rta_real_t diff_left, diff, diff_right, sum;
  unsigned int i;             /* input sample index */

  *abs_min = 1.;
    
  /* auto-correlation */
  rta_correlation_fast(autocorrelation, max_lag, input, input, 
                       input_size - max_lag);
    
  /* diff[0] */
  x = input[0];
  xm = input[max_lag];
  energy = autocorrelation[0] + xm * xm - x * x;
  diff_left = 0.0;
  diff = 0.0;
  diff_right = autocorrelation[0] + energy - 2. * autocorrelation[1];
  sum = 0.0;

  /* diff[1] */
  x = input[1];
  xm = input[1 + max_lag];
  energy += xm * xm - x * x;
  diff_left = diff;
  diff = diff_right;
  diff_right = autocorrelation[0] + energy - 2. * autocorrelation[2];
  sum = diff;

  /* minimum difference search */
  for(i=2; i<max_lag-1 && n_mins < max_mins; i++)
  {
    x = input[i];
    xm = input[i + max_lag];
    energy += xm * xm - x * x;
    diff_left = diff;
    diff = diff_right;
    diff_right = autocorrelation[0] + energy - 2. * autocorrelation[i + 1];
    sum += diff;
    
    /* local minimum */
    if(diff < diff_left && diff < diff_right && sum != 0.)
    {
      const rta_real_t a = diff_left + diff_right - 2.0 * diff;
      const rta_real_t b = 0.5 * (diff_right - diff_left);
      const rta_real_t min = diff - (b * b) / (2.0 * a);
      const rta_real_t min_norm = i * min / sum;
      const rta_real_t lag = i - b / a;
      
      mins[n_mins].min = min_norm;
      mins[n_mins].lag = lag;
      n_mins++;

      if(min_norm < *abs_min)
      {
        *abs_min = min_norm;
        abs_lag = lag;
      }
    }
  }

  /* bias threshold by absolute minimum */
  biased_threshold += *abs_min;

  for(i=0; i<n_mins; i++)
  {
    /* look for first minimum under biased threshold */
    if(mins[i].min < biased_threshold)
    {
      *abs_min = mins[i].min;
      abs_lag = mins[i].lag;
      break;
    }
  }

  /* clip in the [0., 1.] range (rounding errors) */
  if(*abs_min < 0.) 
  {
    *abs_min = 0.;
  }

  return abs_lag;
}

/* Contract: input_size >= 2 * max_lag */
/*           threshold in [0., 1.] == (1. - confidence)^2 */
rta_real_t rta_yin_stride(
  rta_real_t * abs_min,
  rta_real_t * autocorrelation, const unsigned int ac_stride,
  const rta_real_t * input, const unsigned int i_stride,
  const unsigned int input_size,
  const rta_yin_setup_t * yin_setup,
  const unsigned int max_lag, const rta_real_t threshold)
{
  rta_real_t biased_threshold = threshold;
  rta_real_t abs_lag = (rta_real_t) max_lag;
  const unsigned int max_mins = yin_setup->max_mins;
  rta_yin_mins_t * mins = yin_setup->mins;

  unsigned int n_mins = 0;    /* number of minimums */
  rta_real_t x;               /* current input sample */
  rta_real_t xm;              /* (current + max_lag) sample */
  rta_real_t energy;
  rta_real_t diff_left, diff, diff_right, sum;
  unsigned int i,is;          /* input sample index, and with stride */
  unsigned int ac;            /* autocorrelation index */
  const unsigned int max_lag_stride = max_lag * i_stride;

  *abs_min = 1.;
    
  /* auto-correlation */
  rta_correlation_fast_stride(autocorrelation, ac_stride, max_lag, 
                              input, i_stride, input, i_stride,
                              input_size - max_lag);
    
  /* diff[0] */
  x = input[0];
  xm = input[max_lag_stride];
  energy = autocorrelation[0] + xm * xm - x * x;
  diff_left = 0.0;
  diff = 0.0;
  diff_right = autocorrelation[0] + energy - 2. * autocorrelation[ac_stride];
  sum = 0.0;

  /* diff[1] */
  x = input[i_stride];
  xm = input[i_stride + max_lag_stride];
  energy += xm * xm - x * x;
  diff_left = diff;
  diff = diff_right;
  diff_right = autocorrelation[0] + energy - 2. * autocorrelation[2*ac_stride];
  sum = diff;

  /* minimum difference search */
  for(i = 2, is = 2*i_stride, ac = 3*ac_stride;
      i < max_lag - 1 && n_mins < max_mins;
      i++, is += i_stride, ac += ac_stride)
  {
    x = input[is];
    xm = input[is + max_lag_stride];
    energy += xm * xm - x * x;
    diff_left = diff;
    diff = diff_right;
    diff_right = autocorrelation[0] + energy - 2. * autocorrelation[ac];
    sum += diff;
    
    /* local minimum */
    if(diff < diff_left && diff < diff_right && sum != 0.)
    {
      const rta_real_t a = diff_left + diff_right - 2.0 * diff;
      const rta_real_t b = 0.5 * (diff_right - diff_left);
      const rta_real_t min = diff - (b * b) / (2.0 * a);
      const rta_real_t min_norm = i * min / sum;
      const rta_real_t lag = i - b / a;
      
      mins[n_mins].min = min_norm;
      mins[n_mins].lag = lag;
      n_mins++;

      if(min_norm < *abs_min)
      {
        *abs_min = min_norm;
        abs_lag = lag;
      }
    }
  }

  /* bias threshold by absolute minimum */
  biased_threshold += *abs_min;

  for(i=0; i<n_mins; i++)
  {
    /* look for first minimum under biased threshold */
    if(mins[i].min < biased_threshold)
    {
      *abs_min = mins[i].min;
      abs_lag = mins[i].lag;
      break;
    }
  }

  /* clip in the [0., 1.] range (rounding errors) */
  if(*abs_min < 0.) 
  {
    *abs_min = 0.;
  }

  return abs_lag;
}
