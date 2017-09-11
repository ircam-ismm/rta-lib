/**
 * @file   rta_filter.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Fri Aug 29 12:38:46 2008
 * 
 * @brief  Filter types
 * 
 * Copyright (C) 2007 - 2009 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#ifndef _RTA_FILTER_H_
#define _RTA_FILTER_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif


typedef enum
{
  rta_lowpass = 0,
  rta_highpass = 1,
  rta_bandpass_constant_skirt = 2, /* msp name: resonant */
  rta_bandpass_constant_peak = 3, /* msp name: bandpass */
  rta_notch = 4, /* msp name: bandstop */
  rta_allpass = 5,
  rta_peaking = 6, /* msp name: peaknotch */
  rta_lowshelf = 7,
  rta_highshelf = 8
} rta_filter_t;


#ifdef __cplusplus
}
#endif

#endif /* _RTA_FILTER_H_ */
