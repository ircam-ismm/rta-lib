/**
 *
 *
 * @file   rta_mel.h
 * @author Jean-Philippe.Lambert@ircam.fr
 * @date   Fri Jun 15 15:29:25 2007
 * 
 * @brief  Mel conversions (HTK and Auditory Toolbox styles)
 * 
 * Based on Rastamat by Dan Ellis.
 * See http://www.ee.columbia.edu/~dpwe/resources/matlab/rastamat
 *
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#ifndef _RTA_MEL_H_
#define _RTA_MEL_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif

   
typedef enum
{
    rta_mel_slaney = 1,   /**< Slaney-style mel is scaled to be approx
                       * constant E per channel */
    rta_mel_htk = 2       /**< HTK-style is constant max amplitude per
                       * channel */
} rta_mel_t;


/** MEL from melfcc.m */

/** from hz2mel.m */
/**
 * function pointer to avoid tests during conversions
 * rta_real_t parameter in frequency in Hertz
 * rta_real_t return is corresponding mel value
 */
typedef rta_real_t (*rta_hz_to_mel_function) (rta_real_t);

/** 
 * Convert frequencies f (in Hz) to mel 'scale'.
 * Mel fn to match Slaney's Auditory Toolbox mfcc.m
 * 
 * @param freq_in_hz [0.,22050.]
 * 
 * @return corresponding mel value [0.,60.]
 */
rta_real_t rta_hz_to_mel_slaney(rta_real_t freq_in_hz);

/** 
 * Convert frequencies f (in Hz) to mel 'scale'.
 * uses the mel axis defined in the htk_book
 * 
 * @param freq_in_hz [0.,22050.]
 * 
 * @return corresponding mel value [0.,3923.]
 */
rta_real_t rta_hz_to_mel_htk(rta_real_t freq_in_hz);

/** from mel2hz.m */
/**
 * function pointer to avoid tests during conversions
 * rta_real_t parameter is mel value
 * rta_real_t return is corresponding frequency
 */
typedef rta_real_t (*rta_mel_to_hz_function) (rta_real_t);

/** 
 * Convert 'mel scale' frequencies into Hz
 * use the formula from Slaney's mfcc.m
 * 
 * @param freq_in_mel [0.,60.]
 * 
 * @return corresponding frequency [0.,22050.]
 */
rta_real_t rta_mel_to_hz_slaney(rta_real_t freq_in_mel);

/** 
 * Convert 'mel scale' frequencies into Hz
 * use the HTK formula
 * 
 * @param freq_in_mel [0.,3923.]
 * 
 * @return corresponding frequency [0.,22050.]
 */
rta_real_t rta_mel_to_hz_htk(rta_real_t freq_in_mel);

#ifdef __cplusplus
}
#endif

#endif /* _RTA_MEL_H_ */

