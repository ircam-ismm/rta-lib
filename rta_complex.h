/**
 * @file   rta_complex.h
 * @author Jean-Philippe Lambert
 * @date   Mon Sep 10 11:05:09 2007
 * 
 * @brief  rta_complex_t type warper for float, double or long double
 * complex.
 *
 * Default is the same as RTA_REAL_TYPE. Define your RTA_COMPLEX_TYPE to
 * override these.
 * \see rta_configuration.h
 * \see rta_real.h
 * 
 * Copyright (C) 2007 by IRCAM-Centre Georges Pompidou, Paris, France.
 * 
 */
#ifndef _RTA_COMPLEX_H_
#define _RTA_COMPLEX_H_ 1

#include "rta.h"

/** default complex precision is the same as real precision */
#ifndef RTA_COMPLEX_TYPE
#define RTA_COMPLEX_TYPE RTA_REAL_TYPE
#endif

#if (RTA_COMPLEX_TYPE == RTA_FLOAT_TYPE)
#undef rta_complex_t
#define rta_complex_t float complex
#endif

#if (RTA_COMPLEX_TYPE == RTA_DOUBLE_TYPE)
#undef rta_complex_t 
#define rta_complex_t double complex
#endif

#if (RTA_COMPLEX_TYPE == RTA_LONG_DOUBLE_TYPE)
#undef rta_complex_t
#define rta_complex_t long double complex
#endif


#ifdef WIN32

#define rta_cabs cabs
#define rta_cacos cacos
#define rta_cacosh cacosh
#define rta_carg carg
#define rta_casin casin
#define rta_casinh casinh
#define rta_catan catan
#define rta_catanh catanh
#define rta_ccos ccos
#define rta_ccosh ccosh
#define rta_cexp cexp
#define rta_cimag cimag
#define rta_clog clog
#define rta_conj conj
#define rta_cpow cpow
#define rta_cproj cproj
#define rta_creal creal
#define rta_csin csin
#define rta_csinh csinh
#define rta_csqrt csqrt
#define rta_ctan ctan
#define rta_ctanh ctanh

#else

#if (RTA_COMPLEX_TYPE == RTA_FLOAT_TYPE)

#define rta_cabs cabsf
#define rta_cacos cacosf
#define rta_cacosh cacoshf
#define rta_carg cargf
#define rta_casin casinf
#define rta_casinh casinhf
#define rta_catan catanf
#define rta_catanh catanhf
#define rta_ccos ccosf
#define rta_ccosh ccoshf
#define rta_cexp cexpf
#define rta_cimag cimagf
#define rta_clog clogf
#define rta_conj conjf
#define rta_cpow cpowf
#define rta_cproj cprojf
#define rta_creal crealf
#define rta_csin csinf
#define rta_csinh csinhf
#define rta_csqrt csqrtf
#define rta_ctan ctanf
#define rta_ctanh ctanhf

#endif

#if (RTA_COMPLEX_TYPE == RTA_DOUBLE_TYPE)

#define rta_cabs cabs
#define rta_cacos cacos
#define rta_cacosh cacosh
#define rta_carg carg
#define rta_casin casin
#define rta_casinh casinh
#define rta_catan catan
#define rta_catanh catanh
#define rta_ccos ccos
#define rta_ccosh ccosh
#define rta_cexp cexp
#define rta_cimag cimag
#define rta_clog clog
#define rta_conj conj
#define rta_cpow cpow
#define rta_cproj cproj
#define rta_creal creal
#define rta_csin csin
#define rta_csinh csinh
#define rta_csqrt csqrt
#define rta_ctan ctan
#define rta_ctanh ctanh

#endif


#if (RTA_COMPLEX_TYPE == RTA_LONG_DOUBLE_TYPE)

#define rta_cabs cabsl
#define rta_cacos cacosl
#define rta_cacosh cacoshl
#define rta_carg cargl
#define rta_casin casinl
#define rta_casinh casinhl
#define rta_catan catanl
#define rta_catanh catanhl
#define rta_ccos ccosl
#define rta_ccosh ccoshl
#define rta_cexp cexpl
#define rta_cimag cimagl
#define rta_clog clogl
#define rta_conj conjl
#define rta_cpow cpowl
#define rta_cproj cprojl
#define rta_creal creall
#define rta_csin csinl
#define rta_csinh csinhl
#define rta_csqrt csqrtl
#define rta_ctan ctanl
#define rta_ctanh ctanhl

#endif

#endif /* WIN32 */

#include <complex.h>

#endif /* _RTA_COMPLEX_H_ */
