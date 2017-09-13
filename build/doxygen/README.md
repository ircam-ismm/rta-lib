# rta-lib overview {#mainpage}

# Table of Contents

- 1 [Generalities](#generalities)
- 2 [Library configuration](#configuration)
- 3 [Descriptors based on a vector of samples](#descr_vector_samples)
    - 3.1 [Re-sampling](#resampling)
    - 3.2 [Yin](#yin)
    - 3.3 [Gain and integer-float conversion](#gain_integer_float)
    - 3.4 [Pre-emphasis](#pre_emphasis)
    - 3.5 [Windowing](#windowing)
    - 3.6 [Linear predictive coefficients (LPC)](#lpc)
- 4 [Descriptors based on the power spectrum](#descr_power_spectrum)
    - 4.1 [Real Fourier transform](#real_fourier_transform)
    - 4.2 [Complex spectrum to power spectrum](#complex_to_power_spectrum)
    - 4.3 [Statistical moments](#statistical_moments)
        - 4.3.1 [Centroid](#centroid)
        - 4.3.2 [Spread and deviation](#spread_deviation)
        - 4.3.3 [Skewness](#skewness)
        - 4.3.4 [Kurtosis](#kurtosis)
- 5 [Descriptors based on the mel bands](#descr_mel_bands)
    - 5.1 [Mel bands](#mel_bands)
    - 5.2 [Mel-frequency cepstral coefficients (MFCC)](#mfcc)
    - 5.3 [Delta and delta-delta MFCC](#delta_delta_mfcc)
- 6 [Style guide](#style_guide)
    - 6.1 [Prefix](#prefix)
    - 6.2 [Declarations and definitions](#declarations_definitions)
    - 6.3 [Headers](#headers)
    - 6.4 [Real types](#real_types)
    - 6.5 [Memory allocation](#memory_allocation)

<!-- section 1 -->
<a name="generalities"></a>
# 1 Generalities

The _RTA_ library is frame-based, which means that for any vector of
samples, a set of descriptors can be computed without adding any
delay. A noticeable exception to this is the _delta_ (and
_delta-delta_) computation, as it is based on more than one frame.

There are two variants of each function, one being post-fixed with `_stride`.
It allows to directly access interleaved data without copying them (like a
stereophonic samples vector).
However, using a big _stride_ value may lead to a bad usage of the memory cache.

Any index starts at 0.

Some functions require an initialisation before any processing,
in order to pre-calculate what will not depend on the incoming frame.
Every allocation must be done before anything else, outside of the functions
themselves except mentioned otherwise.

Some descriptors can be computed by several functions, and the results
may slightly differ for several reasons: the functions do not rely
on the same algorithms and the signals used for the computation may
differ (due to windowing, filtering, etc.). The auto-correlation from
_yin_ and from the _LPC_ are not the same and there is a
lot of ways to get the energy: from _yin_, as the sum of the
squares of the samples, from _LPC_, or as the first
_MFCC_ coefficient.


<!-- section 2 -->
<a name="configuration"></a>
# 2 Library configuration

A file named `rta_configuration.h` must be present within your sources
in order to use the _RTA_ library.
It is not included within the _RTA_ source files.
An empty file means that the defaults settings are used for the compilation.

Instead of using the `malloc`, `realloc` and `free` functions from the `stdlib.h`,
one can respectively define `rta_malloc`, `rta_realloc` and `rta_free`.
(see [6.5](#memory_allocation))

The floating-point precision can be simple, double or long double, according to
the definition of `RTA_REAL_TYPE` to respectively `RTA_FLOAT_TYPE`,
`RTA_DOUBLE_TYPE` or `RTA_LONG_DOUBLE_TYPE`.
The constants in `rta_float.h` are then redefined according to the
proper type from `float.h`.
The same applies to the functions in `rta_math.h` from `math.h`.
(see [6.4](#real_types))

Note that the long double precision is not supported when using Apple's VecLib
by setting `RTA_USE_VECLIB` to 1.

<!-- INCLUDE FIGURE DATAFLOW -->
<div align="center" style="margin:30px;">
<table style="border:0; background-color:rgba(255,255,255,1);">
  <tr><td>
    <img src="rta_dataflow.jpg" width="600px" style="max-width:100%;"/>
  </td></tr>
  <tr><td align="center">
    _Sound descriptors data-flow_
  </td></tr>
</table>
</div>

<!-- section 3 -->
<a name="descr_vector_samples"></a>
# 3 Descriptors based on a vector of samples

The vector of samples is characterised by its sample-rate and its
size. Moreover, to use sizes larger than those provided by the
sound card (e.g. for _Fourier_ transforms), the hop-size gives the number
of samples between two consecutive vectors that overlap if the hop-size
is smaller than the vector-size.

<a name="resampling"></a>
## 3.1 Re-sampling

Down-sampling a signal is interesting to lower the further computations,
especially for the computational-intensive algorithms, like _yin_.
It can help to concentrate on a range of the spectrum where the information is
pertinent for the analysis: the _LPC_ (which is linear, as its name suggests it)
finds poles and zeros to fit the whole spectrum; to keep the results under 5 kHz,
we simply use a sample-rate of 11 kHz.
If the original signal sample-rate is 44 kHz, we can use the function
`rta_downsample_int_mean` with a factor of 4. This function implies a low-pass
filtering.
The function `rta_downsample_int_remove` can produce aliasing if the original
signal contains information above half of the resulting sample-rate.
Note that the resulting samples vector-size is smaller, according to the given
factor.

<a name="yin"></a>
## 3.2 Yin

The _yin_ algorithm computes the periodicity of a samples
vector, finding its most probable _lag_. To get the period in Hz,
on simply multiply this lag by the input vector sample-rate.

Note that the _yin_ algorithm operates on a non-windowed
samples vector. As it can be computational-intensive, a down-sampling
of the incoming vector is often performed: to track a pitch below 1
kHz, one can use a sampling-rate of 11 kHz, or even 5.5 kHz, depending
on the results quality request. One can then check the _absolute minimum_
found, which gives the _periodicity_ as long as the absolute minimum is
positive:
\f[periodicity = 1 - \sqrt {absolute\_minimum}\f]

Before computing anything, a new _yin_ structure must be allocated and filled
with the `rta_yin_setup_new` function.
It will be released by the function `rta_yin_setup_delete`.
(The auto-correlation result vector must also be allocated beforehand,
like any results vector.)

<a name="gain_integer_float"></a>
## 3.3 Gain and integer-float conversion

The samples are generally coded by floating-points number over 32 bits, within
the range [-1.0, 1.0]. To convert them into 16 bits signed integers (which is
the format used by some systems, like _HTK_), one can apply a simple gain of
\f$2^{15}\f$ by multiplying every sample.

<a name="pre_emphasis"></a>
## 3.4 Pre-emphasis

The pre-emphasis is a simple first-order difference between the current sample
and the previous one (weighted by a factor): \f$s(n) = s(n) - f \times s(n-1)\f$
<br />
It is often used for voice analysis with a factor of 0.97 as it reduces the low
frequencies while raising the high frequencies, thus amplifying the contrast.

<a name="windowing"></a>
## 3.5 Windowing

If the samples vector-size is known, it is possible to pre-calculate
the weights that will be used to apply a given function.
- the function `rta_window_hamming_weights` computes a _Hamming_ window
- the function `rta_window_hann_weights` computes a _von Hann_ window.

These (or any weights vector) can be applied with the `rta_window_apply`
function. <br />
The functions post-fixed with `_in_place` change the input samples
vector values directly.

If the samples vector-size is not known in advance, one can still
apply the window using the `rta_window_rounded_apply` function.
There is no interpolation, then. The weights vector indexes are simply scaled
and rounded: this is efficient but the rounding error may be unacceptable if the
size of the weights vector is too small comparing with the samples vector size.
It is also possible to compute and apply a window on the fly, with the functions
`rta_window_hann_apply` and `rta_window_hamming_apply`.

<a name="lpc"></a>
## 3.6 Linear predictive coefficients (LPC)

The `rta_lpc` function calculates the linear predictive
coefficients (_LPC_) for a samples vector, using an auto-correlation
and a _Levinson-Durbin_ decomposition. Note that the _LPC_
order is one value less than the _LPC_ size.

The first _LPC_ coefficient is always 1 and is often replaced (e.g. in
_HTK_) by the prediction error, which gives the energy of the samples
vector. If the _LPC_ is computed on overlapping samples vectors, they
are often windowed in order for the coefficients to evolve smoothly
from frame to frame, and the energy is then reduced (by a constant
factor depending on the window).


<!-- section 4  -->
<a name="descr_power_spectrum"></a>
# 4 Descriptors based on the power spectrum

Some descriptors are based on the power spectrum of a samples
vector. The samples vector is first windowed. Then a real
_Fourier_ transform is applied, giving a complex spectrum. The
power spectrum is the square of the magnitude of the complex spectrum.

<a name="real_fourier_transform"></a>
## 4.1 Real Fourier transform

Before computing a real _Fourier_ transform, a new real _Fourier_ transform
setup must be allocated and filled with the `rta_fft_real_setup_new` function,
with the type `real_to_complex_1d`.
It will be released by the function `rta_fft_setup_delete`.
The function `rta_fft_execute` applies the _Fourier_ transform to a samples
vector.

The transform size must be a power of 2. If the transform size is
bigger than the actual samples vector, it is then padded with zeros.

By convention (e.g. _HTK_), no scale is applied to this direct
transform. The inverse of the transform size can later be applied to
the inverse transform in order to obtain the identity transform.

<a name="complex_to_power_spectrum"></a>
## 4.2 Complex spectrum to power spectrum

The power spectrum is the square of the magnitude the complex
spectrum.
<!-- % $$ power\_spectrum = |X|^2$$ -->
Its size is half the size of the _Fourier_ transform plus one
(the last element corresponds to the _Nyquist_ frequency). To
get a correspondence between the power spectrum index and the
corresponding frequency, one can apply a simple ratio between the
maximum frequency (which is half of the sample-rate) and the maximum
index (which is half of the _Fourier_ transform size, as all
the indexes start at 0):
\f[frequency = index \frac
  {sample\_rate}
  {transform\_size}
\f]

<a name="statistical_moments"></a>
## 4.3 Statistical moments

The statistical moments can be computed from any samples vector, as
long as the weights are positive, as they represent the probability of
the random variables to appear. The power spectrum conforms to this, as
any value is positive, but not the amplitude spectrum (if not
translated above 0). The same applies for the moments of the mel
bands.

The moments are calculated over the indexes and weighted by the input
values. They will be normalised by the sum of the input values in order
to get the indexes probability. Note that all moments (but the first)
are centred. Any moment above the second can be standardised.

The moments described hereafter describe the power spectrum moments.

<a name="centroid"></a>
### 4.3.1 Centroid

The spectral centroid is the first moment over the indexes weighted by
the vector of power spectrum values. It is computed by
`rta_weighted_moment_1_indexes`. The result unit is
\f$index\f$ (of the power spectrum, starting at 0). This function
returns also the input sum as it can be used in further calculations.
\f[m_1 = centroid = \frac
  {\sum_i i \times input(i)}
  {\sum_i input(i)}
\f]

<a name="spread_deviation"></a>
### 4.3.2 Spread and deviation

The spectral spread is the second central moment over the indexes
weighted by the vector of power spectrum values. It is computed by
`rta_weighted_moment_2_indexes`. The result unit is
\f$index^2\f$ (of the power spectrum, starting at 0).
\f[m_2 = spread = \frac
  {\sum_i (i - centroid )^2 \times input(i)}
  {\sum_i input(i)}
\f]

The standard deviation is \f$std = \sqrt {spread}\f$.

<a name="skewness"></a>
### 4.3.3 Skewness

The spectral skewness is the third standard central moment over the
indexes weighted by the vector of power spectrum values. It is
computed by `rta_std_weighted_moment_3_indexes`. The
result is without unit.
\f[m_{3std} = skewness = \frac
  {\sum_i (i - centroid )^3 \times input(i)}
  {std^3 \sum_i input(i)}
\f]

<a name="kurtosis"></a>
### 4.3.4 Kurtosis

The spectral kurtosis is the fourth standard central moment over the
indexes weighted by the vector of power spectrum values. It is
computed by `rta_std_weighted_moment_4_indexes`. The
result is without unit.
\f[m_{4std} = kurtosis = \frac
  {\sum_i (i - centroid )^4 \times input(i)}
  {std^4 \sum_i input(i)}
\f]

Note that the kurtosis is often defined as the fourth cumulant divided
by the square root of the variance, which gives:
\f[kurtosis = \frac {m_4} {std^4} - 3\f]
This function does not include the `"-3"` term.


<!-- section 5 -->
<a name="descr_mel_bands"></a>
# 5 Descriptors based on the mel bands

The mel scale can be derived from the frequencies in hertz. The
conversion functions are in `rta_mel.h`. They are based on two
slightly different formulas, according to _HTK_ or the _Auditory Toolbox_
with the respective suffix `_htk` or `_slaney`.

The power spectrum is integrated into several bands, according again
to _HTK_ or the _Auditory Toolbox_.
The integration
window peak is 1 for _HTK_, while the sum of any channel is 1
for the _Auditory Toolbox_.

In order to reproduce the results of one of these tools, one must obviously
choose the desired variant among the whole computation process.

<a name="mel_bands"></a>
## 5.1 Mel bands

The mel bands integration is done in the magnitude (\f$abs\f$) or the
power (\f$abs^2\f$) domain, using respectively the function
`rta_spectrum_to_bands_abs` or `rta_spectrum_to_bands_square_abs`.
This respectively gives \f$mel\_bands = weights\_matrix \times spectrum\f$,
or \f$mel\_bands = (weights\_matrix \times \sqrt{spectrum})^2\f$.
The latter is the default (for _HTK_ and _Auditory Toolbox_) but it involves more
computation.

The matrix to multiply the power spectrum vector with, in order to
obtain the mel bands, must be computed beforehand, using the
`rta_spectrum_to_mel_bands_weights` function.

<a name="mfcc"></a>
## 5.2 Mel-frequency cepstral coefficients (MFCC)

First, the logarithm of the mel bands values is taken. In order to
avoid \f$log(0)\f$, one can add a very small value to the mel bands
values before taking the log.

Then a cepstrum is computed for the log of the mel spectrum, using a
discrete cosine transform _DCT_ of type II.
_HTK_ uses an orthogonal but not unitary transform, while the _Auditory Toolbox_
uses an orthogonal and unitary transform.
First, a weights matrix is constructed with the function `rta_dct_weights`
and it is then applied to the mel bands vector using the function `rta_dct`.
The coefficients obtained are the _MFCC_.

The _MFCC_, as any _DCT_ coefficients, are ordered by the order
of importance to model the spectrum. However, one can need to modify
them (for visualisation or further processing) using the liftering
functions in `rta_lifter.h`. These functions are provided for
the _HTK_ or _Auditory Toolbox_ compatibility.

<a name="delta_delta_mfcc"></a>
## 5.3 Delta and delta-delta MFCC

The function `rta_delta` computes a simple linear slope on a
sequence of fixed-rate sampled data (the frames). The _delta_ values
correspond to the mid-point frame (which is not used, by the way).
_It means that the delta values are not those of the last frame_.
Considering the filter-size, which is the size of the sequence taken into
account, the delay (in frames) introduced is half of the filter-size
(rounded down as it is always odd).

Note that the _HTK_ `DELTAWINDOW` variable is not the same as
the filter-size (the same applies for the `ACCWINDOW`
variable):
\f[filter\_size = 1 + 2 \times DELTAWINDOW\f]

Beforehand, a matrix of weights to multiply the sequence vector with
is constructed by the function `rta_delta_weights`.
A normalisation factor, computed by the function `rta_delta_normalization_factor`
gives _delta_ values, which are independent of the filter-size.
The normalisation factor can be applied directly to the weights matrix but
the rounding errors may be unacceptable when using
the simple floating-point precision.

Applying the _delta_ computation again gives the _delta-delta_ values,
_adding a new delay of half of the delta-delta filter-size_.


<!-- section 6 -->
<a name="style_guide"></a>
# 6 Style guide

<a name="prefix"></a>
## 6.1 Prefix

Any file, structure, function, type definition, enumeration
(enumerator and elements) and pre-compiler definition is prefixed with <br />
`rta_`, `RTA_`, `_rta_` or `_RTA_`.
There are two exceptions to this, to respect the common usage:
- `NULL` is defined in `rta_stdlib.h`
- mathematical constants in `rta_math.h` are prefixed with `M_`

<a name="declarations_definition"></a>
## 6.2 Declarations and definitions

A minimum number of `#define` statements is used, to minimise the global effects:
- the `const` keyword is used for values fixed at compilation time
- the `inline` keyword is used for in-lined functions

Moreover, the `static` keyword can be used to limit the definitions to a
particular file. <br />
The pre-compiler definitions are used to provide transparent wrappers for types
(see [6.4](#real_types)), functions (see [6.5](#memory_allocation)) or to
prevent multiple inclusions of files (see [6.3](#headers)).

<a name="headers"></a>
## 6.3 Headers

Any header prevents multiple inclusions, includes the `rta.h`
header and allows a C++ inclusion.
For a file named `rta_filename.h`, the content begins as:
```
#ifndef _RTA_FILENAME_H_
#define _RTA_FILENAME_H_ 1

#include "rta.h"

#ifdef __cplusplus
extern "C" {
#endif
```

and it ends as:
```
#ifdef __cplusplus
}
#endif

#endif /* _RTA_FILENAME_H_ */
```
<a name="real_types"></a>
## 6.4 Real types

Any function within _RTA_ uses the `rta_real_t` type.

`rta_real_t` is determined at compilation time by a
`#define` (not a `typedef`) in `rta_types.h` to
ensure a strict `float` or `double` replacement. (As a
such, `rta_real_t` can be used with `typedef` and
`sizeof` expression.) The functions from the standard header
`math.h` are redefined in `rta_math.h` to provide
corresponding functions for the `rta_real_t` type.

You can define the pre-compiler variable `RTA_REAL_TYPE` to
`RTA_FLOAT_TYPE`, `RTA_DOUBLE_TYPE` or
`RTA_LONG_DOUBLE_TYPE` for using respectively
`float`, `double` or `long double` type as the
effective representation of `rta_real_t`. Note that the
latter is untested. The default is to use `RTA_FLOAT_TYPE` if
`RTA_REAL_TYPE` is not defined into the user's
`rta_configuration.h`.

The same applies to the `rta_complex_t` type.
You can define the pre-compiler variable `RTA_COMPLEX_TYPE` to `RTA_FLOAT_TYPE`,
`RTA_DOUBLE_TYPE` or `RTA_LONG_DOUBLE_TYPE` for using respectively
`float complex`, `double complex` or `long double complex` type as the effective
representation of `rta_complex_t`.
Note that the latter is untested. The default is to use the same type as
`RTA_REAL_TYPE` if `RTA_COMPLEX_TYPE` is not defined into the user's
`rta_configuration.h`.

<a name="memory_allocation"></a>
## 6.5 Memory allocation

The arrays of real values are not allocated within the \rta\ library,
as they are manipulable outside of the library. They must be allocated
beforehand, as no memory allocation is done during a function's call
(except for the followings).

Specific structures are private: they are defined into the `*.c`
files, without the `_t` suffix. As an example, there is a
public declaration of `rta_fft_setup_t` in
`rta_fft.h`:
```
typedef struct rta_fft_setup rta_fft_setup_t;
```
And in `rta_fft.c`, there is a private definition of
`rta_fft_setup`, which may depend on the actual
implementation.

Two functions are provided for any private structure:
- a function post-fixed with `_setup_new` to allocate the structure by using the
`rta_malloc` function
- a function post-fixed with `_setup_delete` to release the structure by using
the `rta_free` function

Any function dealing with memory within _RTA_ uses the functions from
`rta_stdlib.h`. `rta_malloc`, `rta_realloc` and `rta_free` are pre-compiler
definitions that can be provided by the user's `rta_configuration.h`.
They must follow the declarations of respectively `malloc`, `realloc` and `free`
from the standard `stdlib.h`, which they are bound to if no user definition is
given.

