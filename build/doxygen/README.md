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
-6 [Style guide](#style_guide)
  - 6.1 [Prefix](#prefix)
  - 6.2 [Declarations and definitions](#declarations_definitions)
  - 6.3 [Headers](#headers)
  - 6.4 [Real types](#real_types)
  - 6.5 [Memory allocation](#memory_allocation)

<!-- section 1 -->
<a name="generalities"></a>
# 1 Generalities

The _RTA_ library is frame-based, which means that for any vectors of
samples, a set of descriptors can be computed without adding any
delay. A noticeable exception to this is the _delta_ (and
_delta-delta_) computation, as it is based on more than one frame.

There are two variants of each function, one being <strong>\_stride</strong>
post-fixed. It allows to directly access interleaved data without
copying them (like a stereophonic samples vector). However, using a
big _stride_ value may lead to a bad usage of the memory cache.

Any index starts at 0.

Some functions require an initialisation before any processing, in
order to pre-calculate what will not depend on the incoming
frame. Every allocation must be done before anything else, outside of
the functions themselves except mentioned otherwise.

Some descriptors can be computed by several functions, and the results
may slightly differ for several reasons: the functions does not rely
on the same algorithms and the signal used for the computation may
differ (due to windowing, filtering, etc.). The auto-correlation from
_yin_ and from the _LPC_ are not the same and there is a
lot of ways to get the energy: from _yin_, as the sum of the
squares of the samples, from _LPC_, or as the first
_MFCC_ coefficient.


<!-- section 2 -->
<a name="configuration"></a>
# 2 Library configuration

A file named **rta\_configuration.h** must be present within your
sources in order to use the _RTA_ library. It is not included within
the _RTA_ source files. An empty file means that the defaults settings
are used for the compilation.

Instead of using the **malloc**, **realloc** and
**free** functions from the **stdlib.h**, one can
respectively define **rta\_malloc**, **rta\_realloc**
and **rta\_free**. (see [6.5](#memory_allocation))

The floating-point precision can be simple, double or long double,
according to the definition of **RTA\_REAL\_TYPE** to
respectively **RTA\_FLOAT\_TYPE**, **RTA\_DOUBLE\_TYPE** or
**RTA\_LONG\_DOUBLE\_TYPE**. The constants in
**rta\_float.h** are then redefined according to the proper type
from **float.h**. The same applies to the functions in
**rta\_math.h** from **math.h**. (see [6.4](#real_types))

Note that the long double precision is not supported when using the
Apple's VecLib by setting **RTA_USE_VECLIB** to 1.

<!-- INCLUDE FIGURE DATAFLOW -->

<!-- section 3 -->
<a name="descr_vector_samples"></a>
# 3 Descriptors based on a vector of samples

The vector of samples is characterised by its sample-rate and its
size. Moreover, to use sizes larger than those provided by the
sound card (e.g. for Fourier transforms), the hop-size gives the number
of samples between two consecutive vectors that overlap if the hop-size
is smaller than the vector-size.

<a name="resampling"></a>
## 3.1 Re-sampling

Down-sampling a signal is interesting to lower the further
computations, especially for the computational-intensive algorithms,
like _yin_. It can helps to concentrate on a range of the
spectrum where the information is pertinent for the analysis: the
_LPC_ (which is linear, as its name suggests it) finds poles
and zeros to fit the whole spectrum; to keep the results under 5 kHz,
we simply use a sample-rate of 11 kHz. If the original signal
sample-rate is 44 kHz, we can use the function
**rta\_downsample\_int\_mean** with a factor of 4. This function
implies a low-pass filtering. The function
**rta\_downsample\_int\_remove** can produce aliasing if the
original signal contains information above half of the resulting
sample-rate. Note that the resulting samples vector-size is smaller,
according to the given factor.

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
<!-- $$periodicity = 1 - \sqrt{absolute\_minimum}$$ -->

Before computing anything, a new _yin_ structure must be
allocated and filled with the **rta\_yin\_setup_new**
function. It will be released by the function
**rta\_yin\_setup\_delete**. (The auto-correlation result vector
must also be allocated beforehand, like any results vector.)

<a name="gain_integer_float"></a>
## 3.3 Gain and integer-float conversion

The samples are generally coded by floating-points number over 32
bits, within the range [-1.0, 1.0]. To convert them into 16 bits
signed integers (which is the format used by some systems, like
_HTK_), one can apply a simple gain of <!-- $2^{15}$ --> by multiplying
every sample.

<a name="pre_emphasis"></a>
## 3.4 Pre-emphasis

The pre-emphasis is a simple first-order difference between a current
sample and the previous one (weighted by a factor):
<!-- $$s(n) = s(n) - f \times s(n-1)$$ -->

It is often used for voice analysis with a factor of 0.97 as it
reduces the low frequencies while raising the high frequencies, thus
amplifying the contrast.

<a name="windowing"></a>
## 3.5 Windowing

If the samples vector-size is known, it is possible to pre-calculate
the weights that will be used to apply a given function. The function
**rta\_window\_hamming\_weights** computes a _Hamming_
window while the function **rta\_window\_hann\_weights** computes
a _von Hann_ window. These (or any weights vector) can be
applied with the **rta\_window\_apply** function.
The <strong>\_in\_place</strong> post-fixed functions change the input samples
vector values directly.

If the samples vector-size is not known in advance, one can still
apply the window using the **rta\_window\_rounded\_apply**
function. There is no interpolation, then. The weights vector indexes
are simply scaled and rounded: this is efficient but the rounding
error may be unacceptable if the size of the weights vector is too
small comparing with the samples vector size. It is also possible to
compute and apply a window on the fly, with the functions
**rta\_window\_hann\_apply** and
**rta\_window\_hamming\_apply**.

<a name="lpc"></a>
## 3.6 Linear predictive coefficients (LPC)

The **rta\_lpc** function calculates the linear predictive
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

Before computing a real _Fourier_ transform, a new real
_Fourier_ transform setup must be allocated and filled with the
**rta\_fft\_real\_setup\_new** function, with the type
**real\_to\_complex\_1d**. It will be released by the
function **rta\_fft\_setup\_delete**. The function
**rta\_fft\_execute** applies the _Fourier_ transform
to a samples vector.

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
<!-- $$ frequency = index \frac{sample\_rate}{transform\_size}$$ -->

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
**rta\_weighted\_moment\_1\_indexes**. The result unit is
<!-- $index$ --> (of the power spectrum, starting at 0). This function
returns also the input sum as it can be used in further calculations.
<!-- $$ m_1 = centroid = \frac {\sum_i i \times input(i)} -->
<!-- {\sum_i input(i)} $$ -->

<a name="spread_deviation"></a>
### 4.3.2 Spread and deviation

The spectral spread is the second central moment over the indexes
weighted by the vector of power spectrum values. It is computed by
**rta\_weighted\_moment\_2\_indexes**. The result unit is
<!-- $index^2$ --> (of the power spectrum, starting at 0).
<!-- $$ m_2 = spread = \frac {\sum_i (i - centroid )^2 \times input(i)} -->
<!-- {\sum_i input(i)} $$ -->

The standard deviation is <!-- $std = \sqrt{spread}$ -->.

<a name="skewness"></a>
### 4.3.3 Skewness

The spectral skewness is the third standard central moment over the
indexes weighted by the vector of power spectrum values. It is
computed by \texttt{rta\_std\_weighted\_moment\_3\_indexes}. The
result is without unit.
$$ m_{3std} = skewness = \frac {\sum_i (i - centroid )^3 \times input(i)}
{std^3 \sum_i input(i)} $$

<a name="kurtosis"></a>
### 4.3.4 Kurtosis

The spectral kurtosis is the fourth standard central moment over the
indexes weighted by the vector of power spectrum values. It is
computed by **rta\_std\_weighted\_moment\_4\_indexes**. The
result is without unit.
<!-- $$ m_{4std} = kurtosis = \frac {\sum_i (i - centroid )^4 \times input(i)} -->
<!-- {std^4 \sum_i input(i)} $$ -->

Note that the kurtosis is often defined as the fourth cumulant divided
by the square root of the variance, which gives $kurtosis =
<!-- \frac{m_4}{std^4} - 3$. -->
This function does not include the <!-- ``$- 3$'' term -->.


<!-- section 5 -->
<a name="descr_mel_bands"></a>
# 5 Descriptors based on the mel bands

<a name="mel_bands"></a>
## 5.1 Mel bands

<a name="mfcc"></a>
## 5.2 Mel-frequency cepstral coefficients (MFCC)

<a name="delta_delta_mfcc"></a>
## 5.3 Delta and delta-delta MFCC


<!-- section 6 -->
<a name="style_guide"></a>
# 6 Style guide

<a name="prefix"></a>
## 6.1 Prefix

<a name="declarations_definition"></a>
## 6.2 Declarations and definitions

<a name="headers"></a>
## 6.3 Headers

<a name="real_types"></a>
## 6.4 Real types

<a name="memory_allocation"></a>
## 6.5 Memory allocation

