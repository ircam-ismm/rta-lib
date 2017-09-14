/**
 * @file rta_psy.h
 * @author Norbert.Schnell@ircam.fr
 * @ingroup rta_signal
 *
 * @brief yin-based PSOLA analysis
 *
 * @copyright
 * Copyright (C) 2009-2014 by IRCAM – Centre Pompidou, Paris, France.
 * All rights reserved.
 *
 * License (BSD 3-clause)
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met:
 *
 * 1. Redistributions of source code must retain the above copyright notice,
 *    this list of conditions and the following disclaimer.
 * 2. Redistributions in binary form must reproduce the above copyright
 *    notice, this list of conditions and the following disclaimer in the
 *    documentation and/or other materials provided with the distribution.
 * 3. Neither the name of the copyright holder nor the names of its
 *    contributors may be used to endorse or promote products derived from
 *    this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
 * POSSIBILITY OF SUCH DAMAGE.
 */

#ifndef _RTA_PSY_H_
#define _RTA_PSY_H_ 1

#define rta_psy_malloc malloc
#define rta_psy_realloc realloc
#define rta_psy_free free

#define MAX_DOWN_SAMPLING_EXP 3
#define MAX_DOWN_SAMPLING (1 << MAX_DOWN_SAMPLING_EXP)

#define MAX_CANDIDATES 64
#define NUM_TRACKING_STATES 3

typedef struct CandidateSt
{
  double period;
  double pitch;
  double normDiff;
  double forward;
  int forwardIndex;
  double backward;
  int backwardIndex;
} rta_psy_candidate_t;

typedef struct TrackingStateSt
{
  rta_psy_candidate_t candidates[MAX_CANDIDATES];
  int numCandidates;
  double time;
  double energy;
  double ac1;
} rta_psy_tracking_state_t;

typedef struct PsyAnaSt
{
  /* pitch analysis stuff */
  float *inputBuffer; /* downsampled input buffer */
  int maxInputBufferSize; /* maximum size of input buffer */
  int inputBufferSize; /* current size of input buffer */

  int inputTime; /* current input time */
  int inputFill; /* current input fill */
  double outputTime; /* float index of next analysis frame */

  int downSampling;
  int downSamplingExp;
  double downSamplingRatio;
  double downSamplingScaling;

  double absMinPeriod; /* absolute minimum analysis period (= sample rate / minFreq) */
  double absMaxPeriod; /* absolute maximum analysis period (= sample rate / maxFreq) */
  double minPeriod; /* current minimum analysis period (last sure period / ambitus) */
  double maxPeriod; /* current maximum analysis period (last sure period * ambitus) */

  float *corrBuffer; /* temporary correlation coefficient buffer */
  int corrBufferSize; /* current maximum correlation window size (>= max period) */
  int maxCorrBufferSize; /* maximum correlation buffer size */

  double yinThreshold; /* yin normalized difference threshold (default 0.1024) */
  double noiseThreshold; /* yin normalized difference threshold (default 0.3025) */

  rta_psy_tracking_state_t trackingStates[NUM_TRACKING_STATES];
  int trackingIndex;

  double lastPeriod; /* latest period */
  double lastReportPeriod; /* latest reported period */
  double pitchDiffTolerance;

  double referencePitch;
  double referencePeriodScale;

  double sampleRate; /* sample rate */
  double samplePeriod; /* 1 / sample rate */
  int maxInputVectorSize; /* tick size */

  void *receiver;
  int (*callback)(void *obj, double time, double freq, double energy, double ac1, double voiced);
  int numOutput;
} rta_psy_ana_t;

void rta_psy_init(rta_psy_ana_t *self);
void rta_psy_deinit(rta_psy_ana_t *self);
void rta_psy_reset(rta_psy_ana_t *self, double minFreq, double maxFreq, double sampleRate, int maxInputVectorSize, int downSamplingExp);
void rta_psy_set_callback(rta_psy_ana_t *self, void *receiver, int (*callback)(void *receiver, double time, double freq, double energy, double ac1, double voiced));
void rta_psy_set_thresholds(rta_psy_ana_t *self, double yinThreshold, double noiseThreshold);
int rta_psy_calculate_input_vector(rta_psy_ana_t *self, float *in, int vectorSize, int vectorStride);
void rta_psy_finalize(rta_psy_ana_t *self);

#endif  /* _RTA_PSY_H_ */
