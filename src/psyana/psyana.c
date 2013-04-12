/**
 *
 * @file psyana.c
 * @author Norbert.Schnell@ircam.fr
 * 
 * @brief yin-based PSOLA analysis
 * 
 * Copyright (C) 2009-2011 by IMTR IRCAM – Centre Pompidou, Paris, France.
 * All rights reserved.
 * 
 */
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <float.h>

#include "psyana.h"

#ifndef WIN_VERSION

#include "vDSP.h"
#define yinAutocorr(i, c, n, m) conv((i), 1, (i), 1, (c), 1, (n), (m));

#else

static void
yinAutocorr(float *in, float *corr, int n, int m)
{
  int i, j;
  
  for(i = 0; i < n; i++)
  {
    float c = 0.0;
    
    for(j = 0; j < m; j++)
      c += in[j] * in[i + j];
    
    corr[i] = c;
  }
}

#endif

#define REF_FREQ 440.0
#define REF_PITCH 6900.0

/* convert period to MIDI cent and back */
#define periodToMidiCent(x, rp, ps) ((rp) - (1731.23404906676 * log((x) * (ps))))

#define MIN_ANALYSIS_LAG 10
#define MAX_LAG_HIGH 128

static int
estimateCandidates(float *input, float *corrBuffer, double maxPeriod, double minPeriod, double noiseThreshold, CandidateT *candidates)
{
  int numCandidates = 1;
  int nCorr = (int)ceil(maxPeriod);
  double x, xn;
  double energy;
  double diffLeft, diff, diffRight, sum;
  int i;
  
  /* init absolute minimum */
  candidates[0].period = maxPeriod;
  candidates[0].normDiff = 1.0;

  /* auto-correlation */
  yinAutocorr(input, corrBuffer, nCorr, nCorr);
  
  /* diff[0] */
  x = input[0];
  xn = input[nCorr];
  energy = corrBuffer[0] + xn * xn - x * x;
  diffLeft = 0.0;
  diff = 0.0;
  diffRight = corrBuffer[0] + energy - 2 * corrBuffer[1];
  sum = 0.0;
  
  /* diff[1] */
  x = input[1];
  xn = input[1 + nCorr];
  energy += xn * xn - x * x;
  diffLeft = diff;
  diff = diffRight;
  diffRight = corrBuffer[0] + energy - 2 * corrBuffer[2];
  sum = diff;
  
  /* minimum difference search */
  for(i = 2; i < nCorr - 1; i++)
  {
    x = input[i];
    xn = input[i + nCorr];
    energy += xn * xn - x * x;
    diffLeft = diff;
    diff = diffRight;
    diffRight = corrBuffer[0] + energy - 2 * corrBuffer[i + 1];
    sum += diff;
    
    /* local minimum */
    if(diff < diffLeft && diff < diffRight)
    {
      double a = diffLeft + diffRight - 2.0 * diff;
      double b = 0.5 * (diffRight - diffLeft);
      double m = diff - (b * b) / (2.0 * a);
      double period = i - b / a;
      double normDiff = i * m / sum;

      if(period >= minPeriod)
      {
        if(normDiff < 0.0)
          normDiff = 0.0;
        
        /* store candidate */
        if(normDiff <= noiseThreshold)
        {
          candidates[numCandidates].period = period;
          candidates[numCandidates].normDiff = normDiff;
          numCandidates++;
        }
        
        /* get absolute minimum */
        if(normDiff < candidates[0].normDiff)
        {
          candidates[0].period = period;
          candidates[0].normDiff = normDiff;
        }
      }
      
      if(numCandidates >= MAX_CANDIDATES)
        break;
    }
  }
  
  return numCandidates;
}

static void
traceBackAndForth(PsyAnaT *self, TrackingStateT *trackingStates, int trackingIndex, CandidateT *candidate, int candidateIndex)
{
  TrackingStateT *backwardState = trackingStates + (trackingIndex + NUM_TRACKING_STATES - 1) % NUM_TRACKING_STATES;
  double nearestDiff = DBL_MAX;
  int backwardIndex = 0;
  int i;
  
  /* get closest pitch in last state */
  for(i = 1; i < backwardState->numCandidates; i++)
  {
    double pitch = periodToMidiCent(candidate->period, self->referencePitch, self->referencePeriodScale);
    double diff = fabs(pitch - backwardState->candidates[i].pitch);
    
    /* set candidate pitch */
    candidate->pitch = pitch;    
    
    /* exit loop when nearest has passed */
    if(diff > nearestDiff)
      break;
      
    nearestDiff = diff;
    backwardIndex = i;
  }

  candidate->forward = candidate->backward = candidate->normDiff;

  if(nearestDiff <= self->pitchDiffTolerance)
  {
    CandidateT *backwardCandidate = backwardState->candidates + backwardIndex;
    
    candidate->backwardIndex = backwardIndex;
    backwardCandidate->forwardIndex = candidateIndex;

    /* propagate forward */
    candidate->forward += 1.0 * backwardCandidate->forward;
    candidate->forward *= 0.5;

    /* propagate backward */
    backwardCandidate->backward += 1.0 * candidate->normDiff;
    backwardCandidate->backward *= 0.5;
  }
}

static int
getBestCandidate(PsyAnaT *self, TrackingStateT *trackingStates, int trackingIndex)
{
  TrackingStateT *state = trackingStates + trackingIndex;
  double yinThreshold = self->yinThreshold;
  int i;
  
  yinThreshold += state->candidates[0].normDiff;
  
  for(i = 1; i < state->numCandidates; i++)
  {
    CandidateT *candidate = state->candidates + i;
    double normDiff = candidate->normDiff;
    
    /* finish backward tracking */
    if(candidate->forwardIndex > 0)
    {
      TrackingStateT *forwardState = trackingStates + (trackingIndex + NUM_TRACKING_STATES + 1) % NUM_TRACKING_STATES;
      CandidateT *forwardCandidate = forwardState->candidates + candidate->forwardIndex;
      
      candidate->backward += 1.0 * forwardCandidate->backward;
      candidate->backward *= 0.5;
      
      if(candidate->backward < normDiff)
        normDiff = candidate->backward;
    }

    if(candidate->forward < normDiff)
      normDiff = candidate->forward;

    /* classic yin: look for first minimum under threshold */
    if(normDiff < yinThreshold)
      return i;
  }
  
  return 0;
}

static double
estimateAndTraceCandidates(PsyAnaT *self, TrackingStateT *trackingStates, int trackingIndex, double time)
{
  TrackingStateT *state = trackingStates + trackingIndex;
  double yinThreshold = self->yinThreshold;
  double norm = 1.0 / ceil(self->maxPeriod);
  int yinIndex = 0;
  double period;
  int i;
  
  /* calculate candidates (assigns period and normDiff) */
  state->numCandidates = estimateCandidates(self->inputBuffer, self->corrBuffer, self->maxPeriod, self->minPeriod, self->noiseThreshold, state->candidates);

  state->time = time;  
  state->energy = sqrt(self->corrBuffer[0] * norm);
  state->ac1 = self->corrBuffer[1] * norm;
    
  /* bias threshold by absolute minimum */
  yinThreshold += state->candidates[0].normDiff;
  
  /* get yin index */
  for(i = 1; i < state->numCandidates; i++)
  {
    CandidateT *candidate = state->candidates + i;
    double normDiff = candidate->normDiff;
    
    /* init tracing */
    candidate->pitch = DBL_MAX;
    candidate->forward = 0.0;
    candidate->forwardIndex = 0;
    candidate->backward = 0.0;
    candidate->backwardIndex = 0;
    
    /* perform bakward and forward tracing */
    traceBackAndForth(self, trackingStates, trackingIndex, state->candidates + i, i);
    
    if(candidate->forward < normDiff)
      normDiff = candidate->forward;
    
    if(yinIndex == 0 && normDiff < yinThreshold)
      yinIndex = i;
  }
  
  if(state->candidates[yinIndex].normDiff <= self->noiseThreshold)
    period = state->candidates[yinIndex].period;
  else
    period = self->lastPeriod;
  
  self->lastPeriod = period;
  
  return period;
}

static int
reportState(PsyAnaT *self, TrackingStateT *trackingStates, int trackingIndex)
{
  TrackingStateT *state = trackingStates + trackingIndex;
  double reportTime, reportFreq;
  int cont = 1;
  
  if(state->time >= 0.0)
  {
    int candidateIndex = getBestCandidate(self, trackingStates, trackingIndex);
    CandidateT *reportCandidate = state->candidates + candidateIndex;
    double period, voiced;
    
    if(reportCandidate->normDiff < self->noiseThreshold)
      period = reportCandidate->period;
    else
      period = self->lastReportPeriod;

    voiced = 1.0 - sqrt(reportCandidate->normDiff);

    reportTime = state->time * self->downSamplingRatio * self->samplePeriod;
    reportFreq = self->sampleRate / (period * self->downSamplingRatio);

    if(self->callback != NULL)
      cont = (*self->callback)(self->receiver, reportTime, reportFreq, state->energy, state->ac1, voiced);
    
    self->numOutput++;
    self->lastReportPeriod = period;
  }
  
  return cont;
}

static int 
dummyCallback(void *receiver, double time, double freq, double energy, double ac1, double voiced)
{
  return 0;
}

void
psyAna_init(PsyAnaT *self)
{
  /* input buffers */
  self->inputBuffer = NULL;
  self->maxInputBufferSize = 0;
  self->inputBufferSize = 0;
  
  self->inputTime = 0.0;
  self->inputFill = 0;
  
  self->inputTime = 0;
  self->inputFill = 0;
  self->outputTime = 0.0;
  
  self->downSamplingExp = 0;
  self->downSampling = 1; 
  self->downSamplingRatio = 1.0;
  self->downSamplingScaling = 1.0;
  
  /* autocorrelation buffer */
  self->corrBuffer = NULL;
  self->maxCorrBufferSize = 0;
  self->corrBufferSize = 0;

  self->absMaxPeriod = 0.0;
  self->absMinPeriod = 0.0;
  self->maxPeriod = 0.0;
  self->minPeriod = 0.0;
  
  self->yinThreshold = 0.1024;
  self->noiseThreshold = 0.3025;
  
  self->lastPeriod = 0.0;
  self->lastReportPeriod = 0.0;
  
  self->referencePitch = 0.0;
  self->referencePeriodScale = 1.0;
  self->pitchDiffTolerance = 100.0; /* in cent */

  self->sampleRate = 1000.0;
  self->samplePeriod = 1.0;
  self->maxInputVectorSize = 0;
    
  self->receiver = NULL;
  self->callback = dummyCallback;
  self->numOutput = 0;
}

void
psyAna_deinit(PsyAnaT *self)
{
  if(self->inputBuffer != NULL)
    psyAna_free(self->inputBuffer);
  
  if(self->corrBuffer != NULL)
    psyAna_free(self->corrBuffer);
}

void
psyAna_reset(PsyAnaT *self, double minFreq, double maxFreq, double sampleRate, int maxInputVectorSize, int downSamplingExp)
{
  double minMinPeriod = 2.0;
  double maxMaxPeriod = 10000.0;
  double absMinPeriod, absMaxPeriod;
  int i;

  if(downSamplingExp < 0)
    downSamplingExp = 0;
  else if(downSamplingExp > 3)
    downSamplingExp = 3;

  self->downSamplingExp = downSamplingExp;
  self->downSampling = (1 << downSamplingExp); 
  self->downSamplingRatio = (double)self->downSampling;
  self->downSamplingScaling = 1.0 / self->downSamplingRatio;
    
  absMinPeriod = (sampleRate / maxFreq) * self->downSamplingScaling;
  
  if(absMinPeriod < minMinPeriod)
    absMinPeriod = minMinPeriod;

  absMaxPeriod = (sampleRate / minFreq) * self->downSamplingScaling;

  if(absMaxPeriod > maxMaxPeriod)
    absMaxPeriod = maxMaxPeriod;
  
  if(absMinPeriod > absMaxPeriod)
    absMinPeriod = absMaxPeriod;
    
  self->corrBufferSize = (int)ceil(2 * absMaxPeriod) + 2;

  if(self->corrBufferSize > self->maxCorrBufferSize)
  {
    self->corrBuffer = (float *)psyAna_realloc(self->corrBuffer, sizeof(float) * self->corrBufferSize);
    self->maxCorrBufferSize = self->corrBufferSize;
  }

  self->inputBufferSize = 3 * self->corrBufferSize + (maxInputVectorSize >> self->downSamplingExp) + self->downSampling;
  
  if(self->inputBufferSize > self->maxInputBufferSize)
  {
    self->inputBuffer = (float *)psyAna_realloc(self->inputBuffer, sizeof(float) * self->inputBufferSize);
    self->maxInputBufferSize = self->inputBufferSize;
  }

  self->absMinPeriod = absMinPeriod;
  self->absMaxPeriod = absMaxPeriod;

  self->minPeriod = absMinPeriod;
  self->maxPeriod = absMaxPeriod;

  self->sampleRate = sampleRate;
  self->samplePeriod = 1000.0 / sampleRate;
  self->maxInputVectorSize = maxInputVectorSize;

  self->inputTime = 0;
  self->inputFill = 0;
  self->outputTime = 0.0;

  for(i = 0; i < NUM_TRACKING_STATES; i++)
  {
    self->trackingStates[i].numCandidates = 0;
    self->trackingStates[i].time = -DBL_MAX;
    self->trackingStates[i].energy = 0.0;
    self->trackingStates[i].ac1 = 0.0;
  }
  
  self->trackingIndex = 0;
  
  self->lastPeriod = absMaxPeriod; /* start with expected max period */
  self->lastReportPeriod = absMaxPeriod;
  self->referencePitch = REF_PITCH;
  self->referencePeriodScale = REF_FREQ / sampleRate * self->downSamplingRatio;
    
  self->numOutput = 0;
}

void 
psyAna_setCallback(PsyAnaT *self, void *receiver, int (*callback)(void *receiver, double time, double freq, double energy, double ac1, double voiced))
{
  self->receiver = receiver;
  self->callback = callback;
}

void
psyAna_setThresholds(PsyAnaT *self, double yinThreshold, double noiseThreshold)
{
  if(yinThreshold < 0.0)
    yinThreshold = 0.0;
  else if(yinThreshold > 1.0)
    yinThreshold = 1.0;
  
  yinThreshold = 1.0 - yinThreshold;
  
  self->yinThreshold = yinThreshold * yinThreshold;
  
  if(noiseThreshold < 0.0)
    noiseThreshold = 0.0;
  else if(noiseThreshold > 1.0)
    noiseThreshold = 1.0;
  
  noiseThreshold = 1.0 - noiseThreshold;
  
  self->noiseThreshold = noiseThreshold * noiseThreshold;
}

int
psyAna_calculateInputVector(PsyAnaT *self, float *in, int vectorSize, int vectorStride)
{
  int downVectorSize = vectorSize >> self->downSamplingExp;
  float *inputBuffer = self->inputBuffer + self->inputFill;
  double outputTime = self->outputTime;
  int maxTime;
  int i, j;
    
  if(downVectorSize > 0)
  {  
    switch(self->downSamplingExp)
    {
      case 3:
      {
        for(i = 0, j = 0; i < downVectorSize; i++, j += 8 * vectorStride)
          inputBuffer[i] = 0.125 * (in[j] + in[j + vectorStride] + 
                                    in[j + 2 * vectorStride] + in[j + 3 * vectorStride] + 
                                    in[j + 4 * vectorStride] + in[j + 5 * vectorStride] + 
                                    in[j + 6 * vectorStride] + in[j + 7 * vectorStride]);
        break;
      }
        
      case 2:
      {
        for(i = 0, j = 0; i < downVectorSize; i++, j += 4 * vectorStride)
          inputBuffer[i] = 0.25 * (in[j] + in[j + vectorStride] + 
                                   in[j + 2 * vectorStride] + in[j + 3] * vectorStride);
        
        break;
      }
        
      case 1:
      {
        for(i = 0, j = 0; i < downVectorSize; i++, j += 2 * vectorStride)
          inputBuffer[i] = 0.5 * (in[j] + in[j + vectorStride]);
        
        break;
      }
        
      default:
      {
        for(i = 0, j = 0; i < downVectorSize; i++, j += vectorStride)
          inputBuffer[i] = in[j];
      }
    }      
  }
  else
  {
    float sum = 0.0;
    
    for(i = 0; i < vectorSize * vectorStride; i += vectorStride)
      sum += in[i];
        
    inputBuffer[0] = sum / (float)vectorSize;
  }
    
  self->inputFill += downVectorSize;
  maxTime = self->inputTime + self->inputFill - 2 * (int)self->absMaxPeriod;
  
  while(maxTime > (int)outputTime)
  {
    double period = estimateAndTraceCandidates(self, self->trackingStates, self->trackingIndex, outputTime);
    int shift;    
    
    /* advance tracking index */
    self->trackingIndex = (self->trackingIndex + 1) % NUM_TRACKING_STATES;
    
    if(reportState(self, self->trackingStates, self->trackingIndex) == 0)
      return 0;
    
    /* advance time */
    outputTime += period;
    
    /* shift input buffers */
    shift = ((int)(outputTime - self->inputTime) - 1); /* always keep one input sample */
    
    if(shift > 0)
    {
      for(i = 0; i < self->inputFill; i++)
        self->inputBuffer[i] = self->inputBuffer[i + shift];
      
      self->inputFill -= shift;
      self->inputTime += shift;
      maxTime -= shift;
    }
  }
  
  self->outputTime = outputTime;

  return vectorSize;
}
