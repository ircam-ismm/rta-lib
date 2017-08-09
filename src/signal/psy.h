/**
 *
 * @file psy.h
 * @author Norbert.Schnell@ircam.fr
 * 
 * @brief yin-based PSOLA analysis
 * 
 * Copyright (C) 2009-2014 by IRCAM – Centre Pompidou, Paris, France.
 * All rights reserved.
 * 
 */
#define psyAna_malloc malloc
#define psyAna_realloc realloc
#define psyAna_free free

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
} CandidateT;

typedef struct TrackingStateSt
{
  CandidateT candidates[MAX_CANDIDATES];
  int numCandidates;
  double time;
  double energy;
  double ac1;
} TrackingStateT;

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
  
  TrackingStateT trackingStates[NUM_TRACKING_STATES];
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
} PsyAnaT;

void psyAna_init(PsyAnaT *self);
void psyAna_deinit(PsyAnaT *self);
void psyAna_reset(PsyAnaT *self, double minFreq, double maxFreq, double sampleRate, int maxInputVectorSize, int downSamplingExp);
void psyAna_setCallback(PsyAnaT *self, void *receiver, int (*callback)(void *receiver, double time, double freq, double energy, double ac1, double voiced));
void psyAna_setThresholds(PsyAnaT *self, double yinThreshold, double noiseThreshold);
int psyAna_calculateInputVector(PsyAnaT *self, float *in, int vectorSize, int vectorStride);
void psyAna_finalize(PsyAnaT *self);
