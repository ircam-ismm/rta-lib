function [f0, energy, periodicity, ac1_ac0, autocorrelation] = ...
    rta_yin (input, threshold, min_freq, sample_rate, setup)
% function [f0, energy, periodicity, ac1_ac0, autocorrelation] = ...
%    rta_yin (input, threshold, min_freq, sample_rate, setup)
%
% Process the <input> signal through the yin algorithm from the RTA
% library.
%
% <input> is a signal vector which size (input_size) must satisfy
%    input_size > sample_rate/min_freq
%    for good results it should satisfy
%    input_size >= 2 * sample_rate/min_freq
% <threshold> under which to search the minima. It must be in
%    [0., 1.] and  0.1 is a common value. This is the original threshold
%    (and will raise if no minimum is found).
%    threshold == (1. - confidence)^2
% <min_freq> is the minimum f0 searched in Hz
% <sample_rate> is input sample-rate in Hz
% <setup> is created by rta_yin_setup_new(yin_max_mins)
%
% <f0> is the fundamental frequency in Hz found by the yin algorithm
% <energy> of input equals to
%    sqrt(autocorrelation[1]/(sample_rate/min_freq))
%    max_lag = sample_rate/min_freq
% <ac1_ac0> is the ratio between the two first autocorrelation
%    coefficients:
%    ac1_ac0 = autocorrelation[2]/autocorrelation[1]
% <autocorrelation> coefficients up to order 
%    (input_size - (sample_rate/min_freq) - 2)

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
