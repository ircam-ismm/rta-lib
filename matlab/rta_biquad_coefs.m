function [b, a] = rta_biquad_coefs (filter_type, f0, q, gain)
% function [b, a] = rta_biquad_coefs (filter_type, f0, q, gain)
%
% Compute the biquad coefficients b and a (to be used by rta_biquad and
% matlab's filter functions), such as:
% y(n) = b(1) x(n) + b(2) x(n-1) + b(3) x(n-2) 
%                  - a(2) y(n-1) - a(3) y(n-2) 
% Note that a(1) is always 1.
%
% filtertype can be:
%   'lowpass': H(s) = 1 / (s^2 + s/q + 1)
%   'highpass': H(s) = s^2 / (s^2 + s/q + 1)
%   'bandpass_cst_skirt': H(s) = s / (s^2 + s/q + 1)
%        (The peak gain is q*gain) 
%   'bandpass_cst_peak': H(s) = (s/q) / (s^2 + s/q + 1)
%        (The peak gain is gain) 
%   'notch': H(s) = (s^2 + 1) / (s^2 + s/q + 1)
%   'allpass': H(s) = (s^2 - s/q + 1) / (s^2 + s/q + 1)
%   'peaking': H(s) = (s^2 + s*(g/q) + 1) / (s^2 + s/(g*q) + 1),
%        with g = sqrt(gain)
%   'lowshelf': H(s) = g * (s^2 + (sqrt(g)/q)*s + g)/
%                         (g*s^2 + (sqrt(g)/q)*s + 1)
%        with g = sqrt(gain)
%   'highshelf': H(s) = g * (g*s^2 + (sqrt(g)/q)*s + 1)/
%                            (s^2 + (sqrt(g)/q)*s + g)
%        with g = sqrt(gain)
%    
% f0 is the cutoff frequency, normalised by the nyquist frequency.
%
% q must be > 0. and is generally >= 0.5 for audio filtering.    
% q <= 1./sqrt(2.) is the limit for monotonic response for lowpass,
% highpass, lowshelf and highshelf types.
%
% gain is linear.
%
% This coefficients can be used by rta_biquad and matlab's filter
% functions.
    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
