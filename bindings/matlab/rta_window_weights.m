function [weights] = rta_window_weights (size, type, coefficient)
% function [weights] = rta_window_weights (size, type, coefficient)
%
% Create a vector (a <size> x 1 matrix) of <weights>:
%
% <type> can be 'hann' (default) or 'hamming'
% <coefficient> is optional (0.08 by default for a true hamming window)
%    and applies to the 'hamming' <type> only. 
%

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
