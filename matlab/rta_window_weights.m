function [weights] = rta_window_weights (type, size, coefficient)
% function [weights] = rta_window_weights (type, size, coefficient)
%
% Create a vector (a <size> x 1 matrix) of <weights>:
%
% <type> can be 'hann' or 'hamming'
% <coefficient> is optional (0.08 by default for a true hamming window)
%    and applies to the 'hamming' <type> only. 
%
% <setup> is a pointer the the newly created setup
%
% Use rta_yin_setup_delete(setup) to release the memory of this setup
% when every rta_yin calculation is done.

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
