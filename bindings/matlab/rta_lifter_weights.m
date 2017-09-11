function [weights] = rta_lifter_weights (size, coefficient, type, mode)
% function [weights] = rta_lifter_weights (size, coefficient, type, mode)
%
% Create a vector (a <size> x 1 matrix) of <weights>. Only the <size>
% parameter is mandatory.
%
% <coefficient> is optional (0.6 by default, common for 'exp').
%    It must be > 0 for 'sin' (22 is common).
% <type> can be 'slaney' (default, same as 'exp') or 'htk' (same as
%    'sin').
% <mode> can be 'normal' (default) or 'inverse'
%

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
