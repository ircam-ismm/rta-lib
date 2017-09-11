function [weights] = rta_dct_weights (input_size, order, type)
% function [weights] = rta_dct_weights (input_size, order, type)
%
% Create a matrix (a <input_size> x <order> matrix) of <weights>.
%
% <type> can be:
%    'slaney': type II, orthogonal and unitary (Auditory Toolbox like) (default)
%    'htk': type II, orthogonal but not unitary (HTK like)
%

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
