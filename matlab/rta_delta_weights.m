function [weights, normalisation_factor] = rta_delta_weights (size)
% function [weights, normalisation_factor] = rta_delta_weights (size)
%
% Create a vector (a <size> x 1 matrix) of <weights> to be used by the
% rta_delta_apply function.
% 
% The <normalisation_factor> can be used to get the slope of regularly
% sampled input values. One can directly normalise the <weights> vector by
% the <normalisation_factor> but rounding errors may occur.
%

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
