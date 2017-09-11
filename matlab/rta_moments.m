function [moments, input_sum] = rta_moments (input, max_order, type)
%  function [moments, input_sum] = rta_moments (input, max_order, type)
%
% Compute the moments of inertia of the <input> vector up to order <max_order>.
% Moments whose order is > 1 are centered and normalised.
% Moments whose order is > 2 can be standardised by the deviation (square
% root of the second moment), depending on the <type>:
%   'std' (default) means that the moments are standardised
%   'nostd'  means that the moments are not standardised
%       

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
