function [y, final_state] = ...
    rta_onepole (filter_type, f0, x, initial_state, dim)
% function [y, final_state] = ...
%    rta_onepole (filter_type, f0, x, initial_state, dim)
%
% filtertype can be:
%   'lowpass': ylp(n) = f0 * x(n) - (f0 - 1) * y(n-1)
%   'highpass': yhp(n) = x(n) - ylp(n)

% f0 is the cutoff frequency, normalised by the nyquist frequency.
%
% If 'initial_state' is empty or not provided, it is supposed to be [0,0].
%
% Use dimension == 1 for computation along the columns and
% dimension == 2 for computation along the rows. (The default is to
% choose dimension 1 if the number of rows is > 1.)

    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
