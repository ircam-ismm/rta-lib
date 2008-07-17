function [bands] = rta_spectrum_to_bands (spectrum, weights, bounds, type)
% function [bands] = rta_spectrum_to_bands (spectrum, weights, bounds, type)
%
% Apply <weights> to a real power <spectrum> vector. <weights> and <bounds>
% must be created by the rta_window_weights function.
%
% <type> can be:
%   'sqrabs' to integrate the power spectrum in abs^2 domain as:
%       <bands> = (<weights>' * sqrt(<spectrum>)) .^2
%   'abs' to integrate bands is abs domain, as:
%       <bands> = <weights>' * <spectrum>
%    
% If the <spectrum> size is even, only the lower part of the spectrum (plus
% the middle point) is used by this function, that is
% <spectrum>(1:end/2+1); if the <spectrum> size is odd, every point is used.
        
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
