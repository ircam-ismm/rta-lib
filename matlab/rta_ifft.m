function [output] = rta_ifft (input, setup)
% function [output] = rta_ifft (input, setup)
%
% Compute the inverse FFT of the complex <input> signal using the RTA library.
%
% <input> is a complex signal vector which size must not be greater than the
% fft size used to define the setup (see rta_ifft_setup_new).
%
    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
