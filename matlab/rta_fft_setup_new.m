function [setup] = rta_fft_setup_new (real_input, fft_size, scale)
%  function [setup] = rta_fft_setup_new (real_input, fft_size, scale)
%
% Create a new setup for the fft algorithm from the RTA library.
%
% <fft> size must be at least <real_input> size
% <scale> is optional and default value is 1. (no scale)
%
% <setup> is a pointer the the newly created setup
%
% Use rta_fft_setup_delete(setup) to release the memory of this setup
% when every rta_fft calculation is done.

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
