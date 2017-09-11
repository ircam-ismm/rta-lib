function [weights, bounds] = rta_bands_weights (...
    spectrum_size, bands_nb, type, sample_rate, min_freq, max_freq)
% function [weights, bounds] = rta_bands_weights (...
%    spectrum_size, bands_nb, type, sample_rate, min_freq, max_freq)
%
% If <spectrum_size> is even, the value <spectrum_size> / 2 + 1 is used
% instead (to process only the real spectrum); if <spectrum_size> is odd,
% it is inchanged.
%
% Create a matrix (<spectrum_size> x <bands_nb>) of weights to get the
% mel bands from a signal spectrum. See rta_spectrum_to_bands.
%
% <type> can be:
%   'slaney' for Auditory Toolbox style (constant sum per channel) (default)
%   'htk' for HTK style (constant maximum per channel)
% <sample_rate> is the input power spectrum corresponding sample rate in Hertz (default 44100.)
% <min_freq> is the output minimum frequency in Hertz (default 0.)
% <max_freq> is the output maximum frequency in Hertz (default 22050)

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
