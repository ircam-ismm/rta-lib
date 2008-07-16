function [coefficients, error, autocorrelation] = rta_lpc (input, order)
% function [coefficients, error, autocorrelation] = rta_lpc (input, order)
%
% Compute the LPC <coefficients> on the <input> vector. Note that the
% number of coefficients is (<order> + 1), as the <autocorrelation> size.
%
% If the original spectrum is:
%  20 * log( abs( fft( <input>, 512 ) ) )
% Then the LPC estimated spectrum is
%  -20 * log( abs( fft( <coefficients> ./ ( sqrt( abs( <error> ) ) ), 512) ) )
%

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
