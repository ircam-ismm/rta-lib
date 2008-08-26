function [V, M] = rta_var(input, bias, dimension, accuracy)
% function [V, M] = rta_var(input, bias, dimension, accuracy)
%
% Compute the variance V and the mean M along the given dimension of
% input (of size MxN).
%
% V is normalised by (N-1) if bias == 0. This is the default and gives an
% unbiased estimator of the variance of the population from which input
% is drawn. V is normalised by N if bias == 1 and produces the second
% moment of the sample about its mean.
%
% Use dimension == 1 (default) for computation along the comumns and
% dimension == 2 for computation along rows.
%
% accuracy can be
%   'accurate' (default): mean is first calculated, then the variance,
%     according to var(X) = E( (x - mean(x))^2 )
%   'fast': mean and variance are calculated together, and variance is
%     then var(X) = E(x^2) - mean(x)^2
    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
