function [U, S, V] = rta_svd (A)
% function [U, S, V] = rta_svd(A)
%
% Singular Value Decomposition of A such as
% A = U * S * V'
%
% V calculation may be skipped using [U,S] = rta_svd(A)
% U and V calculations may be skipped using S = rta_svd(A)

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
