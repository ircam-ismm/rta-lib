function [output] = rta_delta_apply (input, weights, scale)
% function [output] = rta_delta_apply (input, weights, scale)
%
% Apply a <weights> vector to an <input> vector. <weights> and <input>
% must have the same number of columns. If <scale> is provided, every
% <output> value is normalised by this factor.
%
% The <weights> and <scale> can be computed byt the rta_delta_window
% function.
%
    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
