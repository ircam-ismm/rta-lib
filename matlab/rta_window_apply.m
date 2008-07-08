function [output] = rta_window_apply (input, window)
% function [output] = rta_window_apply (input, window)
%
% Apply a <window> vector to an <input> vector. If <window> and <input>
% sizes differ, then the <window> is rounded (by choosing the closest
% indexes). <output> size is always the same as <input>.
%
% A window can be created by the rta_window_weights function.
%
    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
