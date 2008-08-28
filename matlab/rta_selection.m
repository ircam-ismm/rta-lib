function [value] = rta_selection (input, index, dim)
% function [value] = rta_selection (input, index, dim)
%
% Selection of an index, as if the input was sorted. If the
% given index is not an integer, the weighted mean of the two
% adjacent indexes is returned. The medians of the columns can be
% computed as:
%
% median = rta_selection(input, rows_nb * 0.5, 1)
% where rows_nb is the number of rows
% 
% Use dimension == 1 for computation along the columns and
% dimension == 2 for computation along rows. (The default is to choose
% dimension 1 if the number of rows is > 1.)
%
    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
