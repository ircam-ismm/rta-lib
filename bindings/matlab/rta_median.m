function [value] = rta_median(input, dim)
% function [value] = rta_median(input, dim)
%
% Median of input.
% 
% Use dimension == 1 for computation along the columns and
% dimension == 2 for computation along rows. (The default is to choose
% dimension 1 if the number of rows is > 1.)
%
    
if(nargin < 2)
    if(size(input,1) == 1)
        dim = 2;
    else
        dim = 1;
    end
end

value = rta_selection(input, (1 + size(input, dim)) * 0.5, dim);


%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
