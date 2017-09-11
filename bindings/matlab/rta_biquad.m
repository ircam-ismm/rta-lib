function [y, final_state] = rta_biquad (b, a, x, initial_state, dim, type)
% function [y, final_state] = rta_biquad (b, a, x, initial_state, dim, type)
%
% Compute a biquad filtering on 'x' input as
% y(n) = b(1) x(n) + b(2) x(n-1) + b(3) x(n-2) 
%                  - a(2) x(n-1) - a(3) x(n-2) 
%
% Note that a(1) is not taken into account (and supposed to be 1.)
% 'a' and 'b' must have 3 elements.
%
% If 'initial_state' is empty or not provided, it is supposed to be [0,0].
%
% Use dimension == 1 for computation along the columns and
% dimension == 2 for computation along the rows. (The default is to
% choose dimension 1 if the number of rows is > 1.)
% 
% if 'type' is 'df1', a direct form I is used (this is the default)
%  x           b0                                   y
%  --------+----->----->(  +  )-------->-----+------->
%          |            ^ ^ ^ ^              |
%          V   b1      / /   \ \        -a1  V
%        [x-1]--->----/ /     \ \------<---[y-1]
%          |           /       \             |
%          V   b2     /         \       -a2  V
%        [x-2]--->---/           \-----<---[y-2]
%
% if 'type' is 'df2t', a transposed direct form II is used:
%  x
%  --------+---------------+---------------+
%          |               |               |
%          |b2             |b1             |b0
%          V               V               V            y
%         (+)--->[z-1]--->(+)--->[z-1]--->(+)----+------->
%          ^               ^                     |
%          |-a2            |-a1                  |
%          |               |                     |
%          +---------------+---------------------+

    
['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
