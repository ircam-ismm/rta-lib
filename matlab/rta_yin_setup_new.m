function [setup] = rta_yin_setup_new (yin_max_mins)
% function [setup] = rta_yin_setup_new (yin_max_mins)
%
% Create a new setup for the yin algorithm from the RTA library.
%
% <yin_max_mins> is optional (default is 128)
%
% <setup> is a pointer the the newly created setup
%
% Use rta_yin_setup_delete(setup) to release the memory of this setup
% when every rta_yin calculation is done.


printf(['This file is an help file which relies on a mex file of the same' ...
        ' name (but with a different extension).\n']);
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
