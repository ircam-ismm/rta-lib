function [output, last_sample] = ...
    rta_preemphasis (input, factor, previous_sample)
% function [output, last_sample] = ...
%     rta_preemphasis (input, factor, previous_sample)
%
% Apply preemphasis of <factor> on the <input> signal from the RTA
% library. With i in [0, input_size - 1],
%   <out_sample>[0] = <in_samples>[0] - <factor> * (*<previous_sample>)
%   <out_sample>[i] = <in_samples>[i] - <factor> * <in_samples>[i-1], i>0
% 0.97 is a common <factor> value for voice analysis.
% The <previous_sample> argument is optional (and is 0 by default).
% The <last_sample> is simply taken from <input>.

['This file is an help file which relies on a mex file of the same' ...
 ' name (but with a different extension).\n']
%
% 2008 (C) Ircam - Centre Pompidou
% Jean-Philippe.Lambert@ircam.fr
