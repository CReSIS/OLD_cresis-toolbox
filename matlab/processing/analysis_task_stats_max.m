function [val_idx] = analysis_task_stats_max(param,cmd,data,start_bin,stop_bin)
% [val_idx] = analysis_task_stats_max(param,cmd,data,start_bin,stop_bin)
%
% Max function for use with analysis task command "statistics" that allows
% the maximum value and the index of the maximum value to be returned. This
% function uses nanmax.
%
% Author: John Paden
%
% See also analysis.m

val_idx = zeros(2,size(data,2));
for rline = 1:size(data,2)
  [val_idx(1,rline),val_idx(2,rline)] = nanmax(data(start_bin(rline):stop_bin(rline),rline));
end
