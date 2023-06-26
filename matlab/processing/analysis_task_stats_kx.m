function [data_kx] = analysis_task_stats_kx(param,cmd,data,start_bin,stop_bin)
% [data_kx] = analysis_task_stats_kx(param,cmd,data,start_bin,stop_bin)
%
% Doppler/along-track wavenumber function for use with analysis task
% command "statistics" that allows the Along-track kx wavenumber/Doppler
% FFT to be averaged.
%
% Author: John Paden
%
% See also analysis.m

data_kx = zeros(stop_bin(1)-start_bin(1)+1, cmd.kx);
breaks = 1:cmd.kx:size(data,2);
if size(data,2) - breaks(end) < cmd.kx
  % Remove incomplete end block
  breaks = breaks(end-1);
end
for break_idx = 1:length(breaks)
  data_kx = data_kx + abs(fft(data(start_bin(1):stop_bin(1), (break_idx-1)*cmd.kx + (1:cmd.kx)),[],2)).^2;
end
data_kx = data_kx / length(breaks);
