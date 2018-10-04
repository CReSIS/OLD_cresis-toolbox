function success = preprocess_task(param)
% success = preprocess_task(param)
%
% Cluster enabled task called from preprocess.m.
%
% Example:
% Called from preprocess.m
%
% Author: John Paden
%
% See also: run_preprocess.m, preprocess.m, preprocess_task.m

if strcmpi(param.preprocess.daq.type,'arena')
  success = arena_packet_strip_task(param);
elseif strcmpi(param.preprocess.daq.type,'cresis')
  success = preprocess_task_cresis(param);
else
  error('Invalid param.preprocess.daq.type %s\n', param.preprocess.daq.type);
end


