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

if strcmpi(param.config.daq_type,'arena')
  success = preprocess_task_arena(param);
elseif strcmpi(param.config.daq_type,'cresis')
  success = preprocess_task_cresis(param);
else
  error('Invalid param.config.daq_type %s\n', param.config.daq_type);
end


