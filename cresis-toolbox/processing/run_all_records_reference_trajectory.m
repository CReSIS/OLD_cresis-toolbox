% script run_all_records_reference_trajectory
%
% Script for running records_reference_trajectory on many seasons at once.
%
% Author: John Paden
%
% See also: records_reference_trajectory.m,
% records_reference_trajectory_load.m, run_records_reference_trajectory.m,
% run_all_records_reference_trajectory.m

run_all;

for param_fns_idx = 1:length(param_fns)
  param_fn = param_fns{param_fns_idx};
  param = regexp(param_fn,'(?<radar_name>\w+)_param_(?<season_name>\w+)','names');
  records_reference_trajectory(param);
end
