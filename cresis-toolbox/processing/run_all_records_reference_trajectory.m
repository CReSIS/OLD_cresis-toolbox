% script run_all_records_reference_trajectory
%
% Script for running records_reference_trajectory on many seasons at once.
%
% Author: John Paden
%
% See also: records_reference_trajectory.m,
% records_reference_trajectory_load.m, run_records_reference_trajectory.m,
% run_all_records_reference_trajectory.m

%% User Setup
% =====================================================================

param_override = [];

run_all;

%% Automated Section
% =====================================================================

for param_fns_idx = 1:length(param_fns)
  param_fn = ct_filename_param(param_fns{param_fns_idx});
  params = read_param_xls(param_fn);
  
  params = ct_set_params(params,'cmd.generic',1);
  params = ct_set_params(params,'cmd.generic',0,'cmd.notes','do not process');
  
  % Input checking
  global gRadar;
  if exist('param_override','var')
    param_override = merge_structs(gRadar,param_override);
  else
    param_override = gRadar;
  end
  
  % Process each of the segments
  for param_idx = 1:length(params)
    param = params(param_idx);
    if isfield(param.cmd,'generic') && ~iscell(param.cmd.generic) && ~ischar(param.cmd.generic) && param.cmd.generic
      records_reference_trajectory(param,param_override);
    end
  end
end
