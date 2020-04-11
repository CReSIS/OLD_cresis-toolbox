% script run_gpr_find_bad_records
% 
% runs the script gpr_find_bad_records for all the supplied segments and
% frames.
%
% Author: Nick Holschuh, John Paden
%
% See also: run_gpr_find_bad_records.m, gpr_find_bad_records.m

param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'');
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20191231');

param_override.gpr_find_bad_records.manual_masking = 1;
param_override.gpr_find_bad_records.debug_plot = 0;
param_override.gpr_find_bad_records.bad_vel_threshold = 0.1; % <== OFTEN CHANGED (0.25 default)
param_override.gpr_find_bad_records.bad_heading_diff_threshold = Inf; % <== OFTEN CHANGED (1 default)

%% Automated Section
% =====================================================================

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
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  gpr_find_bad_records(param,param_override);
end

