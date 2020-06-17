% script run_records_bit_mask
% 
% runs the script records_bit_mask for all the supplied segments and
% frames.
%
% Author: Nick Holschuh, John Paden
%
% See also: run_records_bit_mask.m, records_bit_mask.m

param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'');
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20200108_02');
params = ct_set_params(params,'cmd.generic',1);
params = ct_set_params(params,'cmd.generic',0,'cmd.notes','do not process');

param_override.records_bit_mask.manual_masking = 0;
param_override.records_bit_mask.debug_plot = 0;
param_override.records_bit_mask.bad_vel_threshold = 0.1; % <== OFTEN CHANGED (0.25 default)
param_override.records_bit_mask.bad_heading_diff_threshold = Inf; % <== OFTEN CHANGED (1 default)

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
  records_bit_mask(param,param_override);
end

