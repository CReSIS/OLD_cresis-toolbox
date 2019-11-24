% script run_collate_est_nz_records
%
% Runs collate_est_nz_records
% Updates old records using results from collate_est_nz_table
% Adds nyquist_zone_hw and records_mask to the settings fiels in records
%
% General order for processing:
% run_collate_est_nz, run_collate_est_nz_tables, run_collate_est_nz_records
%
% Authors: John Paden, Hara Madhav Talasila

%% USER SETTINGS
% =========================================================================

param_override = [];
param_sheet_name = 'snow_param_2012_Greenland_P3.xls';
param_fn = ct_filename_param(param_sheet_name);
params = read_param_xls(param_fn,'',{'collate_nz_est'});


% % % Enable a specific segment
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120326_01');

% % 
% param_override.collate_est_nz_table.enable_visible_plot = 1; % Set 1/0; Default: 0

%% Automated Section
% =========================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

params_en_idx = 0;
% Process each of the segments
for param_idx =1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  params_en_idx = params_en_idx + 1;
  params_en{params_en_idx} = param;
end

collate_est_nz_records(params_en,param_override);