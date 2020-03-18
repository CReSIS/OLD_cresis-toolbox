% script run_collate_est_nz
%
% Runs collate_est_nz
% Classifies a day_segment into sets based on normalized cross correlation
% between blocks of coherent noise.
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
params = read_param_xls(param_fn,'',{'analysis'});


% % Enable a specific segment
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120414_01');

% % 
param_override.collate_est_nz.enable_visible_plot = 1; % Set 1/0; Default: 0
% param_override.collate_est_nz.enable_extralines_on_plot = 1; % Set 1/0; Default: 0
% param_override.collate_est_nz.enable_verbose = 1; % Set 1/0; Default: 0
% param_override.collate_est_nz.reuse_files = 1; % Set 1/0; Default: 0

% % Disconnect detect limit for max values.
% % Disconnect is detected if a max value for record is less than the limit = median - max( std * mul_factor , threshold);
switch 3
  case 1 % Set only mul_factor (and accept 10 dBm threshold)
    param_override.collate_est_nz.disconnect_mul_factor = 3; % Set a number; Default: 2
  case 2 % Set only threshold (and accept mul_factor = 2)
    param_override.collate_est_nz.disconnect_threshold = 5; % Set a number; Default: 10 dBm
  case 3 % Set HARD THRESHOLD (YOU MIGHT WANT TO USE THIS because mul_factor=0)
    param_override.collate_est_nz.disconnect_mul_factor = 0; % Set a number; Default: 2
    param_override.collate_est_nz.disconnect_threshold = 20; % Set a number; Default: 10 dBm
  case 4 % Set both as you like
    param_override.collate_est_nz.disconnect_mul_factor = 3; % Set a number; Default: 2
    param_override.collate_est_nz.disconnect_threshold = 15; % Set a number; Default: 10 dBm
end


% % To categorize regions, check consecutive diagonal elements of xcorr_norm matrix
% % If the difference in these normalized values is > xcorr tolerance,
% % then those blocks are considered as boundaries of regions.
param_override.collate_est_nz.xcorr_tol = 0.05; % Set a number(0:1); Default: 0.01

% % To recognize similar set of regions
param_override.collate_est_nz.xcorr_region_tol = 0.1; % Set a number(0:1); Default: 0.02

%% Automated Section
% =========================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% Process each of the segments
for param_idx =1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  collate_est_nz(param,param_override);
end