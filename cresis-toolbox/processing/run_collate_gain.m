% script run_collate_gain
%
% Runs collate_gain
%
% Authors: John Paden, Hara Madhav Talasila

%% USER SETTINGS
% =========================================================================

param_override = [];
param_sheet_name = 'rds_param_2019_Greenland_P3.xls';
param_fn = ct_filename_param(param_sheet_name);
params = read_param_xls(param_fn,'',{'analysis_gain' 'analysis'});

param_override.analysis.enable_visible_plot = 1; % Set 1/0; Default: 0

% Directory of gain correction files
param_override.radar.ftg_dir = ...
  'X:\ct_data\ct_tmp\waveform\rds\2019_Greenland_P3\gain_correction_dir';

% To plot the fast time gain curves used in compensation
% param_override.analysis.ftg_plot_en = 0; % Set 1/0; Default: 1

% Not useful to see the raw plots
% param_override.analysis.raw_plot_en = 1; % Set 1/0; Default: 0

param_override.radar.wfs(1).TTL_start = 424;
param_override.radar.wfs(1).TTL_length = 239;
param_override.radar.wfs(2).TTL_start = 424;
param_override.radar.wfs(2).TTL_length = 350;
param_override.radar.wfs(3).TTL_start = 424;
param_override.radar.wfs(3).TTL_length = 739;

param_override.radar.wfs(1).record_start = 1328;
param_override.radar.wfs(1).record_stop = 5776;
param_override.radar.wfs(2).record_start = 1328;
param_override.radar.wfs(2).record_stop = 5776;
param_override.radar.wfs(3).record_start = 1328;
param_override.radar.wfs(3).record_stop = 5776;

param_override.radar.TTL_clock = 1e9/18;

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
  collate_gain;
end
