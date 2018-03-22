% script run_collate_equalization
%
% Script for running collate_equalization which takes outputs from analysis
% surf and finds the time delay, phase, and amplitude differences between
% different wf_adc pairs.
%
% See "Receiver equalization" wiki page for details.
%
% Author: John Paden
%
% See also: collate_equalization

if 0
  %% 2015_Greenland_Polar6 Setup
  
  param_fn = ct_filename_param('rds_param_2017_Antarctica_Polar6.xls');
  params = read_param_xls(param_fn,'20160830_03',{'analysis_surf' 'analysis'});
  params.cmd.generic = 1;
  
  input_fn_dir = 'noise';
  img = 1;
  % zero_surf_bin_override = 11; % Normally not used except for internal layers
  
  debug_level = 1; % <-- TYPICAL SETTINGS ARE 1, 3, and 4

elseif 1
  %% 2018_Greenland_P3 Setup
  
  param_fn = ct_filename_param('rds_param_2018_Greenland_P3.xls');
  params = read_param_xls(param_fn,'20180315_09','analysis');
  params.cmd.generic = 1;
  
  input_fn_dir = 'analysis';
  img = 2;
  params.analysis.surf.ref_wf_adc = 5;
  params.analysis.surf.wf_adc_list = [1:4,6:16];
%   params.analysis.surf.ref_wf_adc = 2;
%   params.analysis.surf.wf_adc_list = [1:4];
%     params.analysis.surf.ref_wf_adc = 2;
%   params.analysis.surf.wf_adc_list = [13:16];
%   params.radar.wfs(8).rx_paths = [8 9 10 11 1 1 2 3 4 5 6 7 12 13 14 15];
  params.analysis.surf.motion_comp.en = true;
%   params.analysis.surf.ref_wf_adc = 3;
%   params.analysis.surf.wf_adc_list = [6:12];
  params.analysis.surf.rlines = [1:650];
  params.analysis.surf.rlines = [900:1450];
  params.analysis.surf.rlines = [1700:2400];
  params.analysis.surf.rlines = [2700:3200];
%   params.analysis.surf.rlines = [3500:4200];
%   params.analysis.surf.rlines = [];
  % zero_surf_bin_override = 11; % Normally not used except for internal layers
  
  debug_level = 1; % <-- TYPICAL SETTINGS ARE 1, 3, and 4

elseif 0
  %% 2017_Greenland_P3 Setup
  
  param_fn = ct_filename_param('rds_param_2014_Greenland_P3.xls');
  params = read_param_xls(param_fn,'20140401_03',{'analysis_surf' 'analysis'});
  params.cmd.generic = 1;
  
  input_fn_dir = 'noise';
  img = 1;
  % zero_surf_bin_override = 11; % Normally not used except for internal layers
  
  debug_level = 1; % <-- TYPICAL SETTINGS ARE 1, 3, and 4

elseif 0
  %% 2016_Antarctica_DC8 Setup
  
  param_fn = ct_filename_param('rds_param_2016_Antarctica_DC8.xls');
  params = read_param_xls(param_fn,'20161004_08',{'analysis_surf' 'analysis'});
  params.cmd.generic = 1;
  
  input_fn_dir = 'noise';
  img = 1;
  % zero_surf_bin_override = 11; % Normally not used except for internal layers
  
  debug_level = 1; % <-- TYPICAL SETTINGS ARE 1, 3, and 4
end

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  fprintf('=====================================================================\n');
  fprintf('%s: %s (%s)\n', 'collate_equalization', param.day_seg, datestr(now,'HH:MM:SS'));
  fprintf('=====================================================================\n');
  
  %param.analysis.surf.rlines = [1000:2000]; % <-- MAY WANT TO HAND MODIFY
  % param.analysis.surf.motion_comp.en = true; % <-- MAY WANT TO HAND MODIFY
  % param.analysis.surf.wf_adc_list = [9:16]; % <-- MAY WANT TO HAND MODIFY
  % param.analysis.surf.ref_wf_adc = 9; % <-- MAY WANT TO HAND MODIFY
  % param.analysis.surf.retrack.en = false; % <-- MAY WANT TO HAND MODIFY
  collate_equalization;

  fprintf('  Done (%s)\n', datestr(now));
end
