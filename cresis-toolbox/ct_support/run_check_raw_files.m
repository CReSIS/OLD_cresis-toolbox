% script run_check_raw_files
%
% Script for running check_raw_files.
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_records.m, create_records.m,
%   run_create_records_sync.m, check_records.m

%% User Setup
% =====================================================================
param_override = [];

% params = read_param_xls(ct_filename_param('rds_param_2019_Greenland_P3.xls'));
params = read_param_xls(ct_filename_param('snow_param_2019_Greenland_P3.xls'));

% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,'cmd.generic',1);
% params = ct_set_params(params,'cmd.generic',0,'cmd.notes','do not process');
% params = ct_set_params(params,'cmd.records',1,'day_seg','20190312_01');

dbstop if error;

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
  if param.cmd.generic
    check_raw_files(param,param_override);
  end
end
