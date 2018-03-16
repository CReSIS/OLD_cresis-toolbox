% script run_create_records_mcords2
%
% Script for running create_records_mcords2 (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_records_mcords2.m, create_records_mcords2.m,
%   create_records_mcords2_sync.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2018_Greenland_P3.xls'));

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');
params = ct_set_params(params,'cmd.records',0);
params = ct_set_params(params,'cmd.records',1,'day_seg','20180315_10');
% params = ct_set_params(params,'cmd.frms',[1]);

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
  if param.cmd.records
    create_records_mcords2(param,param_override);
  end
end
