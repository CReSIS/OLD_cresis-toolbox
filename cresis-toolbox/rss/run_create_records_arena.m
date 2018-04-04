% script run_create_records_arena
%
% Script for running create_records_arena (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_records_arena.m, create_records_arena.m,
%   create_records_arena_sync.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2016_Greenland_TOdtu.xls'));

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');
params = ct_set_params(params,'cmd.create_records',1);
% params = ct_set_params(params,'cmd.create_records',0,'day_seg','20161107_01');
params = ct_set_params(params,'cmd.create_records',0,'day_seg','20161107_02');
% params = ct_set_params(params,'cmd.create_records',0,'day_seg','20161108_02');
% params = ct_set_params(params,'cmd.records',1,'day_seg','20180315_10');
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
  if param.cmd.create_records
    create_records_arena(param,param_override);
  end
end
