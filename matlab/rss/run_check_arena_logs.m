% script run_check_arena_logs
%
% Script for running check_arena_logs
%
% Authors: John Paden
%
% See also: run_check_arena_logs.m, check_arena_logs.m, read_arena_logs.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'');

% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180927_20');

% dbstop if error;

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
ctrl_chain = {};
for param_idx = 1:length(params)
  param = params(param_idx);
  if param.cmd.generic
    check_arena_logs(param,param_override);
  end
end
