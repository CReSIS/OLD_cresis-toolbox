% script run_create_frames
%
% Script for running create_frames (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_create_frames.m, create_frames.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('accum_param_2018_Antarctica_TObas.xls'),'');

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');
% params = ct_set_params(params,'cmd.create_frames',0);
% params = ct_set_params(params,'cmd.create_frames',1,'day_seg','20161101_0[12345]');

param_override.ct_file_lock = 1;

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
    if param.records.frames.mode == 0
      create_frames(param,param_override);
      fprintf('Type dbcont to continue when you are done creating frames for this segment.\n');
      keyboard;
    else
      autogenerate_frames(param,param_override);
    end
  end
end
