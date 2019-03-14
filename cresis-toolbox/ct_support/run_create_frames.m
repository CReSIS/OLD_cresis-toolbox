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

% Syntax for running a specific segment by overriding parameter spreadsheet values
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'cmd.notes','do not process');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20180929_05');

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
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  % Create frames
  if param.records.frames.mode == 0
    fprintf('param.records.frames.mode == 0. Skipping.\n');
  elseif param.records.frames.mode == 1
    create_frames(param,param_override);
  elseif param.records.frames.mode == 2
    autogenerate_frames(param,param_override);
  end
end
