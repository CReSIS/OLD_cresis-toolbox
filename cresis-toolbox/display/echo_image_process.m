function echo_image_process(param,param_override)
% echo_image_process(param,param_override)
%
% FUNCTION DESCRIPTION
%
% param: struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override: parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_echo_image_process.m for how to run this function directly.
%  This function may be called from the run_master.m script using the
%  param spreadsheet and the cmd.generic column.
%
% Authors: FILL_IN_AS_NEEDED
%
% See also: run_echo_image_process.m,echo_image_process.m

%% General Setup
% =====================================================================

param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =====================================================================

% Load frames file
frames = frames_load(param);

param.cmd.frms = frames_param_cmd_frms(param,frames);
