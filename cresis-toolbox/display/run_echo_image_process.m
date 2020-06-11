% script run_echo_image_process
%
% Script for running echo_image_process
%
% Authors: FILL_IN_AS_NEEDED
%
% See also: run_echo_image_process.m, echo_image_process.m

%% User Setup
% =====================================================================
param_override = [];

% params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'));
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'));
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20140313_08');
params = ct_set_params(params,'cmd.frms',[]); % Specify specific frames (or leave empty/undefined to do all frames)


% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO.xls'),'20091228_01');
%params = ct_set_params(params,'cmd.generic',0);
%params = ct_set_params(params,'cmd.generic',1,'day_seg','20091224_01|20091228|20091229|20100101_02|20100103|20100104_01|20100112');
%params = ct_set_params(params,'cmd.frms',[]);

% Set override parameters using ct_set_params, by setting param_override, or by setting FUNCTION_params variable which eventually is copied to param_override.FUNCTION).
% These two examples for setting input parameters are equivalent.
% Example 1
% params = ct_set_params(params,'FUNCTION.param1','first parameters');
% params = ct_set_params(params,'FUNCTION.param2',[2]);
% Example 2
% param_override.FUNCTION.param1 = 'first parameters';
% param_override.FUNCTION.param2 = [2];

% Default rds profile will be called "rds" (default is equal to sys from ct_output_dir(radar_name)...
param_override.echo_image_process.profile = 'AWI_rds';
param_override.echo_image_process.in_path = 'standard'; % This is the default param.echo_image_process.out_path = 'standard_image_process'; % This will be the default out path (append image_process to in_path)

%param_override.sched.type = 'no scheduler'; % Example to override default cluster settings

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
  echo_image_process(param,param_override);
end

% Post process code (only include if necessary)