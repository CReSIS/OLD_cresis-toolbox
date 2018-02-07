% script run_get_heights
%
% Script for running get_heights (usually just used for debugging).
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_get_heights.m, get_heights.m,
%   get_heights_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'));

% Syntax for running a specific segment and frame by overriding parameter spreadsheet values
%params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'),'20161017_02');
params = ct_set_params(params,'cmd.get_heights',0);
params = ct_set_params(params,'cmd.get_heights',1,'day_seg','20161017_02');
params = ct_set_params(params,'cmd.frms',[50]);
params = ct_set_params(params,'cmd.get_heights',1,'day_seg','20161017_01');
params = ct_set_params(params,'cmd.frms',[1]);



% 
% params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'));
% % 
% % % Syntax for running a specific segment and frame by overriding parameter spreadsheet values
% % %params = read_param_xls(ct_filename_param('snow_param_2016_Antarctica_DC8.xls'),'20161017_02');
% params = ct_set_params(params,'cmd.get_heights',0);
% params = ct_set_params(params,'cmd.get_heights',1,'day_seg','20161024_05');
% params = ct_set_params(params,'cmd.frms',[]);



dbstop if error;
param_override.cluster.type = 'torque';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
%param_override.cluster.rerun_only = true;
param_override.cluster.desired_time_per_job  = 70*60;

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
  if param.cmd.get_heights
    ctrl_chain{end+1} = get_heights(param,param_override);
  end
end

cluster_print_chain(ctrl_chain);

[tmp tmp_name] = fileparts(tempname);
chain_fn = fullfile(gRadar.cluster.data_location,sprintf('chain_%s.mat',tmp_name));
fprintf('Saving chain: %s\n', chain_fn);
save(chain_fn,'ctrl_chain');

load(chain_fn);
ctrl_chain = cluster_run(ctrl_chain);
