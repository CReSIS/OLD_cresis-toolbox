% script run_update_layerdata_format
%
% Script for running update_layerdata_format.m
%
% Authors: Jilu Li
%
% See also: update_layerdata_format.m, update_layerdata_format_task.m

%% User Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2012_Greenland_P3.xls'));
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_01');
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20120330_01|20120330_02|20120330_03');
params = ct_set_params(params,'cmd.frms',[1,2]);

param_override.frame_overlap_removal = true;
% param_override.cluster.type = 'debug';
% param_override.cluster.type = 'matlab';
param_override.cluster.type = 'torque';
param_override.cluster.rerun_only = false;

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
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  ctrl_chain{end+1} = update_layerdata_format(param,param_override);
end

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);