% script test_update_img_combine
%
% Script for testing update_img_combine
%
% Authors: John Paden
%
% See also: run_update_img_combine.m, update_img_combine.m

%% Test Setup
% =====================================================================
param_override = [];

params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'20161024_05');

params.cmd.generic = 1;
params.cmd.frms = 1;

test_run = 2;
if test_run == 1
  mode = 'get_heights';
  update_img_combine_param.out_path = 'qlook';
elseif test_run == 2
  mode = 'combine';
  update_img_combine_param.out_path = 'standard';
end

update_img_combine_param.img_comb_mult = inf;
update_img_combine_param.img_comb_bins = 1;
update_img_combine_param.img_comb_layer_params = [];

%% Automated Section
% =====================================================================

param_override = [];
param_override.out_path = fullfile(gRadar.out_path,'tests');
param_override.(mode) = update_img_combine_param;
param_override.update_img_combine.mode = mode;

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
  update_img_combine(param,param_override);
end
