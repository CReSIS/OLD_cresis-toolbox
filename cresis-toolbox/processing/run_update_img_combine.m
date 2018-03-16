% script run_update_img_combine
%
% Script for running update_img_combine
%
% Authors: John Paden
%
% See also: run_update_img_combine.m, update_img_combine.m

%% User Setup
% =====================================================================
params = read_param_xls(ct_filename_param('rds_param_2016_Antarctica_DC8.xls'),'');
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'day_seg','20161024_05');


mode = 'get_heights'; % <== OFTEN CHANGED (get_heights or combine)

update_img_combine_param.out_path = 'qlook';
update_img_combine_param.img_comb_mult = inf; % <== OFTEN CHANGED (inf default)
update_img_combine_param.img_comb_bins = 1; % <== OFTEN CHANGED (1 default)
update_img_combine_param.img_comb_layer_params = struct('name','surface','source','layerdata','layerdata_source','CSARP_post/layerData');% <== OFTEN CHANGED

%% Automated Section
% =====================================================================

param_override = [];
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
