% script run_update_img_combine
%
% Script for running update_img_combine
%
% Authors: John Paden
%
% See also: run_update_img_combine.m, update_img_combine.m

%% User Setup
% =====================================================================
params = read_param_xls(ct_filename_param('accum_param_2017_Greenland_P3.xls'),'');
% params = ct_set_params(params,'cmd.generic',0);

params = ct_set_params(params,'combine.img_comb_weights',[0 30],'day_seg','20170322_01');
params = ct_set_params(params,'combine.img_comb_weights',[0 30],'day_seg','20170322_02');
params = ct_set_params(params,'combine.img_comb_weights',[0 20],'day_seg','20170322_03');
params = ct_set_params(params,'combine.img_comb_weights',[0 20],'day_seg','20170322_04');
% params = ct_set_params(params,'combine.img_comb_weights_mode','auto','day_seg','20170322');
params = ct_set_params(params,'combine.img_comb_weights',[0 15],'day_seg','20170327_04');
params = ct_set_params(params,'combine.img_comb_weights',[0 15],'day_seg','20170328_01');
params = ct_set_params(params,'combine.img_comb_weights',[0 16],'day_seg','20170329_02');
params = ct_set_params(params,'combine.img_comb_weights',[0 18],'day_seg','20170330_01');
params = ct_set_params(params,'combine.img_comb_weights',[0 18],'day_seg','20170330_02');
params = ct_set_params(params,'combine.img_comb_weights',[0 18],'day_seg','20170330_03');
params = ct_set_params(params,'combine.img_comb_weights',[0 18],'day_seg','20170330_05');
params = ct_set_params(params,'combine.img_comb_weights',[0 15],'day_seg','20170331_01');

params = ct_set_params(params,'get_heights.img_comb_weights',[0 30],'day_seg','20170322_01');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 30],'day_seg','20170322_02');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 20],'day_seg','20170322_03');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 20],'day_seg','20170322_04');
% params = ct_set_params(params,'get_heights.img_comb_weights_mode','auto','day_seg','20170322');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 15],'day_seg','20170327_04');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 15],'day_seg','20170328_01');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 16],'day_seg','20170329_02');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 18],'day_seg','20170330_01');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 18],'day_seg','20170330_02');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 18],'day_seg','20170330_03');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 18],'day_seg','20170330_05');
params = ct_set_params(params,'get_heights.img_comb_weights',[0 15],'day_seg','20170331_01');

combine = [];
combine.img_comb_mult = 1.15; % <== OFTEN CHANGED (-inf default)
combine.img_comb_bins = 10; % <== OFTEN CHANGED (1 default)
out_path = 'CSARP_post/standard'; % <== OFTEN CHANGED
update_img_combine_param.mode = 'combine'; % <== OFTEN CHANGED (get_heights or combine)
update_img_combine_param.update_surf = true; % <== Usually false for RDS and true for non-RDS

combine.img_comb_layer_params = struct('name','surface','source','layerdata','layerdata_source','CSARP_post/layerData'); % <== OFTEN CHANGED

%% Automated Section
% =====================================================================

if strcmpi(update_img_combine_param.mode,'get_heights')
  combine.qlook.out_path = out_path;
else
  combine.out_path = out_path;
end

param_override = struct('combine',combine,'get_heights',combine, ...
  'update_img_combine',update_img_combine_param);

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
