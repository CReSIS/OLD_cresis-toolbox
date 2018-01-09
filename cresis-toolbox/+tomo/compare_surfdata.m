function [frm_id,rmse_f,diff_f,med_f,min_f,max_f,total_diff] ...
  = compare_surfData(param,param_override)
% function [rmse_f,diff_f,med_f,min_f,max_f,total_diff] ...
%   = compare_surfData(param,param_override)
%
% Compares layer A (reference) and layer B
%   from a surfData file.
%
% param = struct with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% rmse_f: aggregation of all rmse results
% diff_f: 
%
% Example:
%  See run_compare_surfdata.m for how to run this function directly.
%
% Authors: Victor Berger, Mohanad Al-Ibadi, John Paden
%
% See also: tomo.run_compare_surfdata.m, tomo.compare_surfdata,
%   tomo.surfdata.m

%% General Setup
% =====================================================================

if ~isstruct(param)
  % Functional form
  param();
end
param = merge_structs(param, param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% Load frames file
load(ct_filename_support(param,param.records.frames_fn,'frames'));

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end
% Remove frames that do not exist from param.cmd.frms list
[valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
if length(valid_frms) ~= length(param.cmd.frms)
  bad_mask = ones(size(param.cmd.frms));
  bad_mask(keep_idxs) = 0;
  warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
    param.cmd.frms(find(bad_mask,1)));
  param.cmd.frms = valid_frms;
end

% LAYER A (REFERENCE):
surfdata_ref = param.compare.surfdata_ref;
layer_name_ref = param.compare.layer_name_ref;
% LAYER B (OTHER):
surfdata_other     = param.compare.surfdata_other;
layer_name_other   = param.compare.layer_name_other;
% BOTH lAYERS:
cutoffs            = param.compare.cutoffs;
DOA_trim           = param.compare.DOA_trim;

out_dir = param.compare.out_dir;

display_status_flag = param.compare.display_status_flag;
  
%% Process each frame
% =========================================================================
rmse_f  = zeros(size(param.cmd.frms));
diff_f = zeros(size(param.cmd.frms));
med_f  = zeros(size(param.cmd.frms));
max_df = zeros(size(param.cmd.frms));
min_df = zeros(size(param.cmd.frms));
frm_id = cell(size(param.cmd.frms));
total_diff = [];
h_figure = figure;
for frame_idx = 1:length(param.cmd.frms)
  frm = param.cmd.frms(frame_idx);
  
  frm_id{frame_idx} = sprintf('%s_%03d', param.day_seg, frm);
  
  fn_A = fullfile(ct_filename_out(param,surfdata_ref,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  fn_B = fullfile(ct_filename_out(param,surfdata_other,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  
  A = tomo.surfdata(fn_A);
  B = tomo.surfdata(fn_B);
  
  [rmse,mean_diff,median_diff,min_diff,max_diff,surf_diff] ...
    = A.compare(layer_name_ref, B,layer_name_other,DOA_trim);
  
  rmse_f(frame_idx) = rmse;
  diff_f(frame_idx) = mean_diff;
  med_f(frame_idx)  = median_diff;
  min_f(frame_idx) = min_diff;
  max_f(frame_idx) = max_diff;
  total_diff  = cat(2, total_diff, surf_diff);
  
  if display_status_flag
    fprintf('\n%s_%03d:', param.day_seg, frm);
    fprintf('\n Root mean squared error: %.4f', rmse);
    fprintf('\n Mean difference: %.4f', mean_diff);
    fprintf('\n Median difference: %.4f', median_diff);
    fprintf('\n Minimum difference: %.4f', min_diff);
    fprintf('\n Maximum difference: %.4f\n', max_diff);
    figure(h_figure); clf;
    imagesc(surf_diff); colorbar;
    warning('Run dbcont when ready.\n');
    keyboard;
  end
end

try
  delete(h_figure);
end
