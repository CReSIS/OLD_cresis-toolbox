function labels = tracker_mcmc(data_struct,param)
%UNTITLED7 Summary of this function goes here
%   Detailed explanation goes here
% function mcmc_tracker_2D (params, options, data_struct)
%
% See also: run_tracker_2D.m
%
% Authors: Victor Berger

%% Automated Section
% =====================================================================
% Create param structure array
% =====================================================================

global gRadar;
% clear('param_override');

% Input checking
% if ~exist('params','var')
%   error('Use run_tracker_2D: A struct array of parameters must be passed in\n');
% end
% if exist('param_override','var')
%   param_override = merge_structs(gRadar, param_override);
% else
%   param_override = gRadar;
% end
if ~isfield(param.layer_tracker.track.mcmc,'alg')
  param.layer_tracker.track.mcmc.alg = 'MCMC';
end
if isfield(param.layer_tracker.track.mcmc, 'top_smooth')
  top_smooth = param.layer_tracker.track.mcmc.top_smooth;
else
  top_smooth = 1000;
end

if isfield(param.layer_tracker.track.mcmc, 'bottom_smooth')
  bottom_smooth = param.layer_tracker.track.mcmc.bottom_smooth;
else
  bottom_smooth = 1000;
end

if isfield(param.layer_tracker.track.mcmc, 'top_peak')
  top_peak = param.layer_tracker.track.mcmc.top_peak;
else
  top_peak = 0.5;
end

if isfield(param.layer_tracker.track.mcmc, 'bottom_peak')
  bottom_peak = param.layer_tracker.track.mcmc.bottom_peak;
else
  bottom_peak = 0.5;
end

if isfield(param.layer_tracker.track.mcmc, 'repulsion')
  repulsion = param.layer_tracker.track.mcmc.repulsion;
else
  repulsion = 10;
end

labels = {};

% Load frames file
% load(ct_filename_support(param,'','frames'));
% 
% if isempty(param.cmd.frms)
%   param.cmd.frms = 1:length(frames.frame_idxs);
% end

% Remove frames that do not exist from param.cmd.frms list
% [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
% if length(valid_frms) ~= length(param.cmd.frms)
%   bad_mask = ones(size(param.cmd.frms));
%   bad_mask(keep_idxs) = 0;
%   warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
%     param.cmd.frms(find(bad_mask,1)));
%   param.cmd.frms = valid_frms;
% end

%for frm = param.cmd.frms
  frm = param.layer_tracker.tracker.frm;
  mcmc_tic = tic;
  
  fprintf('\nMCMC: Running frame %s_%03d (%s)\n',param.day_seg, frm, datestr(now,'HH:MM:SS'));
  
%   try
%     big_matrix = data_struct.(sprintf('data_%s_%03d',param.day_seg,frm));
%     data = data_struct.(sprintf('data_%s_%03d',param.day_seg,frm));
%   catch ME
%     fprintf('\nProblem with frame %s_%03d, verify.\n' ,param.day_seg,frm);
%     continue;
%   end
  
  big_matrix = data_struct;
  data = data_struct;
  
  big_matrix.Data = lp(big_matrix.Data);
  big_matrix.Data = 255  * ((big_matrix.Data - min(big_matrix.Data(:))) ...
    ./ (max(big_matrix.Data(:) - min(big_matrix.Data(:)))));
  big_matrix.Data = uint8(repmat(big_matrix.Data, [1 1 3]));
  
  if strcmp(param.layer_tracker.track.mcmc.alg, 'MCMC')
    try
      tic
      [big_matrix.Labels, big_matrix.lower, big_matrix.upper] = ...
        tomo.RJ_MCMC(double(big_matrix.Data(:,:,1)));
      toc
    catch ME
      try
        mex -largeArrayDims RJ_MCMC.cpp
        tic
        [big_matrix.Labels, big_matrix.lower, big_matrix.upper] = ...
          tomo.RJ_MCMC(double(big_matrix.Data(:,:,1)));
        toc
      catch ME
        fprintf('\nProblem executing RJ_MCMC.cpp file. Verify.\n');
        keyboard
      end
    end
  else
    fprintf('\nUnrecognized algorithm (must be ''MCMC'')');
    keyboard
  end
  
  mcmc_toc = toc(mcmc_tic);
  
  if param.layer_tracker.track.debug
    %       figure; imshow(big_matrix.Data); hold on;
    %       plot(big_matrix.Labels(1, :), 'g'); plot(big_matrix.Labels(2, :), 'r');
    %       legend('Ice-surface', 'Ice-bottom');
    %       keyboard
    figure; imagesc(lp(data.Data)); colormap(1 - gray(256));hold on;
    plot(big_matrix.Labels(1, :), 'g'); plot(big_matrix.Labels(2, :), 'r');
    legend('Ice-surface', 'Ice-bottom');
    keyboard
  end
  
  if param.layer_tracker.track.ops_write
    warning('off');
    %% Set write options
    % Load labels into OPS using opsCopyLayers
    copy_param = [];
    copy_param.layer_source.existence_check = false;
    copy_param.layer_dest.existence_check = false;
    
    % Set the source
    copy_param.layer_source.source = 'custom';
    copy_param.layer_dest.source = param.layer_tracker.track.layer_dest_source;
    
    if strcmp(copy_param.layer_dest.source, 'layerdata')
      copy_param.layer_dest.layerdata_source = param.layer_tracker.track.layer_dest_layerdata_source;
      copy_param.layer_dest.echogram_source = param.layer_tracker.track.layer_dest_echogram_source;
    end
    
    copy_param.copy_method = 'fillgaps';
    
    copy_param.gaps_fill.method = 'preserve_gaps';
    copy_param.gaps_fill.method_args = [40 20];
    
    param = merge_structs(param,gRadar);
    
    %% Write surface layer
    % Interpolate from row number to TWTT
    big_matrix.TWTT = interp1(1:length(big_matrix.Time), big_matrix.Time, big_matrix.Labels(1,:));
    copy_param.layer_source.gps_time = big_matrix.GPS_time;
    copy_param.layer_source.twtt = big_matrix.TWTT;
    
    % Load labels into OPS using opsCopyLayers
    copy_param.layer_dest.name = param.layer_tracker.track.mcmc.lyrtop;
    fprintf('\nopsCopyLayers %s (%s) [TOP]', param.day_seg, datestr(now));
    opsCopyLayers(param,copy_param);
    fprintf('\n  Complete (%s)\n', datestr(now));
    
    %% Write bottom layer
    % Interpolate from row number to TWTT
    big_matrix.TWTT = interp1(1:length(big_matrix.Time), big_matrix.Time, big_matrix.Labels(2,:));
    copy_param.layer_source.gps_time = big_matrix.GPS_time;
    copy_param.layer_source.twtt = big_matrix.TWTT;
    
    % Load labels into OPS using opsCopyLayers
    copy_param.layer_dest.name = param.layer_tracker.track.mcmc.lyrbot;
    fprintf('\nopsCopyLayers %s (%s) [BOTTOM]', param.day_seg, datestr(now));
    opsCopyLayers(param,copy_param);
    
    fprintf('\n  Complete (%s)\n', datestr(now));
    warning('on');
  end
  
%   labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).top = ...
%     big_matrix.Labels(1,:);
%   labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).bot = ...
%     big_matrix.Labels(2,:);
%   labels.(sprintf('layer_%s_%03d',param.day_seg,frm)).toc = ...
%     mcmc_toc;
  labels = big_matrix.Labels;
%end
end