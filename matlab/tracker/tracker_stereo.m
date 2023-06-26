function [labels,big_matrix] = tracker_stereo(data_struct,param)
global gRadar;

frm = param.layer_tracker.tracker.frm;
if ~isfield(param.layer_tracker.track.stereo,'alg')
  param.layer_tracker.track.stereo.alg = 'HMM';
end
if isfield(param.layer_tracker.track.stereo, 'top_smooth')
  top_smooth = param.layer_tracker.track.stereo.top_smooth;
else
  top_smooth = 1000;
end

if isfield(param.layer_tracker.track.stereo, 'bottom_smooth')
  bottom_smooth = param.layer_tracker.track.stereo.bottom_smooth;
else
  bottom_smooth = 1000;
end

if isfield(param.layer_tracker.track.stereo, 'top_peak')
  top_peak = param.layer_tracker.track.stereo.top_peak;
else
  top_peak = 0.5;
end

if isfield(param.layer_tracker.track.stereo, 'bottom_peak')
  bottom_peak = param.layer_tracker.track.stereo.bottom_peak;
else
  bottom_peak = 0.5;
end

if isfield(param.layer_tracker.track.stereo, 'repulsion')
  repulsion = param.layer_tracker.track.stereo.repulsion;
else
  repulsion = 10;
end

labels = {};

% Load frames file
load(ct_filename_support(param,'','frames'));

if isempty(param.cmd.frms)
  param.cmd.frms = 1:length(frames.frame_idxs);
end

% % Remove frames that do not exist from param.cmd.frms list
% [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
% if length(valid_frms) ~= length(param.cmd.frms)
%   bad_mask = ones(size(param.cmd.frms));
%   bad_mask(keep_idxs) = 0;
%   warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
%     param.cmd.frms(find(bad_mask,1)));
%   param.cmd.frms = valid_frms;
% end

%for frm = param.cmd.frms
  
  stereo_tic = tic;
  
  fprintf('\nStereo: Running frame %s_%03d (%s)\n',param.day_seg, frm, datestr(now,'HH:MM:SS'));
  
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
  
  if strcmp(param.layer_tracker.track.stereo.alg, 'HMM')
    opts = [top_smooth bottom_smooth top_peak bottom_peak repulsion];
    try
      tic
      [~, big_matrix.Labels] = tomo.stereo(1, double(big_matrix.Data(:,:,1)), opts);
      toc
    catch ME
      try
        mex -largeArrayDims stereo.cpp
        tic
        [~, big_matrix.Labels] = tomo.stereo(1, double(big_matrix.Data(:,:,1)), opts);
        toc
      catch ME
        fprintf('\nProblem executing stereo.cpp file. Verify.\n');
        keyboard
      end
    end
  else
    fprintf('\nUnrecognized algorithm (must be ''HMM'')');
    keyboard
  end
  
  stereo_toc = toc(stereo_tic);
  
  if param.layer_tracker.track.debug
    figure; imagesc(lp(data.Data));colormap(1 - gray(256)); hold on;
    plot(big_matrix.Labels(1, :), 'g'); plot(big_matrix.Labels(2, :), 'r');
    legend('Ice-surface', 'Ice-bottom');
    keyboard
    %     figure; imshow(big_matrix.Data); hold on;
    %     plot(big_matrix.Labels(1, :), 'g'); plot(big_matrix.Labels(2, :), 'r');
    %     legend('Ice-surface', 'Ice-bottom');
    %     keyboard
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
    copy_param.layer_dest.name = param.layer_tracker.track.stereo.lyrtop;
    fprintf('\nopsCopyLayers %s (%s) [TOP]', param.day_seg, datestr(now));
    opsCopyLayers(param,copy_param);
    fprintf('\n  Complete (%s)\n', datestr(now));
    
    %% Write bottom layer
    % Interpolate from row number to TWTT
    big_matrix.TWTT = interp1(1:length(big_matrix.Time), big_matrix.Time, big_matrix.Labels(2,:));
    copy_param.layer_source.gps_time = big_matrix.GPS_time;
    copy_param.layer_source.twtt = big_matrix.TWTT;
    
    % Load labels into OPS using opsCopyLayers
    copy_param.layer_dest.name = param.layer_tracker.track.stereo.lyrbot;
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
%     stereo_toc;
  
%end
labels = big_matrix.Labels;
end

