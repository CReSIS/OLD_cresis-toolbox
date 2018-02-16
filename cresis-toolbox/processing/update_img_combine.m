function update_img_combine(param,param_override)
% update_img_combine(param,param_override)
%
% update_img_combine: Takes img_II echogram data files and combines them using
% difference img_comb* parameters into a combined data file. This function
% is primarily for redoing combining that has been done during get_heights
% and combine_wf_chan so that different img_comb* parameters can be used.
%
% param: struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override: parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Example:
%  See run_update_img_combine.m for how to run this function directly.
%  This function may be called from the run_master.m script using the
%  param spreadsheet and the cmd.generic column.
%
% Authors: John Paden
% See also: run_update_img_combine.m, update_img_combine.m

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

%% Setup processing
% =====================================================================

if ~isfield(param.combine,'img_comb_mult') || isempty(param.combine.img_comb_mult)
  param.combine.img_comb_mult = inf;
end

if ~isfield(param.combine,'img_comb_weights') || isempty(param.combine.img_comb_weights)
  param.combine.img_comb_weights = [];
end

if ~isfield(param.combine,'img_comb_weights_mode') || isempty(param.combine.img_comb_weights_mode)
  param.combine.img_comb_weights_mode = '';
end

if ~isfield(param.combine,'img_comb_bins') || isempty(param.combine.img_comb_bins)
  param.combine.img_comb_bins = 1;
end

if ~isfield(param.get_heights,'img_comb_mult') || isempty(param.get_heights.img_comb_mult)
  param.get_heights.img_comb_mult = inf;
end

if ~isfield(param.get_heights,'img_comb_weights') || isempty(param.get_heights.img_comb_weights)
  param.get_heights.img_comb_weights = [];
end

if ~isfield(param.get_heights,'img_comb_weights_mode') || isempty(param.get_heights.img_comb_weights_mode)
  param.get_heights.img_comb_weights_mode = '';
end

if ~isfield(param.get_heights,'img_comb_bins') || isempty(param.get_heights.img_comb_bins)
  param.get_heights.img_comb_bins = 1;
end
if ~isfield(param.update_img_combine,'update_surf') || isempty(param.update_img_combine.update_surf)
  param.update_img_combine.update_surf = false;
end

% Load frames file
load(ct_filename_support(param,'','frames'));

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

%% Load surface information
if isfield(param.combine,'img_comb_layer_params') && ~isempty(param.combine.img_comb_layer_params)
  param_load_layers = param;
  param_load_layers.cmd.frms = 1:length(frames.frame_idxs);
  
  layers = opsLoadLayers(param_load_layers,param.combine.img_comb_layer_params);
end

if strcmpi(param.update_img_combine.mode,'get_heights')
  out_path = ct_filename_out(param,param.get_heights.qlook.out_path,'');
else
  out_path = ct_filename_out(param,param.combine.out_path,'');
end

difference_report = nan(size(frames.frame_idxs));
for frm = param.cmd.frms
  
  out_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
    param.day_seg, frm));
  fprintf('Combine %s (%s)\n', out_fn, datestr(now));
  
  if strcmpi(param.update_img_combine.mode,'get_heights')
    load(out_fn,'param_get_heights','Surface');
    param.get_heights.imgs = param_get_heights.combine.imgs;
    combine = param.get_heights;
    combine.img_comb = combine.qlook.img_comb;
  elseif strcmpi(param.update_img_combine.mode,'combine')
    load(out_fn,'param_combine','Surface');
    param.combine.imgs = param_combine.combine.imgs;
    combine = param.combine;
  else
    error('Invalid param.update_img_combine.mode %s', param.update_img_combine.mode);
  end
  
  if isempty(combine.img_comb)
    % No image combining is required
    continue;
  end
  
  if length(combine.img_comb) ~= 3*(length(combine.imgs)-1)
    if strcmpi(param.update_img_combine.mode,'get_heights')
      warning('param.get_heights.qlook.img_comb not the right length. There should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
    else
      warning('param.combine.img_comb not the right length. There should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
    end
    keyboard
  end
  
  %% Load each image and then combine with previous image (also trim time<0 values)
  % Call img_combine
  combine.out_path = out_path;
  combine.frm      = frm;
  [Data, Time]     = img_combine(param, combine, layers);
  
  %% Save output
  if strcmpi(param.update_img_combine.mode,'combine')
    % combine_wf_chan file
    save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
      'Elevation','GPS_time','Data','Surface','Bottom', ...
      'param_combine','param_records','param_csarp', ...
      'Roll', 'Pitch', 'Heading');
  elseif strcmpi(param.update_img_combine.mode,'get_heights')
    % get_heights file
    save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
      'Elevation','GPS_time','Data','Surface', ...
      'param_get_heights','param_records', ...
      'Roll', 'Pitch', 'Heading');
  end
  
end

if strcmpi(combine.img_comb_weights_mode,'auto')
  fprintf('  Difference: %.1f\n', lp(nanmean(difference_report)));
  fprintf('  Difference Min: %.1f\n', lp(nanmin(difference_report)));
  fprintf('  Difference Median: %.1f\n', lp(nanmedian(difference_report)));
end
