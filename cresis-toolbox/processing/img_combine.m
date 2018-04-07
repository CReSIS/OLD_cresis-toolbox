function [Data, Time] = img_combine(param, param_mode, layers)
% [Data, Time] = img_combine(param, param_mode, layers)
%
% Blends and combines together individual echogram image files
% "Data_img_II_*" into a combined datafile "Data_*" according to the
% param structure.
%
% param: parameter structure usually loaded from parameter spreadsheet with
%   read_param_xls.m. The parameter structure will either look at the
%   get_heights or combine field depending on the param_mode parameter.
%  .day_seg: segment name to load
%  .load.frm: specifies the frame to load
%  .(param_mode): structure that defines the img_combine parameters
%   .out_path: output path to be updated
%   .img_comb: 3x(length(imgs)-1) vector which describes combination. Each
%     set of 3 numbers controls the merge for the corresponding two images:
%     (1): time after surface return where combine will happen
%     (2): minimum time that combine will occur
%     (3): guard time which specifies how many seconds at the end of the
%          first image that will not be used... this is important because
%          the last samples of img1 will have low signal power and blurred
%          because they will only have captured a portion of the chirp
%          energy (typically this will be set to something close to the
%          pulse duration for img1)
%     (4-6, 7-9, etc.): same fields as above except between images 2 and 3,
%          3 and 4, etc.
%   .img_comb_mult: if defined and nonempty, specifies a surface multiplier
%     for which the combine must happen after (used for dealing with
%     along-track surface clutter aliasing)
%   .img_comb_weights: specifies a weight in dB to apply to each image
%   .img_comb_weights_mode: can be empty or 'auto'. If auto, the intensity
%     difference between the two images in the transition region will be
%     forced to zero by adjusting the second image
%   .img_comb_bins: specifies a positive integer number of bins to use
%     when transitioning from one image to the next (default is 1).
%   .img_comb_layer_params: opsLoadLayers input parameters to load a layer
%     to use for the surface (if "layers" input argument is specified, then
%     this argument is ignored)
%   .trim_time: if true, negative time samples will be removed/trimmed
% param_mode: 'get_heights' or 'combine'
% layers: struct defining the two way travel time to the ice top for each
%   range line, must contain finite values
%  .gps_time: N element vector of GPS time's for the layer (ANSI C
%    standard, seconds since 1970).
%  .twtt: N element vector of two way travel time's for the layer (seconds)
%
% Authors: John Paden, Victor Berger
%
% See also: run_update_img_combine.m, update_img_combine.m

%% Input checks
% =========================================================================

if ~isfield(param.(param_mode), 'img_comb_weights')
  param.(param_mode).img_comb_weights = [];
end
if ~isfield(param.(param_mode), 'img_comb_mult')
  param.(param_mode).img_comb_mult = inf;
end
if ~isfield(param.(param_mode), 'img_comb_bins')
  param.(param_mode).img_comb_bins = 0;
end
if ~isfield(param.(param_mode), 'img_comb_weights_mode')
  param.(param_mode).img_comb_weights_mode = '';
end

%% Setup processing
% =========================================================================
% Output path
img_fn_dir = ct_filename_out(param,param.(param_mode).out_path,'');

if isempty(layers) && isfield(param.(param_mode),'img_comb_layer_params') && ~isempty(param.(param_mode).img_comb_layer_params)
  param_load_layers = param;
  param_load_layers.cmd.frms = param.load.frm;
  layers = opsLoadLayers(param_load_layers,param.(param_mode).img_comb_layer_params);
end

%% Combine images
% =========================================================================
% Check img_comb
if numel(param.(param_mode).imgs) == 1 || isempty(param.(param_mode).img_comb)
  num_imgs = 1;
else
  num_imgs = length(param.(param_mode).imgs);
  if length(param.(param_mode).img_comb) ~= 3*(num_imgs-1)
    error('param.%s.img_comb not the right length. Since it is not empty, there should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.', param_mode);
  end
end
for img = 1:num_imgs
  if length(param.(param_mode).imgs) == 1
    img_fn = fullfile(img_fn_dir, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, param.load.frm));
  else
    img_fn = fullfile(img_fn_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
      img, param.day_seg, param.load.frm));
  end
  
  %% Combine a pair of images: image "img" with image "img-1"
  % Data, Time => combined result
  % append.Data, append.Time => new data to append
  if img == 1
    load(img_fn,'Data','Time','GPS_time');
    if ~isempty(param.(param_mode).img_comb_weights)
      Data = Data*10.^(param.(param_mode).img_comb_weights(img)/10);
    end
    if isfield(param.(param_mode), 'trim_time') && param.(param_mode).trim_time
      first_idx = find(Time >= 0,1,'first');
      if ~isempty(first_idx)
        Time = Time(first_idx:end);
        Data = Data(first_idx:end,:);
      end
    end
    Surface = zeros(size(GPS_time));
    if ~isempty(layers)
      Surface = interp1(layers.gps_time,layers.twtt,GPS_time);
      Surface = interp_finite(Surface,0);
    end
    
    % Report auto difference to stdout
    if strcmpi(param.(param_mode).img_comb_weights_mode,'auto')
      fprintf('%.0f\t%.1f\n', img, lp(difference));
    end
    
  else
    append             = load(img_fn,'Time','Data');
    
    % Interpolate image N onto already loaded data (assumption is that image
    % N-1 always comes before image N)
    dt = Time(2)-Time(1);
    newTime = (Time(1) : dt : append.Time(end)).';
    append.Data = interp1(append.Time,append.Data,newTime,'linear',0);
    
    % Determine guard at end of image 1 that will not be used
    blend_bins = param.(param_mode).img_comb_bins;
    guard_bins = 1 + round(param.(param_mode).img_comb((img-2)*3+3)/dt);
    max_good_time = length(Time)*ones(1,size(Data,2));
    
    % First row of img_bins indicates the start of the blend-region
    if ~(Surface == 0)
      img_bins = round(interp1(newTime, 1:length(newTime), ...
        max(min(Surface * param.(param_mode).img_comb_mult, ...
        Surface + param.(param_mode).img_comb((img-2)*3 + 1)), ...
        param.(param_mode).img_comb((img-2)*3 + 2)), 'linear', 'extrap'));
    else
      img_bins = max_good_time-guard_bins;
    end
    
    % Check to make sure requested time is inside window and just
    % force the combination bin to occur at the second to last bin
    %   img_bins outside the img1 time window will be NaN due to interp1
    %   img_bins inside the img1 time window may still be larger than
    %     the guard allows
    invalid_rlines = find(isnan(img_bins) ...
      | img_bins > max_good_time-guard_bins-blend_bins);
    img_bins(invalid_rlines) = max_good_time(invalid_rlines)-guard_bins-blend_bins;
    
    % Second row of img_bins indicates the end of the blend-region
    img_bins(2,:) = img_bins(1,:) + 1 + blend_bins;
    
    % Estimate difference
    if strcmpi(param.(param_mode).img_comb_weights_mode,'auto')
      newData = zeros(size(append.Data),'single');
      difference = NaN*zeros(1,size(newData,2));
      for rline = 1:size(newData,2)
        trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
        if trans_bins <= size(append.Data,1)
          difference(rline) = mean(Data(trans_bins,rline) ./ append.Data(trans_bins,rline));
        end
      end
      difference = nanmedian(difference);
      % Report auto difference to stdout
      fprintf('\t%.1f', lp(difference));
    else
      difference = 10^(-0/10); % For debugging
    end
    
    % Combine images
    newData = zeros(size(append.Data),'single');
    for rline = 1:size(newData,2)
      trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
      weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
      if trans_bins <= size(append.Data,1)
        newData(:,rline) = [Data(1:img_bins(1,rline),rline); ...
          weights.*Data(trans_bins,rline) ...
          + difference*(1-weights).*append.Data(trans_bins,rline); ...
          difference*append.Data(img_bins(2,rline)+1:end,rline)];
      else
        newData(:,rline) = Data(1:size(newData,1),rline);
      end
    end
    
    % Update the outputs with the combined time axis and image data
    Time = newTime;
    Data = newData;
  end
  
  % Report auto difference to stdout
  if strcmpi(param.(param_mode).img_comb_weights_mode,'auto')
    fprintf('\n');
  end
end
