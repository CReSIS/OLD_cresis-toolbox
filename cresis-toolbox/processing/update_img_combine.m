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
%
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
  for img = 1:length(combine.imgs)
    
    if length(combine.imgs) == 1
      img_fn = fullfile(out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      img_fn = fullfile(out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    if img == 1
      load(img_fn);
      first_idx = find(Time <= 0,1,'last');
      if ~isempty(first_idx)
        Time = Time(first_idx:end);
        Data = Data(first_idx:end,:);
        if ~isempty(combine.img_comb_weights)
          Data = Data*10.^(combine.img_comb_weights(img)/10);
        end
      end
      
      if isfield(param.combine,'img_comb_layer_params') && ~isempty(param.combine.img_comb_layer_params)
        imb_comb_surf = interp1(layers.gps_time,layers.twtt,GPS_time);
      end
      if strcmpi(param.update_img_combine.mode,'get_heights')
        Surface = imb_comb_surf;
      elseif param.update_img_combine.update_surf
        Surface = imb_comb_surf;
      end
    else
      append = load(img_fn,'Time','Data');
      %% Combine images
      % Data,Time => already loaded data
      % append.Data, append.Time => new data to append
      % New_Time, New_Data => Combined result
      
      if ~isempty(combine.img_comb_weights)
        append.Data = append.Data.*10.^(combine.img_comb_weights(img)/10);
      end
      
      % Interpolate image N onto already loaded data (assumption is that image
      % N-1 always comes before image N)
      dt = Time(2)-Time(1);
      New_Time = (Time(1) : dt : append.Time(end)).';
      append.Data = interp1(append.Time,append.Data,New_Time,'linear',0);
      
      % Surface tracking image combine
      %  combine.img_comb(1): time after surface return where
      %    combine will happen
      %  combine.img_comb(2): minimum time that combine will occur
      %  combine.img_comb(3): guard time which specifies how
      %    many seconds at the end of img1 will not be used... this is
      %    important because the last samples of img1 will have low signal
      %    power and blurred because they will only have captured a portion
      %    of the chirp energy (typically this will be set to something
      %    close to the pulse duration for img1)
      %  combine.img_comb(4-6, 7-9, etc.): same fields as above
      %    except between images 2 and 3, 3 and 4, etc.
      
      imb_comb_surf = interp_finite(imb_comb_surf,0);
      % First row of img_bins indicates the start of the blend-region
      img_bins = round(interp1(New_Time, 1:length(New_Time), ...
        min(imb_comb_surf*combine.img_comb_mult, ...
        max(imb_comb_surf+combine.img_comb((img-2)*3+1),combine.img_comb((img-2)*3+2))), 'linear','extrap'));
      
      % Determine guard at end of image 1 that will not be used
      guard_bins = 1 + round(combine.img_comb((img-2)*3+3)/dt);
      
      blend_bins = combine.img_comb_bins;
      
      % Check to make sure requested time is inside window and just
      % force the combination bin to occur at the second to last bin
      %   img_bins outside the img1 time window will be NaN due to interp1
      %   img_bins inside the img1 time window may still be larger than
      %     the guard allows
      max_good_time = length(Time)*ones(1,size(Data,2));
      invalid_rlines = find(isnan(img_bins) ...
        | img_bins > max_good_time-guard_bins-blend_bins);
      img_bins(invalid_rlines) = max_good_time(invalid_rlines)-guard_bins-blend_bins;
      
      % Second row of img_bins indicates the end of the blend-region
      img_bins(2,:) = img_bins(1,:) + 1 + blend_bins;
      
      % Estimate difference
      if strcmpi(combine.img_comb_weights_mode,'auto')
        New_Data = zeros(size(append.Data),'single');
        difference = NaN*zeros(1,size(New_Data,2));
        for rline = 1:size(New_Data,2)
          trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
          if trans_bins <= size(append.Data,1)
            difference(rline) = mean(Data(trans_bins,rline) ./ append.Data(trans_bins,rline));
          end
        end
        difference = nanmedian(difference);
        fprintf('  Difference: %.1f\n', lp(difference));
        difference_report(frm) = difference;
      else
        difference = 10^(-0/10); % For debugging
      end
      
      % Combine images
      New_Data = zeros(size(append.Data),'single');
      for rline = 1:size(New_Data,2)
        trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
        weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
        if trans_bins <= size(append.Data,1)
          New_Data(:,rline) = [Data(1:img_bins(1,rline),rline); ...
            weights.*Data(trans_bins,rline) ...
            + difference*(1-weights).*append.Data(trans_bins,rline); ...
            difference*append.Data(img_bins(2,rline)+1:end,rline)];
        else
          New_Data(:,rline) = Data(1:size(New_Data,1),rline);
        end
      end
      Time = New_Time;
      Data = New_Data;
    end
  end
  
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
