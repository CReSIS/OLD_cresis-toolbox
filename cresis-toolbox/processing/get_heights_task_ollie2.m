function [success surfTimes] = get_heights_task_ollie2(steady_param_file_name)
% [success surfTimes] = get_heights_task_ollie2(steady_param_file_name)
%
% This function generates quick look outputs (default location: CSARP_qlook)
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_get_heights.m, get_heights.m,
%   get_heights_task.m

% =====================================================================
% General Setup
% =====================================================================

load(steady_param_file_name,'steady_param');
param=steady_param;

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

% =====================================================================
% Setup processing
% =====================================================================

% Get WGS84 ellipsoid parameters
physical_constants;

if ~isfield(param.records,'records_fn')
  param.records.records_fn = '';
end
if ~isfield(param.records,'frames_fn')
  param.records.frames_fn = '';
end

% Load frames file
load(ct_filename_support(param,param.records.frames_fn,'frames'));
% Load records file
records_fn = ct_filename_support(param,param.records.records_fn,'records');
records = load(records_fn);

global g_data;
g_data = [];

qlook_out_path = ct_filename_out(param, param.get_heights.qlook.out_path, 'CSARP_qlook');

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

if ~isfield(param.sched,'rerun_only') || isempty(param.sched.rerun_only)
  param.sched.rerun_only = false;
end

if ~isfield(param.get_heights,'ground_based')
  param.get_heights.ground_based = [];
end

if ~isfield(param.get_heights.qlook,'save_format') || isempty(param.get_heights.qlook.save_format)
  param.get_heights.qlook.save_format = '7.3';
end
save_format = sprintf('-v%s',param.get_heights.qlook.save_format);

if isfield(param.get_heights,'deconvolution') ...
    && ~isempty(param.get_heights.deconvolution) ...
    && param.get_heights.deconvolution == 3
  %% Get version information out of the deconvolution file
  out_fn_dir = ct_filename_out(param,'', 'CSARP_noise');
  out_segment_fn_dir = fileparts(out_fn_dir);
  out_segment_fn = fullfile(out_segment_fn_dir,sprintf('deconv_%s.mat', param.day_seg));
  spec = load(out_segment_fn,'param_collate');
  
  param.get_heights.deconvolution_sw_version = spec.param_collate.sw_version;
  param.get_heights.deconvolution_params = spec.param_collate.analysis.specular;
end

if isfield(param.get_heights,'coh_noise_method') ...
    && ~isempty(param.get_heights.coh_noise_method) ...
    && any(param.get_heights.coh_noise_method == [17 19])
  %% Get version information out of the coherent noise file
  
  cdf_fn_dir = fileparts(ct_filename_out(param,param.get_heights.coh_noise_arg{4}, ''));
  cdf_fn = fullfile(cdf_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
  
  tmp = netcdf_to_mat(cdf_fn,[],'^sw_version.*');
  param.get_heights.coh_noise_version = tmp.sw_version;
  tmp = netcdf_to_mat(cdf_fn,[],'^param_collate.*');
  param.get_heights.coh_noise_params = tmp.param_collate;
end

% =====================================================================
% Setup static inputs for get_heights_task
% =====================================================================
task_param = param;
task_param.load.imgs = param.get_heights.imgs;
if isempty(param.get_heights.imgs)
  error('No images specified in param.get_heights.imgs.');
end

%% Check img_comb
if length(param.get_heights.imgs) == 1 || isempty(param.get_heights.qlook.img_comb)
  num_imgs = 1;
else
  num_imgs = length(param.get_heights.imgs);
  if length(param.get_heights.qlook.img_comb) ~= 3*(num_imgs-1)
    warning('param.get_heights.qlook.img_comb not the right length. Since it is not empty, there should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
    keyboard
  end
end

% =====================================================================
%% Loop through all the frames: combine and surface track
[output_dir,radar_type] = ct_output_dir(param.radar_name);
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  % Check digits of proc_mode from frames file and make sure the user has
  % specified to process this frame type
  if ct_proc_frame(frames.proc_mode(frm),param.get_heights.frm_types)
    fprintf('get_heights combine frame %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  %% Output directory
  in_path = fullfile(qlook_out_path, sprintf('ql_data_%03d_01_01',frm));
  
  %% Concatenate blocks for each of the images
  for img = 1:length(param.get_heights.imgs)
    
    filenames = get_filenames(in_path,'','',sprintf('img_%02d.mat', img));
    
    if length(param.get_heights.imgs) == 1
      out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(qlook_out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    
    %% Concatenate all the ql block outputs for this image
    Time = [];
    Latitude = [];
    Longitude = [];
    Elevation = [];
    Roll = [];
    Pitch = [];
    Heading = [];
    GPS_time = [];
    Data = [];
    custom = [];
    if ~param.get_heights.combine_only
      for block_idx = 1:length(filenames)
        qlook_fn = filenames{block_idx};
        tmp = load(qlook_fn);
        time_vector_changed = false;
        if isempty(Time)
          Time = tmp.Time;
        elseif any(size(Time) ~= size(tmp.Time)) || any(Time ~= tmp.Time)
          % Determine the new time axis
          time_vector_changed = true;
          old_time = Time;
          dt = Time(2) - Time(1);
          start_time_diff = round((Time(1) - tmp.Time(1))/dt);
          end_time_diff = round((tmp.Time(end) - Time(end))/dt);
          if start_time_diff > 0
            Time = [Time(1)+dt*(-start_time_diff:-1)'; Time];
          end
          if end_time_diff > 0
            Time = [Time; Time(end)+dt*(1:end_time_diff)'];
          end
        end
        Latitude = [Latitude double(tmp.Latitude)];
        Longitude = [Longitude double(tmp.Longitude)];
        Elevation = [Elevation double(tmp.Elevation)];
        Roll = [Roll double(tmp.Roll)];
        Pitch = [Pitch double(tmp.Pitch)];
        Heading = [Heading double(tmp.Heading)];
        GPS_time = [GPS_time tmp.GPS_time];
        param_records = tmp.param_records;
        param_get_heights = tmp.param_get_heights;
        
        if isfield(tmp,'custom') && ~isempty(tmp.custom)
          % Custom fields are present, so concatenate them on the second dimension
          % - Used for deconvolution waveform index
          fields = fieldnames(tmp.custom);
          if block_idx == 1
            for field_idx = 1:length(fields)
              custom.(fields{field_idx}) = tmp.custom.(fields{field_idx});
            end
          else
            for field_idx = 1:length(fields)
              max_dim = length(size(tmp.custom.(fields{field_idx})));
              custom.(fields{field_idx}) = cat(max_dim,custom.(fields{field_idx}),tmp.custom.(fields{field_idx}));
            end
          end
          
        end
        
        if time_vector_changed
          if strcmpi(radar_type,'fmcw')
            Data = [interp1(old_time,Data,Time,'nearest',0) interp1(tmp.Time,tmp.Data,Time,'nearest',0)];
          else
            Data = [interp1(old_time,Data,Time,'linear',0) interp1(tmp.Time,tmp.Data,Time,'linear',0)];
          end
        else
          Data = [Data tmp.Data];
        end
      end
      %% Save output
      fprintf('  Writing output to %s\n', out_fn);
      Data = single(Data);
      if isempty(custom)
        save(save_format,out_fn,'Time','Latitude','Longitude', ...
          'Elevation','Roll','Pitch','Heading','GPS_time','Data', ...
          'param_get_heights','param_records');
      else
        save(save_format,out_fn,'Time','Latitude','Longitude', ...
          'Elevation','Roll','Pitch','Heading','GPS_time','Data', ...
          'param_get_heights','param_records','custom');
      end
    else
      fprintf('  Reading output %s\n', out_fn);
      load(out_fn);
    end
    
    %% Create temporary output for surface tracker
    if img == 1
      Time_Surface = Time;
      Data_Surface = Data;
    elseif ~isempty(param.get_heights.qlook.img_comb)
      %% Combine image with previous
      
      % Data_Surface,Time_Surface => already loaded data
      % Data, Time => new data to append
      % New_Time, New_Data => Combined result
      
      if Time(end) > Time_Surface(end)
        % Interpolate image N onto already loaded data (assumption is that image
        % N-1 always comes before image N)
        dt = Time_Surface(2)-Time_Surface(1);
        New_Time = (Time_Surface(1) : dt : Time(end)).';
        Data = interp1(Time,Data,New_Time,'linear',0);
        
        % Surface tracking image combine
        %  param.get_heights.qlook.img_comb(1): Not used at this step
        %  param.get_heights.qlook.img_comb(2): Not used at this step
        %  param.get_heights.qlook.img_comb(3): guard time which specifies how
        %    many seconds at the end of img1 will not be used... this is
        %    important because the last samples of img1 will have low signal
        %    power and blurred because they will only have captured a portion
        %    of the chirp energy (typically this will be set to something
        %    close to the pulse duration for img1)
        %  param.get_heights.qlook.img_comb(4-6, 7-9, etc.): same fields as above
        %    except between images 2 and 3, 3 and 4, etc.
        
        % Determine guard at end of image 1 that will not be used
        max_good_time = length(Time_Surface)*ones(1,size(Data_Surface,2));
        guard_bins = 1 + round(param.get_heights.qlook.img_comb((img-2)*3+3)/dt);
        
        % First row of img_bins indicates the start of the blend-region
        img_bins = max_good_time-guard_bins;
        
        % Second row of img_bins indicates the end of the blend-region
        img_bins(2,:) = img_bins(1,:) + 1;
        
        difference = 10^(-0/10);
        
        % Combine waveforms
        New_Data = zeros(size(Data),'single');
        for rline = 1:size(New_Data,2)
          trans_bins = img_bins(1,rline)+1:img_bins(2,rline);
          weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
          New_Data(:,rline) = [Data_Surface(1:img_bins(1,rline),rline); ...
            weights.*Data_Surface(trans_bins,rline) ...
            + difference*(1-weights).*Data(trans_bins,rline); ...
            difference*Data(img_bins(2,rline)+1:end,rline)];
        end
        Time_Surface = New_Time;
        Data_Surface = New_Data;
      end
      
    end
  end
  Time = Time_Surface;
  Data = Data_Surface;
  
  %% Run surface tracker
  surf = param.get_heights.surf;
  if isfield(param.get_heights.surf,'min_bin')
    % Convert time min_bin into range bins
    surf.min_bin = find(Time > param.get_heights.surf.min_bin, 1);
  end
  if isfield(param.get_heights.surf,'max_bin') && ~isempty(param.get_heights.surf.max_bin)
    % Convert time max_bin into range bins
    surf.max_bin = find(Time > param.get_heights.surf.max_bin, 1);
  end
  if isfield(param.get_heights.surf,'max_diff')
    % Convert time max_diff into range bins
    dt = Time(2) - Time(1);
    surf.max_diff = param.get_heights.surf.max_diff/dt;
  end
  
  if ~isfield(surf,'manual')
    surf.manual = false;
  end
  
  
  if isfield(surf,'feedthru')
    %% Optional feed through removal
    
    % Interpolate feed through power levels on to data time axis
    feedthru_threshold = interp1(surf.feedthru.time,surf.feedthru.power_dB,Time);
    feedthru_threshold = interp_finite(feedthru_threshold);
    
    % Remove all data not exceeding feed through threshold power
    for rline=1:size(Data,2)
      Data(:,rline) = Data(:,rline) .* (lp(Data(:,rline)) > feedthru_threshold);
    end
  end
  
  if ~isempty(param.get_heights.ground_based)
    % Hack for ground based radar, surface time is zero
    Surface = param.get_heights.ground_based * ones(1,size(Data,2));
    
  else
    if surf.manual
      [new_surface,pnt] = tracker_snake_simple(Data,surf);
      fprintf('  Press F1 for help\n');
      layer = tracker_snake_manual_gui(lp(Data),pnt);
      
    elseif strcmpi(surf.method,'threshold')
      new_surface = tracker_threshold(Data,surf);
    elseif strcmpi(surf.method,'max')
      new_surface = tracker_max(Data,surf);
    elseif strcmpi(surf.method,'snake')
      new_surface = tracker_snake_simple(Data,surf);
    else
      error('Not a supported surface tracking method.');
    end
    
    %% Apply optional median filter
    if isfield(surf,'medfilt') && ~isempty(surf.medfilt)
      new_surface = medfilt1(new_surface,surf.medfilt);
    end
    
    %% Convert from range bins to two way travel time
    Surface = interp1(1:length(Time), Time, new_surface);
    
    Surface = reshape(Surface, [1 length(Surface)]);
  end
  
  % Reset the "Data" variable in case it was modified during surface
  % tracking
  Data = Data_Surface;
  
  %% Combine images into a single image (also trim time<0 values)
  % Load each image and then combine with previous image
  for img = 1:num_imgs
    if length(param.get_heights.imgs) == 1
      out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(qlook_out_path, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    if img == 1
      load(out_fn);
      first_idx = find(Time <= 0,1,'last');
      if ~isempty(first_idx)
        Time = Time(first_idx:end);
        Data = Data(first_idx:end,:,:);
      end
    else
      append = load(out_fn,'Time','Data');
      %% Combine images
      % Data,Time => already loaded data
      % append.Data, append.Time => new data to append
      % New_Time, New_Data => Combined result
      
      % Interpolate image N onto already loaded data (assumption is that image
      % N-1 always comes before image N)
      dt = Time(2)-Time(1);
      New_Time = (Time(1) : dt : append.Time(end)).';
      append.Data = interp1(append.Time,append.Data,New_Time,'linear',0);
      
      % Surface tracking image combine
      %  param.get_heights.qlook.img_comb(1): time after surface return where
      %    combine will happen
      %  param.get_heights.qlook.img_comb(2): minimum time that combine will occur
      %  param.get_heights.qlook.img_comb(3): guard time which specifies how
      %    many seconds at the end of img1 will not be used... this is
      %    important because the last samples of img1 will have low signal
      %    power and blurred because they will only have captured a portion
      %    of the chirp energy (typically this will be set to something
      %    close to the pulse duration for img1)
      %  param.get_heights.qlook.img_comb(4-6, 7-9, etc.): same fields as above
      %    except between images 2 and 3, 3 and 4, etc.
      
      Surface = interp_finite(Surface,0);
      % First row of img_bins indicates the start of the blend-region
      img_bins = round(interp1(New_Time, 1:length(New_Time), ...
        max(Surface+param.get_heights.qlook.img_comb((img-2)*3+1),param.get_heights.qlook.img_comb((img-2)*3+2)), 'linear','extrap'));
      
      % Determine guard at end of image 1 that will not be used
      guard_bins = 1 + round(param.get_heights.qlook.img_comb((img-2)*3+3)/dt);
      
      % Check to make sure requested time is inside window and just
      % force the combination bin to occur at the second to last bin
      %   img_bins outside the img1 time window will be NaN due to interp1
      %   img_bins inside the img1 time window may still be larger than
      %     the guard allows
      max_good_time = length(Time)*ones(1,size(Data,2));
      invalid_rlines = find(isnan(img_bins) ...
        | img_bins > max_good_time-guard_bins);
      img_bins(invalid_rlines) = max_good_time(invalid_rlines)-guard_bins;
      
      % Second row of img_bins indicates the end of the blend-region
      img_bins(2,:) = img_bins(1,:) + 1;
      
      difference = 10^(-0/10);
      
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
  
  %% Save combined image output
  out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
    param.day_seg, frm));
  fprintf('  Writing output to %s\n', out_fn);
  Data = single(Data);
  if isempty(custom)
    save(save_format,out_fn,'Time','Latitude','Longitude', ...
      'Elevation','Roll','Pitch','Heading','GPS_time','Data','Surface', ...
      'param_get_heights','param_records');
  else
    save(save_format,out_fn,'Time','Latitude','Longitude', ...
      'Elevation','Roll','Pitch','Heading','GPS_time','Data','Surface', ...
      'param_get_heights','param_records','custom');
  end
  
  %% Remove the temporary block files
  if ~param.sched.rerun_only && exist(in_path,'dir')
    rmdir(in_path,'s');
  end
end

if param.get_heights.surf.en
  % Read the "Surface" variable from all the frames that were created
  % by this particular run of get_heights
  
  if isfield(records,'lat')
    num_recs = length(records.lat);
  end
  
  if ~isfield(records,'surface')
    records.surface = zeros(1,num_recs);
  end
  
  if length(records.surface) ~= num_recs
    warning('Debug catch incase records.surface is wrong length... hand-fix to be correct length');
    keyboard
  end
  
  for frm = param.cmd.frms
    if ~ct_proc_frame(frames.proc_mode(frm),param.get_heights.frm_types)
      continue;
    end
    out_fn = fullfile(qlook_out_path, sprintf('Data_%s_%03d.mat', ...
      param.day_seg, frm));
    load(out_fn,'GPS_time','Surface');
    
    if frm < length(frames.frame_idxs)
      recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1);
    else
      recs = frames.frame_idxs(frm):num_recs;
    end
    
    records.surface(recs) = interp1(GPS_time,Surface,records.gps_time(recs),'linear','extrap');
  end
  
  % Store surface information to the records file
  save(records_fn,'-APPEND','-struct','records','surface');
  create_records_aux_files(records_fn);
end

fprintf('Done %s\n', datestr(now));

success = true;

return;