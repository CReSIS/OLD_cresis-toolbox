function success = array_combine_task(param)
% success = array_combine_task(param)
%
% This script combines the chunks and images.
%
% The assumption is that the directories in the input_path are named
% using the following convention:
%   METHOD_FFF
% where
%   METHOD: Name of method
%   FFF: Zero-padded Frame ID
% Examples:
%   standard_001: Contains chunk files from standard array method, frame 1
%
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_array.m, array.m, load_sar_data.m,
% array_proc.m, array_task.m, array_combine_task.m

%% Setup processing
% =====================================================================

% Create output directory path
array_out_dir = ct_filename_out(param, param.array.out_path);
if ~exist(array_out_dir,'dir')
  mkdir(array_out_dir);
end

load(ct_filename_support(param,'','frames')); % Load "frames" variable

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
% Apply presumming
if param.sar.presums > 1
  records.lat = fir_dec(records.lat,param.sar.presums);
  records.lon = fir_dec(records.lon,param.sar.presums);
  records.elev = fir_dec(records.elev,param.sar.presums);
  records.roll = fir_dec(records.roll,param.sar.presums);
  records.pitch = fir_dec(records.pitch,param.sar.presums);
  records.heading = fir_dec(records.heading,param.sar.presums);
  records.gps_time = fir_dec(records.gps_time,param.sar.presums);
  records.surface = fir_dec(records.surface,param.sar.presums);
end
along_track_approx = geodetic_to_along_track(records.lat,records.lon,records.elev);

%% Combine chunks into each frame
% =====================================================================
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  if ct_proc_frame(frames.proc_mode(frm),param.array.frm_types)
    fprintf('array_combine %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  % Temporary output directory for uncombined array processed images
  array_fn_dir = fullfile(ct_filename_out(param, param.array.out_path, 'array_tmp'), ...
    sprintf('array_%03d', frm));

  % Current frame goes from the start record specified in the frames file
  % to the record just before the start record of the next frame.  For
  % the last frame, the stop record is just the last record in the segment.
  start_rec = ceil(frames.frame_idxs(frm)/param.sar.presums);
  if frm < length(frames.frame_idxs)
    stop_rec = ceil((frames.frame_idxs(frm+1)-1)/param.sar.presums);
  else
    stop_rec = length(records.gps_time);
  end
  
  % Determine length of the frame
  frm_dist = along_track_approx(stop_rec) - along_track_approx(start_rec);
  
  % Determine number of chunks for this frame
  num_chunks = round(frm_dist / param.array.chunk_len);

  %% Loop through all the images
  for img = 1:length(param.array.imgs)
    %% Loop through all the chunks and combine into an image
    Latitude = [];
    Longitude = [];
    Elevation = [];
    Roll = [];
    Pitch = [];
    Heading = [];
    GPS_time = [];
    Surface = [];
    Bottom = [];
    Data = [];
    Tomo = [];
    for chunk_idx = 1:num_chunks
      array_fn = fullfile(array_fn_dir, sprintf('img_%02d_chk_%03d.mat', img, chunk_idx));
      tmp = load(array_fn);
      
      if length(tmp.Time) == 1
        % Force length==1 in fast time to just be empty to simplify data
        % handling later. The idea is that a length 1 range line is useless
        % anyway so we might as well simplify things by making it zero length.
        tmp.Time = [];
        tmp.Data = tmp.Data([],:);
      end
      time_vector_changed = false;
      if chunk_idx == 1
        Time = tmp.Time;
      elseif any(size(Time) ~= size(tmp.Time)) || any(Time ~= tmp.Time)
        % Determine the new time axis
        %   Note that even though time axis is aligned with multiples of
        %   dt, there will be rounding errors which need to be dealt with
        %   here.
        time_vector_changed = true;
        if isempty(Time)
          Time = tmp.Time;
          start_time_diff = size(Time,1);
          end_time_diff = 0;
        elseif isempty(tmp.Time)
          start_time_diff = -size(Time,1);
          end_time_diff = 0;
        else
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
      end
      Latitude = [Latitude double(tmp.Latitude)];
      Longitude = [Longitude double(tmp.Longitude)];
      Elevation = [Elevation double(tmp.Elevation)];
      Roll = [Roll double(tmp.Roll)];
      Pitch = [Pitch double(tmp.Pitch)];
      Heading = [Heading double(tmp.Heading)];
      GPS_time = [GPS_time tmp.GPS_time];
      Surface = [Surface double(tmp.Surface)];
      Bottom = [Bottom double(tmp.Bottom)];
      if time_vector_changed
        Data = [[zeros(start_time_diff,size(Data,2)); Data; zeros(end_time_diff,size(Data,2))], ...
          [zeros(-start_time_diff,size(tmp.Data,2)); tmp.Data; zeros(-end_time_diff,size(tmp.Data,2))]];
      else
        Data = [Data tmp.Data];
      end
      param_records = tmp.param_records;
      param_sar = tmp.param_sar;
      if chunk_idx == 1
        param_array = tmp.param_array;
        param_array.array_proc.fcs{1}{1}.x = tmp.param_array.array_proc.fcs{1}{1}.x(:,tmp.param_array.array_proc.lines);
        param_array.array_proc.fcs{1}{1}.y = tmp.param_array.array_proc.fcs{1}{1}.y(:,tmp.param_array.array_proc.lines);
        param_array.array_proc.fcs{1}{1}.z = tmp.param_array.array_proc.fcs{1}{1}.z(:,tmp.param_array.array_proc.lines);
        param_array.array_proc.fcs{1}{1}.origin = tmp.param_array.array_proc.fcs{1}{1}.origin(:,tmp.param_array.array_proc.lines);
      else
        % Concatenate the fcs field
        param_array.array_proc.fcs{1}{1}.x = [param_array.array_proc.fcs{1}{1}.x tmp.param_array.array_proc.fcs{1}{1}.x(:,tmp.param_array.array_proc.lines)];
        param_array.array_proc.fcs{1}{1}.y = [param_array.array_proc.fcs{1}{1}.y tmp.param_array.array_proc.fcs{1}{1}.y(:,tmp.param_array.array_proc.lines)];
        param_array.array_proc.fcs{1}{1}.z = [param_array.array_proc.fcs{1}{1}.z tmp.param_array.array_proc.fcs{1}{1}.z(:,tmp.param_array.array_proc.lines)];
        param_array.array_proc.fcs{1}{1}.origin = [param_array.array_proc.fcs{1}{1}.origin tmp.param_array.array_proc.fcs{1}{1}.origin(:,tmp.param_array.array_proc.lines)];
      end
      if param.array.tomo_en
        %         3D-surface is present so concatenate it too
        %         Tomo = cat(3,Tomo,tmp.Tomo);
        %         Concatenate all the fields under struct Tomo: valR, bins, val, freq
        %         and img.
        if time_vector_changed
          error('Support for chunks with different time vectors not added yet.');
        end
        fields = fieldnames(tmp.Tomo);
        if chunk_idx == 1
          for field_idx = 1:length(fields)
            Tomo.(fields{field_idx}) = tmp.Tomo.(fields{field_idx});
          end
        else
          for field_idx = 1:length(fields)
            max_dim = length(size(tmp.Tomo.(fields{field_idx})));
            Tomo.(fields{field_idx}) = cat(max_dim,Tomo.(fields{field_idx}),tmp.Tomo.(fields{field_idx}));
          end
        end
        
      end
    end
    
    % =====================================================================
    % Save output
    if length(param.array.imgs) == 1 && ~param.array.tomo_en
      out_fn = fullfile(array_out_dir, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(array_out_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    fprintf('  Writing output to %s\n', out_fn);
    if param.ct_file_lock
      file_version = '1L';
    else
      file_version = '1';
    end
    Data = single(Data);
    if ~param.array.tomo_en
      % Do not save 3D surface
      save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
        'Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_array','param_records','param_sar', ...
        'Roll', 'Pitch', 'Heading','file_version');
    else
      % Save 3D surface
      save('-v7.3',out_fn,'Tomo','Time','Latitude', ...
        'Longitude','Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_array','param_records','param_sar', ...
        'Roll', 'Pitch', 'Heading','file_version');
    end
  end
  
  %% Delete temporary files now that all combined files are created
  if 0 % HACK: NEED TO REMOVE THE "if 0"
    for img = 1:length(param.array.imgs)
      % Determine where breaks in processing blocks are going to occur
      for chunk_idx = 1:num_chunks
        array_fn = fullfile(array_fn_dir, sprintf('img_%02d_chk_%03d.mat', img, chunk_idx));
        delete(array_fn);
      end
    end
    % Attempt to remove METHOD_FFF directory since it is no longer
    % needed.
    try
      rmdir(array_fn_dir);
    end
  end
  
  %% Combine images
  [output_dir,radar_type] = ct_output_dir(param.radar_name);
  if ~isempty(param.array.img_comb) || param.array.tomo_en || (length(param.array.imgs) == 1 && ~strcmpi(radar_type,'deramp')))
    % Combine images into a single image and/or trim invalid times with
    % img_comb_trim
    img_combine_param = param;
    img_combine_param.load.frm = frm;
    surf_layer.gps_time = GPS_time;
    surf_layer.twtt = Surface;
    [Data, Time] = img_combine(img_combine_param, 'array', surf_layer);
    
    %% Save combined image output
    if length(param.array.imgs) == 1 || ~isempty(param.array.img_comb)
      % A combined file should be created
      out_fn = fullfile(array_out_dir, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
      fprintf('  Writing output to %s\n', out_fn);
      % Note that image combining here never includes "Tomo" variable. Use
      % tomo.run_collate.m to create the combined image with the "Tomo" variable.
      save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
        'Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_array','param_records','param_sar', ...
        'Roll', 'Pitch', 'Heading','file_version');
    end
  end
end

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
