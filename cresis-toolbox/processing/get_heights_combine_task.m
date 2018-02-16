function success = get_heights_combine_task(param)
% success = get_heights_combine_task(param)
%
% Task for combining get heights temporary outputs into images and running
% surface tracking.
%
% param = struct with processing parameters
%
% Example:
%  This function is usually called via the cluster batch setup in
%  get_heights.
%
% Authors: John Paden
%
% See also: run_master.m, master.m, run_get_heights.m, get_heights.m,
%   get_heights_task.m

%% Input Checks
% =====================================================================

% Check img_comb
if numel(param.get_heights.imgs) == 1 || isempty(param.get_heights.qlook.img_comb)
  num_imgs = 1;
else
  num_imgs = length(param.get_heights.imgs);
  if length(param.get_heights.qlook.img_comb) ~= 3*(num_imgs-1)
    error('param.get_heights.qlook.img_comb not the right length. Since it is not empty, there should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
  end
end

%% Setup Processing
% =====================================================================

load(ct_filename_support(param,'','frames')); % Load "frames" variable

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);

% Quick look radar echogram output directory
qlook_out_dir = ct_filename_out(param, param.get_heights.qlook.out_path, 'CSARP_qlook');

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
  
  % recs: Determine the records for this frame
  if frm < length(frames.frame_idxs)
    recs = frames.frame_idxs(frm):frames.frame_idxs(frm+1)-1;
  else
    recs = frames.frame_idxs(frm):length(records.gps_time);
  end
  
  %% Output directory
  in_fn_dir = fullfile(qlook_out_dir, sprintf('ql_data_%03d_01_01',frm));
  
  %% Concatenate blocks for each of the images
  for img = 1:length(param.get_heights.imgs)
    
    if length(param.get_heights.imgs) == 1
      out_fn = fullfile(qlook_out_dir, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(qlook_out_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
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
    
    % Determine where breaks in processing blocks are going to occur
    %   Rename variables for readability
    block_size = param.get_heights.block_size(1);
    breaks = 1:block_size:length(recs)-0.5*block_size;
    
    % Load each block
    for break_idx = 1:length(breaks)
      
      % Determine the records for this block
      if break_idx < length(breaks)
        cur_recs_keep = [recs(breaks(break_idx)) recs(breaks(break_idx+1)-1)];
      else
        cur_recs_keep = [recs(breaks(break_idx)) recs(end)];
      end
      
      in_fn_name = sprintf('qlook_img_%02d_%d_%d.mat',img,cur_recs_keep(1),cur_recs_keep(end));
      in_fn = fullfile(in_fn_dir,in_fn_name);
      
      tmp = load(in_fn);
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
        if break_idx == 1
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
      save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
        'Elevation','Roll','Pitch','Heading','GPS_time','Data', ...
        'param_get_heights','param_records');
    else
      save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
        'Elevation','Roll','Pitch','Heading','GPS_time','Data', ...
        'param_get_heights','param_records','custom');
    end
    
    %% Create temporary output for surface tracker
    if img == 1
      Time_Surface = Time;
      Data_Surface = Data;
    end
  end
  
  if ~isempty(param.get_heights.qlook.img_comb)
    %% Combine image with previous
    if Time(end) > Time_Surface(end)
      clear combine; 
      combine.imgs     = param.get_heights.imgs;
      combine.out_path = qlook_out_path;
      combine.img_comb = param.get_heights.qlook.img_comb;
      combine.frm  = frm;
      combine.Data = Data_Surface;
      combine.Time = Time_Surface;
      [Data, Time] = img_combine(param, combine);
      Data_Surface = Data;
      Time_Surface = Time;
    end
  end
  
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
    if strcmpi(surf.method,'threshold')
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
  Surface = interp_finite(Surface,0);
  
  %% Combine images into a single image (also trim time<0 values)
  clear combine;
  combine.frm           = frm;
  combine.imgs          = param.get_heights.imgs;
  combine.out_path      = qlook_out_dir;
  combine.img_comb      = param.get_heights.qlook.img_comb;
  combine.img_comb_surf = Surface;
  combine.trim_time     = true;
  [Data, Time]          = img_combine(param, combine);  
  
  %% Save combined image output
  out_fn = fullfile(qlook_out_dir, sprintf('Data_%s_%03d.mat', ...
    param.day_seg, frm));
  fprintf('  Writing output to %s\n', out_fn);
  Data = single(Data);
  if isempty(custom)
    save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
      'Elevation','Roll','Pitch','Heading','GPS_time','Data','Surface', ...
      'param_get_heights','param_records');
  else
    save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
      'Elevation','Roll','Pitch','Heading','GPS_time','Data','Surface', ...
      'param_get_heights','param_records','custom');
  end
  
end

if param.get_heights.surf.en
  % Read the "Surface" variable from all the frames that were created
  % by this particular run of get_heights
  
  if isfield(records,'gps_time')
    num_recs = length(records.gps_time);
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
    out_fn = fullfile(qlook_out_dir, sprintf('Data_%s_%03d.mat', ...
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

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;

return;

