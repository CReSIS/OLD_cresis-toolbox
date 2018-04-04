function success = combine_wf_chan_combine_task(param)
% success = combine_wf_chan_combine_task(param)
%
% This script combines the receive channels and outputs the result
% for each waveform. It also combines the waveforms. It takes in
% f-k files one directory at a time and:
%  1. combines the receive channels
%  2. concatenates the results
%  3. square-law detects the data, abs()^2
%  4. takes incoherent averages (multilooks data)
%  5. saves the result in a new directory
%  6. Combines the waveforms
%
% The assumption is that the directories in the input_path are named
% using the following convention:
%   PROC-TYPE-STRING_data_#{_SUBAPERTURE-STRING}
% where
%   PROC-TYPE-STRING can be 'fk','tdbp', or 'pc' for f-k migrated,time domain
%   back projected,and pulse compressed respectively ('fk' and tdbp supported)
%   _data_ is always present
%   #, \d+: one or more numbers
%   _SUBAPERTURE-STRING, {_[mp]\d\.\d}: optional subaperture string
% Examples:
%   fk_data_01_01: f-k migrated, frame 1, subaperture 1
%   fk_data_04_02: f-k migrated, frame 4, subaperture 2
%   fk_data_01_03: f-k migrated, frame 1, subaperture 3
%   pc_data_01: pulse compressed only, frame 1
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
% See also: run_master.m, master.m, run_combine_wf_chan.m, combine_wf_chan.m,
%   combine_wf_chan_task.m

%% Setup processing
% =====================================================================

% Create output directory path
combine_out_dir = ct_filename_out(param, param.combine.out_path);
if ~exist(combine_out_dir,'dir')
  mkdir(combine_out_dir);
end

load(ct_filename_support(param,'','frames')); % Load "frames" variable

% Load records file
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);
% Apply presumming
if param.csarp.presums > 1
  records.lat = fir_dec(records.lat,param.csarp.presums);
  records.lon = fir_dec(records.lon,param.csarp.presums);
  records.elev = fir_dec(records.elev,param.csarp.presums);
  records.roll = fir_dec(records.roll,param.csarp.presums);
  records.pitch = fir_dec(records.pitch,param.csarp.presums);
  records.heading = fir_dec(records.heading,param.csarp.presums);
  records.gps_time = fir_dec(records.gps_time,param.csarp.presums);
  records.surface = fir_dec(records.surface,param.csarp.presums);
end
along_track_approx = geodetic_to_along_track(records.lat,records.lon,records.elev);

%% Combine chunks into each frame
% =====================================================================
for frm_idx = 1:length(param.cmd.frms);
  frm = param.cmd.frms(frm_idx);
  
  if ct_proc_frame(frames.proc_mode(frm),param.combine.frm_types)
    fprintf('combine %s_%03i (%i of %i) %s\n', param.day_seg, frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
  else
    fprintf('Skipping frame %s_%03i (no process frame)\n', param.day_seg, frm);
    continue;
  end
  
  % Temporary output directory for uncombined array processed images
  array_fn_dir = fullfile(ct_filename_out(param, param.combine.array_path), ...
    sprintf('%s_%03d', param.combine.method, frm));

  % Current frame goes from the start record specified in the frames file
  % to the record just before the start record of the next frame.  For
  % the last frame, the stop record is just the last record in the segment.
  start_rec = ceil(frames.frame_idxs(frm)/param.csarp.presums);
  if frm < length(frames.frame_idxs)
    stop_rec = ceil((frames.frame_idxs(frm+1)-1)/param.csarp.presums);
  else
    stop_rec = length(records.gps_time);
  end
  
  % Determine length of the frame
  frm_dist = along_track_approx(stop_rec) - along_track_approx(start_rec);
  
  % Determine number of chunks for this frame
  num_chunks = round(frm_dist / param.combine.chunk_len);

  %% Loop through all the images
  for img = 1:length(param.combine.imgs)
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
    Topography = [];
    for chunk_idx = 1:num_chunks
      array_fn = fullfile(array_fn_dir, sprintf('img_%02d_chk_%03d.mat', img, chunk_idx));
      tmp = load(array_fn);
      Time = tmp.Time;
      Latitude = [Latitude double(tmp.Latitude)];
      Longitude = [Longitude double(tmp.Longitude)];
      Elevation = [Elevation double(tmp.Elevation)];
      Roll = [Roll double(tmp.Roll)];
      Pitch = [Pitch double(tmp.Pitch)];
      Heading = [Heading double(tmp.Heading)];
      GPS_time = [GPS_time tmp.GPS_time];
      Surface = [Surface double(tmp.Surface)];
      Bottom = [Bottom double(tmp.Bottom)];
      Data = [Data tmp.Data];
      param_records = tmp.param_records;
      param_csarp = tmp.param_csarp;
      if chunk_idx == 1
        param_combine = tmp.param_combine;
        param_combine.array_param.fcs{1}{1}.x = tmp.param_combine.array_param.fcs{1}{1}.x(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.y = tmp.param_combine.array_param.fcs{1}{1}.y(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.z = tmp.param_combine.array_param.fcs{1}{1}.z(:,tmp.param_combine.array_param.lines);
        param_combine.array_param.fcs{1}{1}.origin = tmp.param_combine.array_param.fcs{1}{1}.origin(:,tmp.param_combine.array_param.lines);
      else
        % Concatenate the fcs field
        param_combine.array_param.fcs{1}{1}.x = [param_combine.array_param.fcs{1}{1}.x tmp.param_combine.array_param.fcs{1}{1}.x(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.y = [param_combine.array_param.fcs{1}{1}.y tmp.param_combine.array_param.fcs{1}{1}.y(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.z = [param_combine.array_param.fcs{1}{1}.z tmp.param_combine.array_param.fcs{1}{1}.z(:,tmp.param_combine.array_param.lines)];
        param_combine.array_param.fcs{1}{1}.origin = [param_combine.array_param.fcs{1}{1}.origin tmp.param_combine.array_param.fcs{1}{1}.origin(:,tmp.param_combine.array_param.lines)];
      end
      if isfield(tmp,'Topography')
        %         3D-surface is present so concatenate it too
        %         Topography = cat(3,Topography,tmp.Topography);
        %         Concatenate all the fields under struct Topography: valR, bins, val, freq
        %         and img.
        fields = fieldnames(tmp.Topography);
        if chunk_idx == 1
          for field_idx = 1:length(fields)
            Topography.(fields{field_idx}) = tmp.Topography.(fields{field_idx});
          end
        else
          for field_idx = 1:length(fields)
            max_dim = length(size(tmp.Topography.(fields{field_idx})));
            Topography.(fields{field_idx}) = cat(max_dim,Topography.(fields{field_idx}),tmp.Topography.(fields{field_idx}));
          end
        end
        
      end
    end
    
    % =====================================================================
    % Save output
    if length(param.combine.imgs) == 1
      out_fn = fullfile(combine_out_dir, sprintf('Data_%s_%03d.mat', ...
        param.day_seg, frm));
    else
      out_fn = fullfile(combine_out_dir, sprintf('Data_img_%02d_%s_%03d.mat', ...
        img, param.day_seg, frm));
    end
    fprintf('  Writing output to %s\n', out_fn);
    if isempty(Topography)
      % Do not save 3D surface
      save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
        'Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_combine','param_records','param_csarp', ...
        'Roll', 'Pitch', 'Heading');
    else
      % Save 3D surface
      save('-v7.3',out_fn,'Topography','Time','Latitude', ...
        'Longitude','Elevation','GPS_time','Data','Surface','Bottom', ...
        'param_combine','param_records','param_csarp', ...
        'Roll', 'Pitch', 'Heading');
    end
  end
  
  if isempty(param.combine.img_comb)
    % No image combining is required
    continue;
  end
  
  if length(param.combine.img_comb) ~= 3*(length(param.combine.imgs)-1)
    warning('param.combine.img_comb not the right length. There should be 3 entries for each image combination interface ([Tpd second image for surface saturation, -inf for second image blank, Tpd first image to avoid roll off] is typical). Set correctly here and update param spreadsheet before dbcont.');
    keyboard
  end
  
  %% Combine images
  if isempty(param.combine.img_comb_layer_params)
    layers.gps_time = GPS_time;
    layers.twtt = interp_finite(Surface,0);
  end
  param.load.frm        = frm;
  [Data, Time]          = img_combine(param, 'combine', layers);
  
  %% Save output
  out_fn = fullfile(combine_out_dir, sprintf('Data_%s_%03d.mat', ...
    param.day_seg, frm));
  fprintf('  Writing output to %s\n', out_fn);
  save('-v7.3',out_fn,'Time','Latitude','Longitude', ...
    'Elevation','GPS_time','Data','Surface','Bottom', ...
    'param_combine','param_records','param_csarp', ...
    'Roll', 'Pitch', 'Heading');
end

fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
