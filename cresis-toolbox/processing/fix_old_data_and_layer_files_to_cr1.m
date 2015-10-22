% Fix data and layers files to the new CR1 format
%
% Converts old files to the new CR1 format that DO NOT HAVE records,
% frames, etc. (i.e. radars such as ICARDS, ACORDS)
% 1. Run this script to create dataset to load into param spreadsheet.
%    Creates .csv file for this purpose.
% 2. This script creates records and frames files from posted products
% 3. This script then copies frames, echograms, and layer data to cr1

old_data_dir = '/cresis/scratch2/mdce/'
radar_name = 'acords'
season_name = '2003_Greenland_P3'

% old_data_dir = '/cresis/scratch2/mdce/'
% radar_name = 'acords'
% season_name = '2004_Antarctica_P3chile'

param_fn_dir = '/users/paden/scripts/branch/params-cr1/'

% =========================================================================
%% Automated Section

base_dir = fullfile(old_data_dir,radar_name,season_name)

layer_dir = fullfile(base_dir,'CSARP_layerData')

data_dir = fullfile(base_dir,'CSARP_standard')

day_segs_dir = get_filenames(data_dir,'[0-9]','','',struct('type','d'));
day_segs_dir = day_segs_dir(1:end);

param = [];
param_idx = 0;
for old_param_idx = 1:length(day_segs_dir)
  [~,day_seg] = fileparts(day_segs_dir{old_param_idx});
  fprintf('Convert to CR1 %s\n', day_segs_dir{old_param_idx});
  data_fns = get_filenames(day_segs_dir{old_param_idx},'Data_','','.mat');
  GPS_time = [];
  Latitude = [];
  Longitude = [];
  Elevation = [];
  Surface = [];
  Frame = [];
  Data_Size = [];
  for data_idx=1:length(data_fns)
    fprintf('  %s\n', data_fns{data_idx});
    data = load(data_fns{data_idx},'GPS_time','Latitude','Longitude','Elevation','Surface','Data');
    frm = str2double(data_fns{data_idx}(end-6:end-4));
    if frm ~= data_idx
      keyboard;
    end
    GPS_time = cat(2,GPS_time,data.GPS_time);
    Latitude = cat(2,Latitude,data.Latitude);
    Longitude = cat(2,Longitude,data.Longitude);
    Elevation = cat(2,Elevation,data.Elevation);
    Surface = cat(2,Surface,data.Surface);
    Frame = cat(2,Frame,frm*ones(size(data.GPS_time)));
    Data_Size = cat(2,Data_Size,size(data.Data,1)*ones(size(data.GPS_time)));
  end
  
  outlier_expected = medfilt1(Latitude,3);
  outlier_check = abs(outlier_expected - Latitude) > 0.1;
  if any(outlier_check)
    warning('Found %d outliers in Latitude', sum(outlier_check));
    find(outlier_check)
  end
  lat_outlier_check = outlier_check;
  Latitude(outlier_check) = outlier_expected(outlier_check);
  
  outlier_expected = medfilt1(Longitude,3);
  outlier_check = abs(outlier_expected - Longitude) > 0.1;
  for outlier_idx = find(outlier_check)
    Longitude(outlier_idx) = interp1(GPS_time(outlier_idx+[-1 1]), ...
      Longitude(outlier_idx+[-1 1]),GPS_time(outlier_idx));
  end
  if any(outlier_check)
    warning('Found %d outliers in Longitude', sum(outlier_check));
    find(outlier_check)
  end
  if any(lat_outlier_check ~= outlier_check)
    warning('Latitude and longitude errors usually line up. Did not. Probably an error in one that did not get caught');
    keyboard
  end
  
  outlier_expected = medfilt1(Elevation,3);
  outlier_check = abs(outlier_expected - Elevation) > 100;
  Elevation(outlier_check) = outlier_expected(outlier_check);
  if any(outlier_check)
    warning('Found %d outliers in Elevation. Not expecting these.', sum(outlier_check));
    find(outlier_check)
    keyboard
  end
  
  if strcmp(day_seg,'20021210_01')
    % Temporary hack to get it to work
    figure(1); clf;
    plot(GPS_time);
    GPS_time(3856:end) = GPS_time(3856:end) - diff(GPS_time(3856+[-1 0]));
    GPS_time(4132:end) = GPS_time(4132:end) - diff(GPS_time(4132+[-1 0]));
    hold on;
    plot(GPS_time,'r');
    hold off;
    keyboard
  end
  
  % Remove all negative jumps in time
  bad_mask = zeros(size(GPS_time));
  cur_time = GPS_time(1);
  good_state = 1;
  for gps_idx = 2:length(GPS_time)
    if GPS_time(gps_idx) <= cur_time
      if good_state == 1
      end
      bad_mask(gps_idx) = 1;
      good_state = 0;
    else
      if good_state == 0
      end
      cur_time = GPS_time(gps_idx);
      good_state = 1;
    end
  end
  GPS_time = GPS_time(~bad_mask);
  Latitude = Latitude(~bad_mask);
  Longitude = Longitude(~bad_mask);
  Elevation = Elevation(~bad_mask);
  Surface = Surface(~bad_mask);
  Data_Size = Data_Size(~bad_mask);
  Frame_Good = Frame(~bad_mask);
  
  % plot(diff(GPS_time))
  along_track = geodetic_to_along_track(Latitude,Longitude,Elevation);
  
  % Check negative time jump
  diff_GPS_time = diff(GPS_time);
  neg_gps_idxs = find(diff_GPS_time < 0)
  diff_GPS_time(neg_gps_idxs)
  big_gps_idxs = find(diff_GPS_time > 8)
  diff_GPS_time(big_gps_idxs)
  round(diff_GPS_time(big_gps_idxs))
  fprintf('min: %f max: %f\n', min(diff_GPS_time), max(diff_GPS_time))
  
  frame_idxs = [1 find(diff(Frame_Good) > 0)]
  
  diff_along_track = diff(along_track);
  big_along_idxs = find(diff_along_track > 1e3)
  round(diff_along_track(big_along_idxs)/1e3)
  
  if length(GPS_time) ~= length(Latitude) ...
      || length(Longitude) ~= length(Latitude) ...
      || length(Longitude) ~= length(Elevation) ...
      || length(Surface) ~= length(Elevation)
    keyboard;
  end
  
  seg_boundaries = intersect(big_along_idxs,big_gps_idxs)
  error_idxs = setdiff(big_along_idxs,big_gps_idxs)
  if ~isempty(error_idxs)
    diff_GPS_time(error_idxs)
    diff_along_track(error_idxs)
    warning('Big GPS/along track jump(s) that do not line up');
    keyboard
    if strcmp(day_seg,'20021210_01')
      % Temporary hack to get it to work
      seg_boundaries = big_along_idxs;
    end
    if any(strcmp(day_seg,{'20030509_01','20030512_01'}))
      % Temporary hack to get it to work
      seg_boundaries = big_along_idxs;
    end
  end
  
  data_size_change_idxs = find(abs(diff(Data_Size)) > 0)
  seg_boundaries = union(data_size_change_idxs,seg_boundaries);
  
  cur_param_idx = 0;
  seg_boundaries = [0 seg_boundaries];
  for segment = 1:length(seg_boundaries)
    param_idx = param_idx + 1;
    cur_param_idx = cur_param_idx + 1;
    if strcmpi(radar_name,'wise')
      param(param_idx).day_seg = day_seg;
    else
      day_seg(end-1:end) = sprintf('%02d',cur_param_idx);
      param(param_idx).day_seg = day_seg;
    end
    param(param_idx).radar_name = radar_name;
    param(param_idx).season_name = season_name;
    fprintf('Creating segment %s\n', param(param_idx).day_seg)
    if segment == length(seg_boundaries)
      seg_idxs = seg_boundaries(segment)+1 : length(GPS_time);
    else
      seg_idxs = seg_boundaries(segment)+1 : seg_boundaries(segment+1);
    end
    %% Create records file
    records = [];
    records.gps_time = GPS_time(seg_idxs);
    records.lat = Latitude(seg_idxs);
    records.lon = Longitude(seg_idxs);
    records.elev = Elevation(seg_idxs);
    records.roll = zeros(size(records.lat));
    records.pitch = zeros(size(records.lat));
    records.heading = zeros(size(records.lat));
    records.surface = Surface(seg_idxs);
    records_fn = ct_filename_support(param(param_idx),'','records');
    records_fn_dir = fileparts(records_fn);
    if ~exist(records_fn_dir,'dir')
      mkdir(records_fn_dir);
    end
    records.notes = sprintf('Records file artificially created from data files on %s', datestr(now));
    records.ver = 3;
    records.gps_source = 'NMEA-field';
    records.param_records.vectors.gps.time_offset = 0;
    save(records_fn,'-struct','records');
    
    %% Create frames
    seg_along_track = along_track(seg_idxs);
    seg_along_track = seg_along_track - seg_along_track(1);
    frames.frame_idxs = [];
    frame_along_track_length = 50e3;
    along_track_breaks = 0:frame_along_track_length:seg_along_track(end);
    if seg_along_track(end) - along_track_breaks(end) < frame_along_track_length/2 ...
        && length(along_track_breaks) > 1
      along_track_breaks = along_track_breaks(1:end-1);
    end
    cur_frame_idx = 0;
    frames = struct('frame_idxs',[]);
    data_fn_dir = ct_filename_out(param(param_idx),'','CSARP_standard');
    if ~exist(data_fn_dir,'dir')
      mkdir(data_fn_dir);
    end
    layer_fn_dir = ct_filename_out(param(param_idx),'','CSARP_layerData');
    if ~exist(layer_fn_dir,'dir')
      mkdir(layer_fn_dir);
    end
    for frm = 1:length(along_track_breaks)
      if frm == length(along_track_breaks)
        frame_idxs = cur_frame_idx+1 : length(seg_along_track);
      else
        frame_idxs = cur_frame_idx+1 : find(seg_along_track < along_track_breaks(frm+1),1,'last');
      end
      frames.frame_idxs(frm) = frame_idxs(1);
      fprintf('  Creating frame %s_%03d from %d-%d\n', ...
        param(param_idx).day_seg, frm, frame_idxs([1 end]));
      if strcmpi('20021210_03_007',sprintf('%s_%03d',param(param_idx).day_seg, frm))
        keyboard
      end
      % Load original file data
      old_data_frms = unique(Frame_Good(seg_idxs(frame_idxs)));
      old_data_frms = old_data_frms(1):old_data_frms(end); % Just in case any frames were completely bad... we still will load them so that indexing works out
      old = [];
      old.GPS_time = [];
      old.Latitude = [];
      old.Longitude = [];
      old.Elevation = [];
      old.Surface = [];
      old.Bottom = [];
      old.Surface_Quality = [];
      old.Bottom_Quality = [];
      old.Surface_Manual = [];
      old.Bottom_Manual = [];
      old.Data = [];
      old.Frame = [];
      for data_idx = old_data_frms
        %fprintf('  %s\n', data_fns{data_idx});
        data = load(data_fns{data_idx});
        old_frm = str2double(data_fns{data_idx}(end-6:end-4));
        if old_frm ~= data_idx
          keyboard;
        end
        [~,day_seg] = fileparts(day_segs_dir{old_param_idx});
        [~,old_data_fn] = fileparts(data_fns{data_idx});
        layer_fn = fullfile(layer_dir,day_seg,[old_data_fn '.mat']);
        if exist(layer_fn,'file')
          old_layer = load(layer_fn);
          if length(old_layer.layerData{1}.value{2}.data) ...
              ~= length(data.Surface)
            keyboard
          end
          old.Surface = cat(2,old.Surface,old_layer.layerData{1}.value{2}.data);
          old.Bottom = cat(2,old.Bottom,old_layer.layerData{2}.value{2}.data);
          old.Surface_Quality = cat(2,old.Surface_Quality,old_layer.layerData{1}.quality);
          old.Bottom_Quality = cat(2,old.Bottom_Quality,old_layer.layerData{2}.quality);
          old.Surface_Manual = cat(2,old.Surface_Manual,old_layer.layerData{1}.value{1}.data);
          old.Bottom_Manual = cat(2,old.Bottom_Manual,old_layer.layerData{2}.value{1}.data);
        else
          old.Surface = cat(2,old.Surface,data.Surface);
          old.Bottom = cat(2,old.Bottom,data.Bottom);
          old.Surface_Quality = cat(2,old.Surface_Quality,ones(size(data.Surface)));
          old.Bottom_Quality = cat(2,old.Bottom_Quality,ones(size(data.Bottom)));
          old.Surface_Manual = cat(2,old.Surface_Manual,NaN*data.Surface);
          old.Bottom_Manual = cat(2,old.Bottom_Manual,NaN*data.Bottom);
        end
        
        old.GPS_time = cat(2,old.GPS_time,data.GPS_time);
        old.Latitude = cat(2,old.Latitude,data.Latitude);
        old.Longitude = cat(2,old.Longitude,data.Longitude);
        old.Elevation = cat(2,old.Elevation,data.Elevation);
        old.Data = cat(2,old.Data,data.Data);
        old.Frame = cat(2,old.Frame,frm*ones(size(data.GPS_time)));
      end
      old.Time = data.Time;
      old.Depth = data.Depth;
      old.param_records = records.param_records;
      old.param_records.gps_source = records.gps_source;
      
      % Remove bad indices (bad_mask)
      loaded_idxs = find(Frame >= old_data_frms(1) & Frame <= old_data_frms(end));
      old.GPS_time = old.GPS_time(~bad_mask(loaded_idxs));
      old.Latitude = old.Latitude(~bad_mask(loaded_idxs));
      old.Longitude = old.Longitude(~bad_mask(loaded_idxs));
      old.Elevation = old.Elevation(~bad_mask(loaded_idxs));
      old.Surface = old.Surface(~bad_mask(loaded_idxs));
      old.Bottom = old.Bottom(~bad_mask(loaded_idxs));
      old.Surface_Quality = old.Surface_Quality(~bad_mask(loaded_idxs));
      old.Bottom_Quality = old.Bottom_Quality(~bad_mask(loaded_idxs));
      old.Surface_Manual = old.Surface_Manual(~bad_mask(loaded_idxs));
      old.Bottom_Manual = old.Bottom_Manual(~bad_mask(loaded_idxs));
      old.Data = old.Data(:,~bad_mask(loaded_idxs));
      
      % Select indices for this frame: seg_idxs(frame_idxs)
      keep_idxs = old.GPS_time >= GPS_time(seg_idxs(frame_idxs(1))) ...
        & old.GPS_time <= GPS_time(seg_idxs(frame_idxs(end)));
      old.GPS_time = old.GPS_time(keep_idxs);
      old.Latitude = old.Latitude(keep_idxs);
      old.Longitude = old.Longitude(keep_idxs);
      old.Elevation = old.Elevation(keep_idxs);
      old.Surface = old.Surface(keep_idxs);
      old.Bottom = old.Bottom(keep_idxs);
      old.Surface_Quality = old.Surface_Quality(keep_idxs);
      old.Bottom_Quality = old.Bottom_Quality(keep_idxs);
      old.Surface_Manual = old.Surface_Manual(keep_idxs);
      old.Bottom_Manual = old.Bottom_Manual(keep_idxs);
      old.Data = old.Data(:,keep_idxs);
      
      old.Latitude = Latitude(seg_idxs(frame_idxs));
      old.Longitude = Longitude(seg_idxs(frame_idxs));

      % Create layer file
      layer_fn = fullfile(layer_fn_dir,sprintf('Data_%s_%03d.mat', ...
        param(param_idx).day_seg,frm));
      layer.GPS_time = old.GPS_time;
      layer.Latitude = old.Latitude;
      layer.Longitude = old.Longitude;
      layer.Elevation = old.Elevation;
      layer.layerData{1}.value{1}.data = old.Surface_Manual;
      layer.layerData{1}.value{2}.data = old.Surface;
      layer.layerData{1}.quality = old.Surface_Quality;
      layer.layerData{2}.value{1}.data = old.Bottom_Manual;
      layer.layerData{2}.value{2}.data = old.Bottom;
      layer.layerData{2}.quality = old.Bottom_Quality;
      fprintf('   %s\n', layer_fn);
      save(layer_fn, '-struct', 'layer');
      
      old = rmfield(old,'Surface_Quality');
      old = rmfield(old,'Bottom_Quality');
      old = rmfield(old,'Surface_Manual');
      old = rmfield(old,'Bottom_Manual');
      old = rmfield(old,'Frame');
      % Create data file
      data_fn = fullfile(data_fn_dir,sprintf('Data_%s_%03d.mat', ...
        param(param_idx).day_seg,frm));
      fprintf('   %s\n', data_fn);
      save(data_fn,'-struct','old');
      
      cur_frame_idx = frame_idxs(end);
    end
    % Create frames file
    frames.proc_mode = zeros(size(frames.frame_idxs));
    frames_fn = ct_filename_support(param(param_idx),'','frames');
    frames_fn_dir = fileparts(frames_fn);
    if ~exist(frames_fn_dir,'dir')
      mkdir(frames_fn_dir);
    end
    save(frames_fn,'frames');
  end
  
  %plot(GPS_time)
end

%% Save param file as csv
param_fn = fullfile(param_fn_dir,sprintf('%s_param_%s.csv', ...
  ct_output_dir(radar_name), season_name));
fprintf('\nWriting %s\n', param_fn);
fid = fopen(param_fn,'w');
fprintf(fid,'Version,3.0\n');
fprintf(fid,'Radar,%s\n', radar_name);
fprintf(fid,'Season,%s\n', season_name);
fprintf(fid,'Date,\n');
fprintf(fid,'YYYYMMDD,Segment\n');
for param_idx = 1:length(param)
  fprintf(fid,'%s,%s\n', param(param_idx).day_seg(1:8), param(param_idx).day_seg(end-1:end));
end
fclose(fid);

return;
