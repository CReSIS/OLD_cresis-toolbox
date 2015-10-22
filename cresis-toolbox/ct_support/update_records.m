function update_records(param)
% update_records(param)
%
% This function updates the records according to the param struct.
% The commands are:
% 1. Update GPS offset time
% 2. Update layers from ATM
% 3. Update layers from layer files
%
% param
%  .records_fns = for default leave empty or undefined, otherwise it is
%    a cell array of records to be considered for processing
%  .segments = cell array of day_seg strings which will be used to
%    filter the records_fns, only segments in this list will be
%    updated. If the field is not declared then all files in records_fns
%    are processed
%  .gps_time
%   .update_en = boolean enabling this update function
%   .time_offset = this number of seconds is added to the GPS time
%     recorded by the radar to correct the radar's GPS time
%     2009_Antarctica_DC8 was 1 second
%     2010_Greenland_DC8 was -14 seconds
%     2010_Greenland_P3 was 1 second
%     2010_Antarctica_DC8 was 1 second
%     2011_Greenland_P3 varies between 0 and 1 seconds
%     If left empty, the current value stored in the records file is used
%   .path = path to GPS files
%  .layers
%   .update_en = boolean enabling this update function
%   .source_type = 'ops'; % 'ops' or 'layerData'
%   .path = ''; % Required for 'layerData' (path to layer files)
%      Empty for default (equivalent to 'layerData')
%      'CSARP_post/layerData' to use posted layer data
%   .layer_names = {'surface' 'bottom'}; % Required for 'ops'
%   .gaps_dist = [300 60]; % Required for 'ops'
%  .save_changes = boolean, must be true for any changes to be saved
%
% Examples: See run_update_records.m
%
% Author: John Paden
%
% See also: run_update_records.m

update_data_files_tstart = tic;


%% Get Parameter Sheet information (uses generic flag on command worksheet)
params = read_param_xls(param.param_fn,[],'post');

param.radar_name = params(1).radar_name;
param.season_name = params(1).season_name;

if any(strcmpi(param.radar_name,{'kuband','snow','kuband2','snow2','kaband3','kuband3','snow3'}))
  
  for param_idx = 1:length(params)
    param.day_seg = params(param_idx).day_seg;
    if ~isempty(param.skip_phrase) ...
        && ~isempty(strfind(lower(params(param_idx).cmd.notes),param.skip_phrase)) ...
        || ~params(param_idx).cmd.generic
      continue;
    end
    
    fprintf('Updating record file for %s (%.1f sec)\n',...
      param.day_seg, toc(update_data_files_tstart));
    
    records_fn = ct_filename_support(param,'','records');
    records = load(records_fn);
    
    if param.gps_source.update_en
      % Update GPS file (corrects time offset too)
      
      gps_fn = ct_filename_support(param,param.gps_source.path,'gps',1);
      
      fprintf('  Using gps file %s\n', gps_fn);
      gps = load(gps_fn);
      
      tmp_records = records;
      delta_offset = params(param_idx).vectors.gps.time_offset - records.param_records.vectors.gps.time_offset
      records.param_records.vectors.gps.time_offset = params(param_idx).vectors.gps.time_offset;
      records.gps_time = records.gps_time + delta_offset;
      records.param_records.vectors.gps.time_offset = params(param_idx).vectors.gps.time_offset;
      records.lat = interp1(gps.gps_time, gps.lat, records.gps_time);
      records.lon = mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),records.gps_time)*180/pi+180, 360)-180;
      records.elev = interp1(gps.gps_time, gps.elev, records.gps_time);
      records.roll = interp1(gps.gps_time, gps.roll, records.gps_time);
      records.pitch = interp1(gps.gps_time, gps.pitch, records.gps_time);
      records.heading = mod(interp1(gps.gps_time,unwrap(gps.heading),records.gps_time)+pi,2*pi)-pi;
      records.gps_source = gps.gps_source;
    end
    
    if param.save_changes
      % Save outputs
      fprintf('  Saving records %s\n', records_fn);
      save(records_fn,'-APPEND','-struct','records','gps_time','lat','lon','elev','roll','pitch','heading','gps_source','param_records');
      create_records_aux_files(records_fn);
    else
      fprintf('  Not saving information (TEST MODE)\n');
    end
  end
  return;
end

if param.layers.update_en
  error('Please use opsCopyLayers.m to update surface variable in records file.');

  if ~isfield(param.layers, 'source_type')
    param.layers.source_type = 'layerData';
  end
  
  physical_constants;
  
  if strcmpi(param.layers.source_type, 'layerData')
    %% Load all the layer data from layerData files
    layer_path = ct_filename_out(param,param.layers.path,'CSARP_layerData',1);
    layer_fns = get_filenames(layer_path,'Data','','.mat',struct('recursive',1));
    
    fprintf('Loading layer files (%.1f sec)\n', toc(update_data_files_tstart));
    surface = [];
    bottom = [];
    elev = [];
    GPS_time = [];
    for file_idx = 1:length(layer_fns)
      tmp = load(layer_fns{file_idx});
      surface = [surface tmp.layerData{1}.value{2}.data];
      bottom = [bottom tmp.layerData{2}.value{2}.data];
      elev = [elev tmp.Elevation];
      GPS_time = [GPS_time tmp.GPS_time];
      if length(tmp.GPS_time) ~= length(tmp.layerData{1}.value{2}.data)
        fprintf('GPS_time dimensions do not match surface dimensions (ERROR!)\n')
        fprintf('  %s\n', layer_fns{file_idx});
        return;
      end
    end
    
    % Load of layer files may not be in order of the time that the data
    % was collected, so we fix that here since we will be interpolating
    % using the time axis later.
    [GPS_time,sorting_idxs] = unique(GPS_time);
    surface = surface(sorting_idxs);
    bottom = bottom(sorting_idxs);
    elev = elev(sorting_idxs);
    
  end
end

for param_idx = 1:length(params)
  param.day_seg = params(param_idx).day_seg;
  if ~isempty(param.skip_phrase) ...
      && ~isempty(strfind(lower(params(param_idx).cmd.notes),param.skip_phrase)) ...
      || ~params(param_idx).cmd.generic
    continue;
  end
  
  fprintf('Updating record file for %s (%.1f sec)\n',...
    param.day_seg, toc(update_data_files_tstart));
  
  records_fn = ct_filename_support(param,'','records');
  records = load(records_fn);
  
  if param.layers.update_en
    
    if strcmpi(param.layers.source_type, 'ops')
      %% Load the layer data from OPS
      
      fprintf('  OPS query to get layer data (%.1f sec)\n', toc(update_data_files_tstart))
      ops_sys = ct_output_dir(params(param_idx).radar_name);
      ops_param = [];
      ops_param.properties.location = params(param_idx).post.ops.location;
      ops_param.properties.season = params(param_idx).season_name;
      ops_param.properties.segment = params(param_idx).day_seg;
      ops_param.properties.return_geom = 'geog';
      ops_layer = {};
      for layer_idx = 1:length(param.layers.layer_names)
        ops_param.properties.lyr_name = param.layers.layer_names{layer_idx};
        [~,ops_layer{layer_idx}] = opsGetLayerPoints(ops_sys,ops_param);
        ops_layer{layer_idx} = ops_layer{layer_idx}.properties;
      end
      
      fprintf('  Interpolating on to records, can take ~20 minutes (%.1f sec)\n', toc(update_data_files_tstart))
      master.GPS_time = records.gps_time;
      master.Latitude = records.lat;
      master.Longitude = records.lon;
      master.Elevation = records.elev;
      lay = opsInterpLayersToMasterGPSTime(master,ops_layer,param.layers.gaps_dist);

      % Store output in variables to match layerData ingest
      GPS_time = lay.GPS_time;
      surface = lay.layerData{1}.value{2}.data;
      bottom = lay.layerData{2}.value{2}.data;
      elev = lay.Elevation;
    end

    fprintf('  Writing new data to records fields (%.1f sec)\n', toc(update_data_files_tstart))
    % Update surface and bottom variables from layer files
    % 1. Interpolate layer file variables onto record timestamps
    new_elev = interp1(GPS_time,elev,records.gps_time);
    new_surface = interp1(GPS_time,surface,records.gps_time);
    new_bottom = interp1(GPS_time,bottom,records.gps_time);
    
    tmp = records; % Store for debugging purposes
    % 2. Update surface
    if isfield(records,'surface')
      records.surface(isfinite(new_surface)) ...
        = new_surface(isfinite(new_surface));
    else
      records.surface = new_surface;
    end
    
    % 3. Update bottom
    if isfield(records,'bottom')
      records.bottom(isfinite(new_bottom)) ...
        = new_bottom(isfinite(new_bottom));
    else
      records.bottom = new_bottom;
    end
  end
  
  if param.gps_source.update_en
    % Update GPS file (corrects time offset too)
    gps_fn = ct_filename_support(param,param.gps_source.path,'gps',1);
    
    fprintf('  Using gps file %s\n', gps_fn);
    gps = load(gps_fn);
    
    tmp_records = records;
    if any(strcmpi(param.radar_name,{'accum2','mcrds'}))
      % Determine time offset delta and apply to radar time
      delta_offset = params(param_idx).vectors.gps.time_offset - records.param_records.vectors.gps.time_offset
      records.param_records.vectors.gps.time_offset = params(param_idx).vectors.gps.time_offset;
      records.raw.radar_time = records.raw.radar_time + delta_offset;

      % Find the section of the GPS file that is to be used for this record
      guard_time = 5;
      good_idxs = find(gps.comp_time >= records.raw.comp_time(1)-guard_time ...
        & gps.comp_time <= records.raw.comp_time(end)+guard_time);
      good_radar_time = gps.radar_time(good_idxs);
      good_sync_gps_time = gps.sync_gps_time(good_idxs);
      % From these good indexes, remove any repeat radar times (usually caused
      % by there being more than one NMEA string every 1 PPS
      good_idxs = 1+find(diff(good_radar_time) ~= 0);
      good_radar_time = good_radar_time(good_idxs);
      good_sync_gps_time = good_sync_gps_time(good_idxs);
      % Interpolate gps.sync_gps_time to records.gps_time using gps.radar_time
      % and records.radar_time
      records.gps_time = interp1(good_radar_time, good_sync_gps_time, ...
        records.raw.radar_time,'linear','extrap');
      
      % Interpolate to find new trajectory
      records.lat = double(interp1(gps.gps_time,gps.lat,records.gps_time));
      records.lon = double(mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),records.gps_time)*180/pi+180, 360)-180);
      records.elev = double(interp1(gps.gps_time,gps.elev,records.gps_time));
      records.roll = double(interp1(gps.gps_time,gps.roll,records.gps_time));
      records.pitch = double(interp1(gps.gps_time,gps.pitch,records.gps_time));
      records.heading = double(mod(interp1(gps.gps_time,unwrap(gps.heading),records.gps_time)+pi,2*pi)-pi);
      records.gps_source = gps.gps_source;
    else
      % Determine time offset delta and apply to radar time
      delta_offset = params(param_idx).vectors.gps.time_offset - records.param_records.vectors.gps.time_offset
      records.param_records.vectors.gps.time_offset = params(param_idx).vectors.gps.time_offset;
      records.gps_time = records.gps_time + delta_offset;
      
      % Interpolate to find new trajectory
      records.lat = interp1(gps.gps_time, gps.lat, records.gps_time);
      records.lon = mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),records.gps_time)*180/pi+180, 360)-180;
      records.elev = interp1(gps.gps_time, gps.elev, records.gps_time);
      records.roll = interp1(gps.gps_time, gps.roll, records.gps_time);
      records.pitch = interp1(gps.gps_time, gps.pitch, records.gps_time);
      records.heading = mod(interp1(gps.gps_time,unwrap(gps.heading),records.gps_time)+pi,2*pi)-pi;
      records.gps_source = gps.gps_source;
    end

  end
  
  if param.gps_time.update_en
    error('Do not use this');
    % Update GPS file (corrects time offset too)
    gps_fn = ct_filename_support(param,param.gps_time.path,'gps',1);
    
    fprintf('  Using gps file %s\n', gps_fn);
    gps = load(gps_fn);
    if isempty(param.gps_time.time_offset)
      fprintf('  Using old time offset: %.2f\n', param_records.gps.time_offset);
      time_offset = param_records.gps.time_offset;
    else
      fprintf('  Old time %.2f sec to new time %.2f\n', ...
        param_records.gps.time_offset, param.gps_time.time_offset);
      time_offset = param.gps_time.time_offset;
    end
    
    tmp_records = records;
    delta_offset = time_offset - param_records.gps.time_offset;
    param_records.gps.time_offset = time_offset;
    records.time = records.time + delta_offset;
    gps.utc_time_sod = epoch_to_sod(gps.gps_time - utc_leap_seconds(gps.gps_time(1)), param.day_seg(1:8));
    records.gps_time = interp1(gps.utc_time_sod,gps.gps_time,records.time);
    records.lat = interp1(gps.gps_time, gps.lat, records.gps_time);
    records.lon = mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),records.gps_time)*180/pi+180, 360)-180;
    records.elev = interp1(gps.gps_time, gps.elev, records.gps_time);
    records.roll = interp1(gps.gps_time, gps.roll, records.gps_time);
    records.pitch = interp1(gps.gps_time, gps.pitch, records.gps_time);
    records.heading = mod(interp1(gps.gps_time,unwrap(gps.heading),records.gps_time)+pi,2*pi)-pi;
    records.gps_source = gps.gps_source;
  end
  
  if param.gps_time.special_en
    % Update GPS file (corrects time offset too)
    gps_fn = ct_filename_support(param,param.gps_time.path,'gps',1)
    
    fprintf('  Using gps file %s\n', gps_fn);
    fprintf('  Old time %.2f sec to new time %.2f\n', ...
      param_records.gps.time_offset, param.gps_time.time_offset);
    gps = load(gps_fn);
    
    % Get the GPS seconds of day to sync to radar
    utc_time_datenum = epoch_to_datenum(gps.gps_time - utc_leap_seconds(gps.gps_time(1)));
    [year month day hour minute sec] = datevec(utc_time_datenum);
    
    UTC_sod = (day-day(1))*86400+hour*3600+minute*60+sec;  % GPS seconds of day
    
    tmp_records = records;
    records.lat = double(interp1(UTC_sod,gps.lat,records.time));
    records.lon = double(mod(interp1(UTC_sod,unwrap(gps.lon/180*pi),records.time)*180/pi+180, 360)-180);
    records.elev = double(interp1(UTC_sod,gps.elev,records.time));
    records.roll = double(interp1(UTC_sod,gps.roll,records.time));
    records.pitch = double(interp1(UTC_sod,gps.pitch,records.time));
    records.heading = double(mod(interp1(UTC_sod,unwrap(gps.heading),records.time)+pi,2*pi)-pi);
    records.gps_time = interp1(UTC_sod,gps.gps_time,records.time);
    records.gps_source = gps.gps_source;
  end
  
  if param.save_changes
    % Save outputs
    fprintf('  Saving records %s\n', records_fn);
    save(records_fn,'-struct','records');
    create_records_aux_files(records_fn);
  else
    fprintf('  Not saving information (TEST MODE)\n');
  end
end

return;



