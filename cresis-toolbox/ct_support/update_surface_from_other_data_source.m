% script update_surface_from_other_data_source
%
% Converts layers from one radar or ATM lidar and applies them to another
% radar. Uses param spreadsheet to determine which frames and segments
% to update.
%
% Author: John Paden

% =====================================================================
%% User Settings
% =====================================================================

params = read_param_xls(ct_filename_param('snow_param_2014_Greenland_P3.xls'));

sync_param.surface_type = 'radar'; % 'radar' or 'lidar'
sync_param.location = 'arctic'; % 'arctic' or 'antarctic'

sync_param.lidar.maxgap = 500;

sync_param.radar.radar_name = 'kuband';
sync_param.radar.data_type = 'qlook';
sync_param.manual = false;
sync_param.output_type = 'qlook';
sync_param.echo_type = 'qlook';
sync_param.range_offset = 0;

% =====================================================================
%% Automated Section
% =====================================================================

fprintf('=======================================================\n');
physical_constants;

for param_idx = 1:length(params)
  param = params(param_idx);
  if ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  fprintf('Updating surface %s\n', param.day_seg);
  
  %% Load the surface data source specified by sync_param.surface_type
  % surf.gps_time: seconds since Jan 1, 1970
  % surf.surface: surface height in meters referenced to WGS-84
  
  if strcmp(sync_param.surface_type,'lidar')
    %% Load LIDAR surface to replace with
    atm_dir = fullfile(sync_param.lidar.base_dir,param.season_name);
    
    global gRadar;
    atm_fns = get_filenames_atm(sync_param.location,param.day_seg(1:8),gRadar.data_support_path);
    %atm_fns = get_filenames(atm_dir,param.day_seg(3:8),'smooth_nadir','seg',struct('recursive',1));
    
    lidar = read_lidar_atm(atm_fns);
    
    % Remove NAN's from LIDAR Data
    good_lidar_idxs = ~isnan(lidar.gps_time);
    lidar.gps_time = lidar.gps_time(good_lidar_idxs);
    lidar.surface = lidar.surface(good_lidar_idxs);
    lidar.lat = lidar.lat(good_lidar_idxs);
    lidar.lon = lidar.lon(good_lidar_idxs);
    if ~isempty(isnan(lidar.surface))
      good_lidar_idxs = ~isnan(lidar.surface);
      lidar.gps_time = lidar.gps_time(good_lidar_idxs);
      lidar.surface = lidar.surface(good_lidar_idxs);
      lidar.lat = lidar.lat(good_lidar_idxs);
      lidar.lon = lidar.lon(good_lidar_idxs);
    end
    
    % ============= FIND LARGE GAPS IN LIDAR =============
    
    
    % Find Gaps Based on Along_Track Distance
    lidar.along_track = geodetic_to_along_track(lidar.lat,lidar.lon,lidar.surface);
    lidar.lidar_gaps = (diff(lidar.along_track) > sync_param.lidar.maxgap);
    
    % Prepare Data for Interp1 if gaps exist
    if ~isempty(lidar.lidar_gaps)
      % Find Start Indexes of Gaps
      lidar.lidar_gap_start_idxs = strfind([0 lidar.lidar_gaps 0],[0 1]);
      lidar.lidar_gap_start_idxs = lidar.lidar_gap_start_idxs-1;
      if lidar.lidar_gap_start_idxs(1) == 0
        lidar.lidar_gap_start_idxs(1) = 1;
      end
      
      % Find End Indexes of Gaps
      lidar.lidar_gap_end_idxs = strfind([0 lidar.lidar_gaps 0],[1 0]);
      if lidar.lidar_gap_end_idxs(end) > length(lidar.lidar_gaps)
        lidar.lidar_gap_end_idxs(end) = length(lidar.lidar_gaps);
      end
      
      % Add NaN's to lidar at Gap_IDXS (Start & End)
      lidar.surface([lidar.lidar_gap_start_idxs lidar.lidar_gap_end_idxs]) = NaN;
    end
    
    % ============= SET UP RADAR DATA FOR SYNC =============
    
    %Test for and remove duplicates in LIDAR
    if length(unique(lidar.gps_time)) < length(lidar.gps_time)
      [lidar.gps_time lidar.unique_idxs] = unique(lidar.gps_time);
      lidar.surface = lidar.surface(lidar.unique_idxs);
    end
    
    surf = lidar;
    
  elseif strcmp(sync_param.surface_type,'radar')
    %% Load RADAR surface to replace with
    
    % Get all the layer files from the corresponding day
    sync_param.radar.season_name = param.season_name;
    seg_dirs = get_filenames(ct_filename_out(sync_param.radar,sync_param.radar.data_type,'',1),param.day_seg(1:8),'','',struct('type','d'));
    
    %% Load each of the altimeter data files, store Elevation, Surface, GPS_time
    surf.surface = [];
    surf.gps_time = [];
    for dir_idx = 1:length(seg_dirs)
      fprintf('Loading segment %s\n', seg_dirs{dir_idx});
      % Get layer files in this segment directory
      layer_files = get_filenames(seg_dirs{dir_idx},'Data_','','.mat');
      
      %% Loop through each of the layer files and load it
      for layer_idx=1:length(layer_files)
        % fprintf('  File %i of %i\n', layer_idx, length(layer_files));
        warning off;
        layer = load(layer_files{layer_idx},'layerData','Surface','GPS_time','Elevation');
        warning on;
        if isfield(layer,'layerData')
          surf.gps_time = cat(2,surf.gps_time,layer.GPS_time);
          surf.surface = cat(2,surf.surface,layer.Elevation - layer.layerData{1}.value{2}.data*c/2);
        else
          surf.gps_time = cat(2,surf.gps_time,layer.GPS_time);
          surf.surface = cat(2,surf.surface,layer.Elevation - reshape(layer.Surface,[1 length(layer.Surface)])*c/2);
        end
      end
      
    end
    
    %% REALLY need to check for gaps in the data here like we do for lidar data source
    
  end

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
  
  %% Loop through each frame and correct its corresponding layer file
  for frm = param.cmd.frms
    fprintf('Processing frame %s%03d %s\n', param.day_seg, frm, datestr(now));
    
    layer_fn = fullfile(ct_filename_out(param,sync_param.output_type),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    
    warning off;
    layer = load(layer_fn,'layerData','Surface','GPS_time','Elevation');
    warning on;
    
    if isfield(layer,'layerData')
      old_surface = layer.layerData{1}.value{2}.data;
    else
      old_surface = layer.Surface;
    end
    
    %% Interpolate new surface onto old surface GPS time stamps
    warning off;
    new_Surface = interp1(surf.gps_time,surf.surface,layer.GPS_time,'linear');
    warning on;
    new_Surface = (layer.Elevation - new_Surface + sync_param.range_offset) / (c/2);
    
    % Get the indexes of newsurface to replace.
    newsurface_replace_idxs = isnan(new_Surface);
    
    % Recommended param.range_offset value
    recommended_offset = sync_param.range_offset + median(new_Surface(~isnan(new_Surface)) - old_surface(~isnan(new_Surface)))*3e8/2
    
    % Fill in the newsurface at the replacement points with radar surface.
    new_Surface(newsurface_replace_idxs)  = interp1(layer.GPS_time,old_surface,...
      layer.GPS_time(newsurface_replace_idxs));
    
    %% Save the output if the user wants to
    if sync_param.manual
      
      %% Load the echogram and plotting
      echo_fn = fullfile(ct_filename_out(param,sync_param.echo_type),sprintf('Data_%s_%03d.mat',param.day_seg,frm))

      echo = load(echo_fn);
      figure(1); clf;
      imagesc(echo.GPS_time,echo.Time*1e6,lp(echo.Data));
      colormap(1-gray(256));
      hold on;
      plot(layer.GPS_time,old_surface*1e6,'b','LineWidth',2);
      
      % ============= INSERT RADAR INTO LIDAR GAPS =============
      
      plot(layer.GPS_time,new_Surface*1e6,'r--','LineWidth',2);
      hold off;
      ylabel('propagation time (us)')
      xlabel('GPS time (sec)');
      
      cmd = input('    (1) Skip, (2) Correct layer data: ');
    else
      % Automated, just update all files
      cmd = 2;
    end
    if cmd == 2
      fprintf('    Correcting layer data\n');
      if isfield(layer,'layerData')
        layer.layerData{1}.value{2}.data = new_Surface;
        save(layer_fn,'-append','-struct','layer','layerData');
      else
        layer.Surface = new_Surface;
        save(layer_fn,'-append','-struct','layer','Surface');
      end
    end

  end
  
end

return
