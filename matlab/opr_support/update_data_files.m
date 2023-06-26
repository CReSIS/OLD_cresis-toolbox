% script update_data_files
%
% Used to update GPS OR layer information in data files (e.g. qlook, csarp-combined,
% standard, layerData, etc). These are the outputs of the qlook processing for the
% FMCW and accum radars and the outputs of the combine_chan_wf for
% MCoRDS.
%
% GPS:
% Should ONLY be used with field products. All other products require
% the data to be processed with the GPS data.
% IMPORTANT: make_layer_files SHOULD BE RUN AFTER THIS FUNCTION with
% do_not_overwrite to true and update_gps to 'no-layer-correction'
%
% Layer:
% Updates SAR processed files (useful for doing 3d-imaging without having
% to reprocess all the files)
%
% Author: John Paden

% ======================================================================
% User Settings
% ======================================================================

% Parameters spreadsheet to use for updating
%   1. Segment and frame list are taken from the parameter sheet
%   2. For GPS update, GPS time offsets are pulled from the parameter sheet
% param.param_fn = ct_filename_param('kuband_param_2015_Greenland_LC130.xls');
% param.param_fn = ct_filename_param('rds_param_2015_Greenland_LC130.xls');
% param.param_fn = ct_filename_param('snow_param_2015_Greenland_LC130.xls');
param.param_fn = ct_filename_param('snow_param_2017_Greenland_P3.xls');

% param.skip_phrase: All segments will be skipped with this phrase in their
%   verification field. "do not process" is the standard. Leave this field
%   blank to do all segments.
param.skip_phrase = '';

% param.types: cell array of strings specifying output directories to
%   update, e.g. {'csarp-combined','qlook','mvdr','CSARP_post/mvdr'}
% param.types = {'CSARP_post/qlook','CSARP_post/mvdr','CSARP_post/standard'};
param.types = {'qlook'};

% param.update_mode: string 'gps' or 'layer'
param.update_mode = 'gps';

% param.save_en: Logical, For debugging purposes, you can turn the file save on/off
param.save_en = true;

% param.debug_level: default is 1, anything higher is for debugging
param.debug_level = 1;

% ======================================================================
% Automated Section
% ======================================================================

update_data_files_tstart = tic;
physical_constants;

if strcmpi(param.update_mode,'gps')
  %% =====================================================================
  % GPS Update Section
  % ======================================================================
  
  % params: structure of parameters for each segment
  params = read_param_xls(param.param_fn);
  params = ct_set_params(params,'cmd.generic',0);
  params = ct_set_params(params,'cmd.generic',1,'day_seg','20170426_02');
%   params = ct_set_params(params,'cmd.frms',[1:5]);
  
  param.radar_name = params(1).radar_name;
  param.season_name = params(1).season_name;
  
  for param_idx = 1:length(params)
    param.day_seg = params(param_idx).day_seg;
    if ~isempty(param.skip_phrase) ...
        && ~isempty(strfind(lower(params(param_idx).cmd.notes),param.skip_phrase)) ...
        || ~params(param_idx).cmd.generic
      continue;
    end
    fprintf('Updating data files for %s (%.1f sec)\n',...
      param.day_seg, toc(update_data_files_tstart));
    
    gps_fn = ct_filename_support(param,'','gps',true);
    gps = load(gps_fn);
    
    for type_idx = 1:length(param.types)
      
      data_dir = ct_filename_out(param,param.types{type_idx});
      
      if ~exist(data_dir,'dir')
        continue;
      end
      
      fprintf('  Output type %d of %d %s (%.1f sec)\n', type_idx, ...
        length(param.types), param.types{type_idx}, ...
        toc(update_data_files_tstart));
      
      data_fns = get_filenames(data_dir,'Data','','.mat');
      for data_idx = 1:length(data_fns)
        % If the frames list is not empty, then we need to check if this
        % frame is in the list. If it is not, we skip ahead.
        [data_fn_path data_fn_name] = fileparts(data_fns{data_idx});
        if data_fn_name(6) == 'i'
          frm_str = data_fn_name(13:27);
          img_idx = str2double(data_fn_name(10:11));
        else
          frm_str = data_fn_name(6:20);
          img_idx = 1;
        end
        frm = str2double(frm_str(end-2:end));
        if ~isempty(params(param_idx).cmd.frms)
          if isempty(find(frm == params(param_idx).cmd.frms))
            continue;
          end
        end
        fprintf('    File %d of %d %s (%.1f sec)\n', data_idx, length(data_fns), ...
          data_fns{data_idx}, toc(update_data_files_tstart));
        clear param_records;
        if param.debug_level > 1
          tmp = load(data_fns{data_idx});
        end
        warning off;
        clear GPS_time param_records param_qlook param_array Elevation_Correction Time;
        load(data_fns{data_idx},'GPS_time','param_records','param_qlook','param_array','param_sar','Elevation_Correction','Time');
        warning on;
        if isfield(param_records,'vectors')
          warning('Old parameter format with "vectors" field found in records file. Updating format.');
          param_records.records.file.start_idx = param_records.vectors.file.start_idx;
          param_records.records.file.stop_idx = param_records.vectors.file.stop_idx;
          param_records.records.file.base_dir = param_records.vectors.file.base_dir;
          param_records.records.file.board_folder_name = param_records.vectors.file.adc_folder_name;
          param_records.records.file.prefix = param_records.vectors.file.file_prefix;
          param_records.records.gps = param_records.vectors.gps;
          
          param_records.records.file.version = param_records.records.file_version;
          param_records.records.frames.geotiff_fn = param_records.records.geotiff_fn;
          param_records.records.frames.mode = param_records.records.frame_mode;
          param_records = rmfield(param_records,'vectors');
          try
            param_records.records.file = rmfield(param_records.records.file,'adcs');
          end
          try
            param_records.records.file = rmfield(param_records.records.file,'adc_headers');
          end
          try
            param_records.records.gps = rmfield(param_records.records.gps,'verification');
          end
          try
            param_records.records.gps = rmfield(param_records.records.gps,'fn');
          end
          if isfield(param_records.records.gps,'utc_time_halved') && isempty(param_records.records.gps.utc_time_halved)
            param_records.records.gps = rmfield(param_records.records.gps,'utc_time_halved');
          end
          try
            param_records.records = rmfield(param_records.records,'records_fn');
          end
          try
            param_records.records = rmfield(param_records.records,'force_all');
          end
          try
            param_records.records = rmfield(param_records.records,'frames_fn');
          end
          try
            param_records.records = rmfield(param_records.records,'geotiff_fn');
          end
          try
            param_records.records = rmfield(param_records.records,'frame_mode');
          end
          try
            param_records.records = rmfield(param_records.records,'debug_level');
          end
          try
            param_records.records = rmfield(param_records.records,'file_version');
          end
          try
            param_records.records = rmfield(param_records.records,'tmp_fn_uses_adc_folder_name');
          end
          try
            param_records.records = rmfield(param_records.records,'manual_time_correct');
          end
        end
        
        % Desired minus the old
        delta_offset = params(param_idx).records.gps.time_offset - param_records.records.gps.time_offset;
        if delta_offset ~= 0
          fprintf('    Delta GPS time offset is %.1f sec\n', delta_offset);
        end
        GPS_time = GPS_time + delta_offset;
        param_records.records.gps.time_offset = params(param_idx).records.gps.time_offset;
        param_records.gps_source = gps.gps_source;
        
        if exist('param_qlook','var')
          qlook_based = true;
          new_param = param_qlook;
        elseif exist('param_sar','var')
          qlook_based = false;
          new_param = param_sar;
        else
          error('Could not find param_qlook or param_sar in file.');
        end
        
        new_param.radar.lever_arm_fh = params(param_idx).radar.lever_arm_fh;
        if isempty(params(param_idx).radar.lever_arm_fh)
          lever_arm_available = false;
        else
          lever_arm_available = true;
          % Default values to use
          wf = abs(new_param.load.imgs{img_idx}(1,1));
          adc = abs(new_param.load.imgs{img_idx}(1,2));
          lever_arm_fh = new_param.radar.lever_arm_fh;
          trajectory_param = struct('rx_path', new_param.radar.wfs(wf).rx_paths(adc), ...
            'tx_weights', new_param.radar.wfs(wf).tx_weights, 'lever_arm_fh', lever_arm_fh);
          for tmp_wf_adc_idx = 2:size(new_param.load.imgs{img_idx},1)
            tmp_wf = abs(new_param.load.imgs{img_idx}(tmp_wf_adc_idx,1));
            tmp_adc = abs(new_param.load.imgs{img_idx}(tmp_wf_adc_idx,2));
            trajectory_param.rx_path(tmp_wf_adc_idx) = new_param.radar.wfs(tmp_wf).rx_paths(tmp_adc);
          end
        end
        
        new_gps.lat = interp1(gps.gps_time, gps.lat, GPS_time);
        new_gps.lon = interp1(gps.gps_time, gps.lon, GPS_time);
        new_gps.elev = interp1(gps.gps_time, gps.elev, GPS_time);
        new_gps.roll = interp1(gps.gps_time, gps.roll, GPS_time);
        new_gps.pitch = interp1(gps.gps_time, gps.pitch, GPS_time);
        new_gps.heading = interp1(gps.gps_time, gps.heading, GPS_time);
        new_gps.gps_source = gps.gps_source;
        if lever_arm_available
          trajectory_param.gps_source = new_gps.gps_source;
          trajectory_param.radar_name = params(param_idx).radar_name;
          trajectory_param.season_name = params(param_idx).season_name;
          new_gps = trajectory_with_leverarm(new_gps,trajectory_param);
        end
        Latitude = new_gps.lat;
        Longitude = new_gps.lon;
        if exist('Elevation_Correction','var')
          dt = Time(2)-Time(1);
          Elevation = new_gps.elev + Elevation_Correction*dt*c/2;
        else
          Elevation = new_gps.elev;
        end
        Roll = new_gps.roll;
        Pitch = new_gps.pitch;
        Heading = new_gps.heading;
        
        if param.debug_level > 1
          [x,y,z] = geodetic2ecef(Latitude/180*pi,Longitude/180*pi,Elevation,WGS84.ellipsoid);
          [tmp.x,tmp.y,tmp.z] = geodetic2ecef(tmp.Latitude/180*pi,tmp.Longitude/180*pi,tmp.Elevation,WGS84.ellipsoid);
          figure(1); clf;
          plot(sqrt( (x-tmp.x).^2 + (y-tmp.y).^2 + (z-tmp.z).^2 ));
          figure(2); clf;
          plot(Latitude);
          hold on;
          plot(tmp.Latitude,'r');
          hold off;
          figure(3); clf;
          plot(Longitude);
          hold on;
          plot(tmp.Longitude,'r');
          hold off;
          figure(4); clf;
          plot(Elevation);
          hold on;
          plot(tmp.Elevation,'r');
          hold off;
          figure(5); clf;
          imagesc(lp(tmp.Data));
          keyboard
        end
        if param.save_en
          if qlook_based % Quick look (qlook) or SAR file
            param_qlook = new_param;
            save(data_fns{data_idx},'GPS_time','Latitude','Longitude','Elevation','Roll','Pitch','Heading','param_records','param_qlook','-APPEND');
          else
            param_sar = new_param;
            save(data_fns{data_idx},'GPS_time','Latitude','Longitude','Elevation','Roll','Pitch','Heading','param_records','param_sar','-APPEND');
          end
        else
          fprintf('      Not saving information (TEST MODE)\n');
        end
        
      end
    end
  end
  
elseif strcmpi(param.update_mode,'layer')
  %% =====================================================================
  % Layer Update Section
  % ======================================================================
  error('This mode not supported');
  physical_constants;
  
  % Load all the layer data
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
  [GPS_time,sorting_idxs] = sort(GPS_time);
  surface = surface(sorting_idxs);
  bottom = bottom(sorting_idxs);
  elev = elev(sorting_idxs);
  
  reset_val_param_array_wf_chan = [];
  for file_idx = 1:length(fns)
    fn = fns{file_idx};
    
    % =======================================================================
    % Select only the frames files that are in the frames list
    % =======================================================================
    [fn_dir fn_name] = fileparts(fn);
    if fn_name(6) ~= 'i'
      % fn_name Format: Data_20110507_02_001
      frm_id = fn_name(6:end);
      frm = str2double(frm_id(end-2:end));
    else
      % fn_name Format: Data_img_01_20110507_02_001
      frm_id = fn_name(13:end);
      frm = str2double(frm_id(end-2:end));
    end
    if ~isempty(param.frms) && ~any(frm == param.frms)
      continue;
    end
    fprintf('Loading data file %s\n', fn);
    
    if param.gps.update_en
      param_array_wf_chan = reset_val_param_array_wf_chan;
      clear param_records;
      clear param_vectors;
      clear Bottom;
      load(fn);
      update_data = false;
      if exist('param_records','var')
        % Use records
        records_based = true;
      else
        % Use vectors
        records_based = false;
      end
      if exist('Bottom','var')
        update_bottom = true;
      else
        update_bottom = false;
      end
      
      if 1
        if records_based
          fprintf('  Old time %.2f sec to new time %.2f\n', ...
            param_records.gps.time_offset, new_time_offset);
          delta_offset = new_time_offset - param_records.gps.time_offset;
          param_records.gps.time_offset = new_time_offset;
        else
          fprintf('  Old time %.2f sec to new time %.2f\n', ...
            param_records.gps.time_offset, new_time_offset);
          delta_offset = new_time_offset - param_records.gps.time_offset;
          param_records.gps.time_offset = new_time_offset;
        end
        
        GPS_time = GPS_time + delta_offset;
        Latitude = interp1(gps.gps_time, gps.lat, GPS_time);
        Longitude = interp1(gps.gps_time, gps.lon, GPS_time);
        if ~exist('param_array_wf_chan','var') || isempty(param_array_wf_chan)
          fprintf('Old or qlook data files, does not contain elevation compensation info\n');
          fprintf('Qlook data files do not have elevation compensation\n');
          figure;
          plot(Elevation);
          title('Elevation Compensation Applied???');
          fprintf('Look at elevation plot\n');
          fprintf('Set elevation comp variables to correct values and then dbcont:\n');
          fprintf('  reset_val_param_array_wf_chan.elev_comp.en = 1\n');
          fprintf('  reset_val_param_array_wf_chan.elev_comp.en = 0\n');
          keyboard
          param_array_wf_chan = reset_val_param_array_wf_chan;
        end
        if param_array_wf_chan.elev_comp.en
          % Re-compensate for elevations changes with this new information
          Wrong_Elevation = interp1(gps.gps_time, gps.elev, GPS_time - delta_offset);
          Correct_Elevation = interp1(gps.gps_time, gps.elev, GPS_time);
          
          % dz = Range bin sample spacing (this Depth axis is air dielectric)
          dz = Depth(2)-Depth(1);
          % dt = Range bin time spacing
          dt = Time(2)-Time(1);
          % dBins = Amount of compensation to apply to each range line
          dBins = round((Correct_Elevation - Wrong_Elevation)/dz);
          if any(dBins ~= 0)
            update_data = true;
            % Circular shift data to compensate for elevation changes (note
            % that shift has range bin resolution so shift is rounded... this
            % is done for speed)
            for rline = 1:length(Elevation)
              if dBins(rline) > 0
                % Move target closer to radar
                Data(1:end-dBins(rline),rline) = Data(1+dBins(rline):end,rline);
                Data(end-dBins(rline)+1:end,rline) = 0;
              elseif dBins(rline) < 0
                % Move target further from radar
                Data(1-dBins(rline):end,rline) = Data(1:end+dBins(rline),rline);
                Data(1:-dBins(rline),rline) = 0;
              end
            end
            Surface = Surface - dBins*dt;
            if update_bottom
              Bottom = Bottom - dBins*dt;
            end
          end
        else
          Elevation = interp1(gps.gps_time, gps.elev, GPS_time);
        end
        
      else
        % Special fix for 20110411
        %   This is needed when the data files have seconds of day instead of
        %   seconds since the Jan 1, 1970 epoch.
        tmp = load(fn);
        Latitude = double(interp1(UTC_sod,gps.lat,GPS_time));
        Longitude = double(mod(interp1(UTC_sod,unwrap(gps.lon/180*pi),GPS_time)*180/pi+180, 360)-180);
        Elevation = double(interp1(UTC_sod,gps.elev,GPS_time));
        GPS_time = interp1(UTC_sod,gps.gps_time,GPS_time);
      end
      
      if param.save_en
        fprintf('  Saving GPS info\n');
        if records_based
          save(fn,'GPS_time','Latitude','Longitude','Elevation','param_records','-APPEND');
        else
          save(fn,'GPS_time','Latitude','Longitude','Elevation','param_vectors','-APPEND');
        end
        if update_data
          fprintf('  Saving Data/Layer info\n');
          if update_bottom
            save(fn,'Data','Surface','Bottom','-APPEND');
          else
            save(fn,'Data','Surface','-APPEND');
          end
        end
      else
        fprintf('  Not saving information (TEST MODE)\n');
      end
    end
    
  end
end

return;
