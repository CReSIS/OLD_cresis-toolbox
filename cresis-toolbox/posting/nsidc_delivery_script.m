% Scripts nsidc_delivery_script
%
% Delivery L1B and L2 data to nsidc
% For L1B
%        Input: .mat files
%        Output: .premet, .spatial and netcdf files
% For L2
%        INput: .csv file
%        Output: .premet, .spatial and ascii files
%
% See also: type "nsidc_help.m"
%
% Author: Yi Zhu, John Paden

%% User Settings
% User defined directory  
USER_SPECIFIED_DIRECTORY_BASE = '/cresis/snfs1/scratch/jliwestc/nsidc/';

% Hardcoded local version ID (our local version number, check with NSIDC
% before changing if re-sending data)
premet_param.version_id = '002';  

% MCF version ID (only change if NSIDC requests change)
mcf_version_id_L1B = '002';
mcf_version_id_L2 = '001';

% Read rds spreadsheet
params_fn = ct_filename_param('rds_param_2014_Antarctica_DC8.xls');
% params_fn = ct_filename_param('rds_param_2015_Greenland_L130.xls');
skip_phrase = 'do not process';

% Post L1B data
L1B_cmd = false;
% data_dir_L1 = fullfile('CSARP_post','qlook'); % Non-RDS
data_dir_L1 = fullfile('CSARP_post','mvdr'); % RDS

% Post L2 data
L2_cmd = true;
data_dir_L2 = fullfile('CSARP_post','csv');

if L2_cmd
  premet_param.version_id = '001';
end

%% Automated Section
fprintf('===============================================\n');
fprintf('NSIDC Delivery Script\n\n');

% Control parameters for netcdf file
netcdf_param = [];
netcdf_param(end+1).mat_name = 'GPS_time';
netcdf_param(end).cdf_name = 'time';
netcdf_param(end).dim_name = {'time'};
% netcdf_param(end).eval = @(x) (epoch_to_sod(x - utc_leap_seconds(x(1)),param.day_seg(1:8)));
% This is done later
netcdf_param(end).attributes = {'units' 'TO BE FILLED IN LATER', ...
  'calendar', 'standard', ...
  'long_name', 'Time of day UTC', ...
  'standard_name', 'time', ...
  'axis', 'T'};

netcdf_param(end+1).mat_name = 'Time';
netcdf_param(end).cdf_name = 'fasttime';
netcdf_param(end).dim_name = {'fasttime_dim'};
netcdf_param(end).eval = @(x) (x*1e6);
netcdf_param(end).attributes = {'units' 'microseconds', ...
  'long_name', '2-way travel time in useconds', ...
  'standard_name', 'time', ...
  'positive', 'down', ...
  'axis', 'Z'};

netcdf_param(end+1).mat_name = 'Elevation';
netcdf_param(end).cdf_name = 'altitude';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'standard_name','height', ...
  'units', 'meters', ...
  'positive', 'up', ...
  'long_name', 'Altitude of antenna above nominal sea level (WGS84)'};

netcdf_param(end+1).mat_name = 'Latitude';
netcdf_param(end).cdf_name = 'lat';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'standard_name','latitude', ...
  'units', 'degrees_north', ...
  'long_name', 'Latitude of sample', ...
  'axis', 'Y'};

netcdf_param(end+1).mat_name = 'Longitude';
netcdf_param(end).cdf_name = 'lon';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'standard_name','longitude', ...
  'units', 'degrees_east', ...
  'long_name', 'Longitude of sample', ...
  'axis', 'X'};

netcdf_param(end+1).mat_name = 'Data';
netcdf_param(end).cdf_name = 'amplitude';
netcdf_param(end).dim_name = {'fasttime_dim' 'time'};
netcdf_param(end).eval = @(x) (10*log10(x));
netcdf_param(end).attributes = {'units', 'counts in dB', ...
  'long_name', 'Amplitude of low/high gain merged radar reflection after processing', ...
  'coordinates', 'time travel channel'};

netcdf_param(end+1).mat_name = 'Heading';
netcdf_param(end).cdf_name = 'heading';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).eval = @(x) (x*180/pi);
netcdf_param(end).attributes = {'standard_name','heading', ...
  'units', 'degrees', ...
  'long_name', 'Heading of the platform.  Positive is clockwise from above.  Zero is true north.'};

netcdf_param(end+1).mat_name = 'Pitch';
netcdf_param(end).cdf_name = 'pitch';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).eval = @(x) (x*180/pi);
netcdf_param(end).attributes = {'standard_name','pitch', ...
  'units', 'degrees', ...
  'long_name', 'Pitch of the platform.  Positive is nose up.  Zero is horizontal.'};

netcdf_param(end+1).mat_name = 'Roll';
netcdf_param(end).cdf_name = 'roll';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).eval = @(x) (x*180/pi);
netcdf_param(end).attributes = {'standard_name','roll', ...
  'units', 'degrees', ...
  'long_name', 'Roll of the platform.  Positive is right wing down.  Zero is horizontal.'};

netcdf_param(end+1).mat_name = 'Surface';
netcdf_param(end).cdf_name = 'Surface';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'seconds', ...
  'long_name', 'Two way travel time to surface used during processing. Not the final picked surface.'};

netcdf_param(end+1).mat_name = 'Bottom';
netcdf_param(end).cdf_name = 'Bottom';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'seconds', ...
  'long_name', 'Two way travel time to bottom used during processing. Not the final picked bottom.'};

netcdf_param(end+1).mat_name = 'Elevation_Correction';
netcdf_param(end).cdf_name = 'Elevation_Correction';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'range bins', ...
  'long_name', 'Represents the number of zeros that were inserted during elevation compensation for each range line to simulate near-level flight. These zeros are not included in the truncation noise statistics.'};

netcdf_param(end+1).mat_name = 'Truncate_Bins';
netcdf_param(end).cdf_name = 'Truncate_Bins';
netcdf_param(end).dim_name = {'fasttime_dim'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Indices into the original fasttime vector for which the data "amplitude" are available.'};

netcdf_param(end+1).mat_name = 'Truncate_Mean';
netcdf_param(end).cdf_name = 'Truncate_Mean';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Represents a mean of the noise power for the truncated range bins before the surface return. When no range bins were truncated before the surface return the value is NaN.'};

netcdf_param(end+1).mat_name = 'Truncate_Median';
netcdf_param(end).cdf_name = 'Truncate_Median';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Represents a median of the noise power for the truncated range bins before the surface return. When no range bins were truncated before the surface return the value is NaN.'};

netcdf_param(end+1).mat_name = 'Truncate_Std_Dev';
netcdf_param(end).cdf_name = 'Truncate_Std_Dev';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Represents a standard deviation of the noise power for the truncated range bins before the surface return. When no range bins were truncated before the surface return the value is NaN.'};


%% Main loop for each segment
params = read_param_xls(params_fn);
for param_idx = 1:length(params)
  
  param = params(param_idx);
  if ~isempty(skip_phrase) ...
      && ~isempty(strfind(lower(param.cmd.notes),skip_phrase)) ...
      || ~params(param_idx).cmd.generic
    continue;
  end
  
  USER_SPECIFIED_DIRECTORY = fullfile(USER_SPECIFIED_DIRECTORY_BASE, ...
    param.season_name,ct_output_dir(param.radar_name));
  if ~exist(USER_SPECIFIED_DIRECTORY,'dir')
    mkdir(USER_SPECIFIED_DIRECTORY);
  end

  fprintf('NSIDC prep %s\n', param.day_seg);
  
  % Extract year(2013), location(GR) and platform(P3) information from
  % seanson name in param
  [season_name_year location] = strtok(param.season_name,'_');
  [location platform] = strtok(location,'_');
  platform = platform(2:end);
  
  if strcmpi(location,'Greenland')
    location = 'GR';
  elseif strcmpi(location,'Antarctica')
    location = 'AN';
  else
    error('Unsupported location %s\n', location);
  end
  
  % Construct the necessary element: AircraftID and platform_short_name
  if strcmpi(platform,'P3')
    premet_param.nsidc_platform_short_name = 'P-3B';        % change to meet the valids file
    premet_param.nsidc_aircraft_id = 'N426NA';
  elseif strcmpi(platform,'DC8')
    premet_param.nsidc_platform_short_name = 'DC-8';
    premet_param.nsidc_aircraft_id = 'N817NA';
  else
    error('Unsupported platform %s\n', platform);
  end
  
  % Construct the necessary element: ThemeID
  premet_param.nsidc_season_name = sprintf('%s_%s_NASA', season_name_year, location);
  data_files_dir = sprintf('Data_%s_%s', season_name_year, location);
  % Construct the necessary element: instrument_short_name and sensor_short_name
  radar_sensor = ct_output_dir(param.radar_name);
  if strcmpi(radar_sensor,'accum')
    premet_param.nsidc_instrument_short_name = 'Accumulation Radar';
    premet_param.nsidc_sensor_short_name = 'Accumulation Radar';
    radar_type = 'IRACC';
  elseif strcmpi(radar_sensor,'kuband')
    premet_param.nsidc_instrument_short_name = 'Ku-Band Radar';
    premet_param.nsidc_sensor_short_name = 'Ku-Band Radar';
    radar_type = 'IRKUB';
  elseif strcmpi(radar_sensor,'rds')
    premet_param.nsidc_instrument_short_name = 'MCoRDS';
    premet_param.nsidc_sensor_short_name = 'MCoRDS';
    radar_type = 'IRMCR';
  elseif strcmpi(radar_sensor,'snow')
    premet_param.nsidc_instrument_short_name = 'Snow Radar';
    premet_param.nsidc_sensor_short_name = 'Snow Radar';
    radar_type = 'IRSNO';
  else
    error('Unsupported sensor %s\n', radar_sensor);
  end
  
  %% For L2 data
  if L2_cmd

    data_type = '2';
    
    % Input L2 .csv file (The path here should be self-defined)
    csv_fn = fullfile(ct_filename_out(param,'',data_dir_L2,2), ...
      sprintf('Data_%s.csv',param.day_seg));
    
    % Output filename for .premet, e.g. IRMCRL2_20130402_01.csv.premet
    out_fn_premet_L2 = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
      sprintf('%s%s_Files', radar_type, data_type),'Spatial_Premet_Files', ...
      sprintf('%s%s_%s.csv.premet',radar_type,data_type,param.day_seg));
    pre_spatial_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
      sprintf('%s%s_Files', radar_type, data_type),'Spatial_Premet_Files');
    if ~exist(pre_spatial_dir,'dir')
      mkdir(pre_spatial_dir);
    end
    
    % Output filename for .spatial, e.g. IRMCRL2_20130402_01.csv.spatial
    out_fn_spatial_L2 = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
      sprintf('%s%s_Files', radar_type, data_type),'Spatial_Premet_Files', ...
      sprintf('%s%s_%s.csv.spatial',radar_type,data_type,param.day_seg));

    % Output txt/csv directory
    data_txt_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
      sprintf('%s%s_Files', radar_type, data_type),data_files_dir);
    if ~exist(data_txt_dir,'dir')
      mkdir(data_txt_dir);
    end
    
    % Create output pdr and met files directories
    out_pdr_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
      sprintf('%s%s_Files', radar_type, data_type),'output','pdrs');
    if ~exist(out_pdr_dir,'dir')
      mkdir(out_pdr_dir);
    end   
    out_met_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
      sprintf('%s%s_Files', radar_type, data_type),'output','mets');
    if ~exist(out_met_dir,'dir')
      mkdir(out_met_dir);
    end
    
    new_csv_fn = fullfile(data_txt_dir, ...
      sprintf('%s%s_%s.csv',radar_type,data_type,param.day_seg));
    
    fprintf('  Original csv: %s\n', csv_fn);
    fprintf('  premet: %s\n', out_fn_premet_L2);
    fprintf('  spatial: %s\n', out_fn_spatial_L2);
    fprintf('  New csv %s\n', new_csv_fn);
    
    if ~exist(csv_fn)
      warning('Skipping: csv file %s is missing\n',csv_fn);
      keyboard
      continue;
    end
    
    % Copy csv file
    [flag,message,~] = copyfile(csv_fn,new_csv_fn);
      
    % Create the data filename without extension
    premet_param.data_fn_name = sprintf('%s_%s', ...
      sprintf('%s%s', radar_type, data_type), param.day_seg);
    
    % Create .premet, .spatial and .txt file for L2 data
    nsidc_create_premet_L2(csv_fn,out_fn_premet_L2,premet_param);
    nsidc_create_spatial_L2(csv_fn,out_fn_spatial_L2);
    
    % Create MetGen config file (this actually only needs to be run once for the whole season)
    nsidc_change_metgen(USER_SPECIFIED_DIRECTORY, mcf_version_id_L2, ...
      pre_spatial_dir, data_txt_dir, radar_type, data_type);
    
  end
  
  %% For L1B data 
  if L1B_cmd
    
    data_type = '1B';
    
    % Find the frames information
    frames_fn = ct_filename_support(param,'','frames');
    load(frames_fn);
    if isempty(param.cmd.frms)
      param.cmd.frms = 1:length(frames.frame_idxs);
    end
    
    % Find the records information for updating global attributes
    records_fn = ct_filename_support(param,'','records');
    records = load(records_fn);
    
    % Remove frames that do not exist from param.cmd.frms list
    [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
    if length(valid_frms) ~= length(param.cmd.frms)
      bad_mask = ones(size(param.cmd.frms));
      bad_mask(keep_idxs) = 0;
      warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
        param.cmd.frms(find(bad_mask,1)));
      param.cmd.frms = valid_frms;
    end
    
    % Create premet, spatial, netcdf files for each frame
    param.records = records;
    for frm_idx = 1:length(param.cmd.frms)
      frm = param.cmd.frms(frm_idx);
      
      fprintf('  Converting frame %d (%d of %d) (%s)\n', frm, frm_idx, ...
        length(param.cmd.frms), datestr(now,'HH:MM:SS'));
      
      % Input L1B .mat file (The path here should be self-defined)
      echo_fn = fullfile(ct_filename_out(param,data_dir_L1), ...
        sprintf('Data_%s_%03d.mat',param.day_seg,frm));
      
      % Output filename for .premet, e.g. IRMCR1B_V01_20130402_01_001.premet
      out_fn_premet_L1B = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
        sprintf('%s%s_Files', radar_type, data_type),'Spatial_Premet_Files', ...
        sprintf('%s%s_%s_%03d.premet',radar_type, data_type, param.day_seg, frm));
      pre_spatial_dir = fileparts(out_fn_premet_L1B);
      if ~exist(pre_spatial_dir,'dir')
        mkdir(pre_spatial_dir);
      end
      
      % Output filename for .spatial, e.g. IRMCR1B_V01_20130402_01_001.spatial
      out_fn_spatial_L1B = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
        sprintf('%s%s_Files', radar_type, data_type),'Spatial_Premet_Files', ...
        sprintf('%s%s_%s_%03d.spatial',radar_type, data_type, param.day_seg, frm));
      
      % Output filename for .nc, e.g. IRMCR1B_V01_20130402_01_001.nc
      out_fn_netcdf = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
        sprintf('%s%s_Files', radar_type, data_type),data_files_dir, ...
        sprintf('%s%s_%s_%03d.nc',radar_type, data_type, param.day_seg, frm));
      out_data_dir = fileparts(out_fn_netcdf);
      if ~exist(out_data_dir,'dir')
        mkdir(out_data_dir);
      end
      
      % Create output pdr and met files directories
      out_pdr_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
        sprintf('%s%s_Files', radar_type, data_type),'output','pdrs');
      if ~exist(out_pdr_dir,'dir')
        mkdir(out_pdr_dir);
      end
      out_met_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
        sprintf('%s%s_Files', radar_type, data_type),'output','mets');
      if ~exist(out_met_dir,'dir')
        mkdir(out_met_dir);
      end
      
      % Create Maps, Echogram and Picked Echogram image filenames
      map_fn = fullfile(ct_filename_out(param,'post','',1),'images', ...
        param.day_seg, sprintf('%s_%03d_0maps.jpg',param.day_seg,frm));
      new_map_fn = fullfile(out_data_dir, ...
        sprintf('%s%s_%s_%03d_%s.jpg',radar_type,data_type,param.day_seg,frm,'Map'));

      echogram_fn = fullfile(ct_filename_out(param,'post','',1),'images', ...
        param.day_seg, sprintf('%s_%03d_1echo.jpg',param.day_seg,frm));
      new_echogram_fn = fullfile(out_data_dir, ...
        sprintf('%s%s_%s_%03d_%s.jpg',radar_type,data_type,param.day_seg,frm,'Echogram'));
      
      echogram_picks_fn = fullfile(ct_filename_out(param,'post','',1),'images', ...
        param.day_seg, sprintf('%s_%03d_2echo_picks.jpg',param.day_seg,frm));
      new_echogram_picks_fn = fullfile(out_data_dir, ...
        sprintf('%s%s_%s_%03d_%s.jpg',radar_type,data_type,param.day_seg,frm,'Echogram_Picks'));
      
      % Check for the existence of image files
      if ~exist(map_fn,'file')
        warning('Skipping: missing map file %s', map_fn);
        keyboard;
        continue;
      end
      if ~exist(echogram_fn,'file')
        warning('Skipping: missing map file %s', echogram_fn);
        keyboard;
        continue;
      end
      if strcmpi(radar_type,'IRMCR')
        if ~exist(echogram_picks_fn,'file')
          warning('Skipping: missing pick file %s', echogram_fn);
          keyboard;
          continue;
        end
      end
      
      fprintf('   input mat: %s\n', echo_fn);
      fprintf('   premet: %s\n', out_fn_premet_L1B);
      fprintf('   spatial: %s\n', out_fn_spatial_L1B);
      fprintf('   nc: %s\n', out_fn_netcdf);
      fprintf('   NetCDF and images are copying to %s\n', out_data_dir);
      
      % Copy image files
      [flag,message,~] = copyfile(map_fn,new_map_fn);
      [flag,message,~] = copyfile(echogram_fn,new_echogram_fn);
      if strcmpi(radar_type,'IRMCR')
        [flag,message,~] = copyfile(echogram_picks_fn,new_echogram_picks_fn);
      end
      
      % Create .premet, .spatial files for L1B data
      nsidc_create_premet_L1B(echo_fn,out_fn_premet_L1B,premet_param);
      nsidc_create_spatial_L1B(echo_fn, out_fn_spatial_L1B);
      
      %% Create .nc files
      mdata = load(echo_fn);
      if isfield(mdata,'Truncate_Bins')
        % Special modification for compressed echograms
        mdata.Time = mdata.Time(mdata.Truncate_Bins);
      end

      % Netcdf time is seconds of day, so we need to indicate which day it
      % is taken relative to
      netcdf_param(1).attributes{2} = sprintf('seconds since %s-%s-%s 00:00:00', ...
        param.day_seg(1:4), param.day_seg(5:6), param.day_seg(7:8));
      netcdf_param(find(strcmp('GPS_time',{netcdf_param.mat_name}))).eval = @(x) (epoch_to_sod(x - utc_leap_seconds(x(1)),param.day_seg)); 
      netcdf_from_mat(out_fn_netcdf,mdata,netcdf_param);

      % Add global attributes
      nsidc_netcdf_attributes(out_fn_netcdf,echo_fn,param);
      
      % Create MetGen config file (this actually only needs to be run once for the whole season)
      nsidc_change_metgen(USER_SPECIFIED_DIRECTORY, mcf_version_id_L1B, pre_spatial_dir, out_data_dir, radar_type, data_type);
      
    end
  end
end

return;
