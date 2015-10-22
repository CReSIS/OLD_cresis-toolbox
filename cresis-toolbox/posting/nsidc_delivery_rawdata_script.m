% Scripts nsidc_delivery_rawdata_script
%
% Delivery raw data of snow radar to nsidc
% For L0
%        Input: .dat or .bin files
%        Output: .premet, .spatial and netcdf files
%
% See also: type "nsidc_help.m"
%
% Author: Yi Zhu, John Paden, Jilu Li

%% User Settings
% User defined directory  
USER_SPECIFIED_DIRECTORY_BASE = '/cresis/snfs1/scratch/jliwestc/nsidc/';

% Hardcoded local version ID (our local version number, check with NSIDC
% before changing if re-sending data)
premet_param.version_id = '001';  

% MCF version ID (currently have 001~006 correspoinding to 6 file versions for snow radar)
mcf_version_id_L0 = '001';  % get overwirtten later according to record spreadsheet

% Read rds spreadsheet
params_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls');
skip_phrase = 'do not process';

% Post L0 data (raw data)
L0_cmd = true;
if L0_cmd
  data_type = '0';
end
data_dir_L0 = '/cresis/snfs1/data/SnowRadar/2014_Greenland_P3/'; 

%% Automated Section
fprintf('===============================================\n');
fprintf('NSIDC Delivery Script\n\n');

%% Extract information from param spreadsheets
params = read_param_xls(params_fn);
param = params(1);
USER_SPECIFIED_DIRECTORY = fullfile(USER_SPECIFIED_DIRECTORY_BASE, ...
  param.season_name,ct_output_dir(param.radar_name));
if ~exist(USER_SPECIFIED_DIRECTORY,'dir')
  mkdir(USER_SPECIFIED_DIRECTORY);
end
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
elseif strcmpi(platform,'C130')
  premet_param.nsidc_platform_short_name = 'C-130';
  premet_param.nsidc_aircraft_id = 'N436NA';
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

% Create directory for .premet and .spatial files
pre_spatial_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
        sprintf('%s%s_Files', radar_type, data_type),'Spatial_Premet_Files');
if ~exist(pre_spatial_dir,'dir')
  mkdir(pre_spatial_dir);
end

% Create directory for raw data files
out_data_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
  sprintf('%s%s_Files', radar_type, data_type),data_files_dir);
if ~exist(out_data_dir,'dir')
  mkdir(out_data_dir);
end

% Create directory for output pdr files 
out_pdr_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
  sprintf('%s%s_Files', radar_type, data_type),'output','pdrs');
if ~exist(out_pdr_dir,'dir')
  mkdir(out_pdr_dir);
end

% Create directory for output met files 
out_met_dir = fullfile(ct_filename_out(param,USER_SPECIFIED_DIRECTORY,'',1), ...
  sprintf('%s%s_Files', radar_type, data_type),'output','mets');
if ~exist(out_met_dir,'dir')
  mkdir(out_met_dir);
end

%% Main loop for each segment
for param_idx = 1:length(params)  
  param = params(param_idx);
  if ~isempty(skip_phrase) ...
      && ~isempty(strfind(lower(param.cmd.notes),skip_phrase)) ...
      || ~params(param_idx).cmd.generic
    continue;
  end 
  fprintf('NSIDC prep %s\n', param.day_seg);
  mcf_version_id_L0 = sprintf('%03d',param.records.file_version');
  
  %% For L0 data 
  if L0_cmd
    fname_prefix = sprintf('%s%s_%02d_', radar_type, data_type,param.records.file_version);
    
    % Get raw data files    
    [base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param);
        
    % Create premet, spatial files for each raw data file
    for file_idx = 1:length(file_idxs)
      raw_fn = fns{file_idxs(file_idx)};
      fname = fname_info_fmcw(raw_fn); 
      raw_fn0 = fname.dir_info.name;
%       fprintf('  File index %d/%d (filename idx %d) (%s)\n', ...
%         file_idx, length(file_idxs), fname.file_idx, datestr(now,'HH:MM:SS'));
      if strcmpi(raw_fn0(1:4),'data')
        fname_midfix = sprintf('%s_%s%s_%s',raw_fn0(5:6),raw_fn0(12:15),raw_fn0(8:11),raw_fn0(end-7:end-4));
        fname_ext = raw_fn0(end-3:end);
      elseif strcmpi(raw_fn0(1:5),'snow0')
        fname_midfix = sprintf('%s_%s_%s',raw_fn0(5:6),raw_fn0(8:15),raw_fn0(end-7:end-4));
        fname_ext = raw_fn0(end-3:end);
      elseif strcmpi(raw_fn0(1:5),'snow_')
        fname_midfix = sprintf('%s_%s_%s',raw_fn0(6:7),raw_fn0(9:16),raw_fn0(end-7:end-4));
        fname_ext = raw_fn0(end-3:end);
      elseif strcmpi(raw_fn0(1:5),'snow3')
        fname_midfix = sprintf('%s_%s_%s',raw_fn0(7:8),raw_fn0(10:17),raw_fn0(end-7:end-4));
        fname_ext = raw_fn0(end-3:end);
      end
      new_raw_fn0 = sprintf('%s%s%s',fname_prefix,fname_midfix,fname_ext);
      new_raw_fn = fullfile(out_data_dir,new_raw_fn0);

      % Output filename for .premet, e.g. IRSNO0_00_20130402_001.premet
      out_fn_premet_L0 = fullfile(pre_spatial_dir,sprintf('%s%s%s.premet',fname_prefix,fname_midfix,fname_ext));
                                
      fprintf('   raw data file: %s\n', raw_fn);
      fprintf('   premet: %s\n', out_fn_premet_L0);
      fprintf('   raw data files are copying to %s\n', out_data_dir);
      
      % Copy raw data files
      [flag,message,~] = copyfile(raw_fn,new_raw_fn);
      
      % Create .premet, files for L0 data
      nsidc_create_premet_L0(raw_fn,out_fn_premet_L0,premet_param);
           
      % Create MetGen config file (this actually only needs to be run once for the whole season)
%       nsidc_change_metgen(USER_SPECIFIED_DIRECTORY, mcf_version_id_L0, pre_spatial_dir, out_data_dir, radar_type, data_type);
      
    end
  end
end

return;
