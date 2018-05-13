% script run_create_segment_raw_file_list_v2_online.m
%
% Runs script run_create_segment_raw_file_list_v2_online.m
%
% Instructions:
% 1. Add in sections from run_create_segment_raw_file_list_v2
%    If already done, just ensure enables are setup correctly.
% 2. Enable GPS copying and make GPS if needed
%    Update GPS paths if necessary
%    Update make_gps_* script name and update the script if necessary
% 4. Run the script
%
% Author: John Paden

while 1
  counter_correction_en = false;
  online_mode = 2; % Set to 1 to read all files, 2 to read only last file
  pause_time = 180;
  day_string = '20180512';
  mcords5_accum_en = 0;
  kuband3_en = 0;
  mcords3_en = 0;
  uwb_en = 1;
  snow3_en = 0;
  gps_en = 0; % Remember to add a new entry into make_gps_*.m file that you are using  
  
  %% User Settings that should not generally be changed
  % You may have to set to false to read some of the results from this function when it was first written (should always be true)
  tmp_fn_uses_adc_folder_name = true;
  
  MIN_SEG_SIZE = 2;
  MAX_TIME_GAP = 1000/75;
  MIN_PRF = 100;
  
  %% Accumulation Radar (mcords5-accum)
  if mcords5_accum_en
    param = [];
    param.radar_name = 'mcords5';
    param.clk = 1.6e9/8;
    adcs = 1:1;
    raw_file_suffix = '.bin';
    reuse_tmp_files = true; % Set to false if you want to overwrite current results
    file_prefix_override = ''; % most of the time
    counter_correction_en = true;
    presum_bug_fixed = false; % Seasons from 2015 Greenland Polar6 onward should be set to true, except for 2017 Antarctica Basler which uses the cresis DDS with this bug
    union_time_epri_gaps = true;
    
    % Parameters below this point OFTEN NEEDS TO BE CHANGED
    param.season_name = '2018_Greenland_P3';
    base_dir = '/process/accum/';
    param.adc_folder_name = 'chan%d';
    file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
    day_string = '20180403'; % Only used for stdout print of the vectors worksheet
        
    try
      create_segment_raw_file_list_v2;
    catch ME
      warning(ME.getReport())
    end
  end
  
  %% Ku-band 3
  if kuband3_en
    param = [];
    param.radar_name = 'kuband3';
    param.clk = 125e6;
    adcs = 1;
    param.file_version = 5; % 3 for 2013 Gr, 5 for after that
    raw_file_suffix = '.bin';
    reuse_tmp_files = true; % Set to false if you want to overwrite current results
    file_prefix_override = '';
    counter_correction_en = true;
    
    % Parameters below this point OFTEN NEEDS TO BE CHANGED
    param.season_name = '2016_Antarctica_DC8';
    base_dir = '/process/fmcw/kuband/';
    param.adc_folder_name = '';
    file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
    
    try
      create_segment_raw_file_list_v2;
    catch ME
      warning(ME.getReport())
    end
  end
  
  %% RDS: MCoRDS 3
  if mcords3_en
    param = [];
    param.radar_name = 'mcords3';
    param.clk = 150e6;
    adcs = [1 5];
    raw_file_suffix = '.bin';
    reuse_tmp_files = true; % Set to false if you want to overwrite current results
    file_prefix_override = ''; % most of the time
    counter_correction_en = true;
    presum_bug_fixed = false;
    union_time_epri_gaps = true;
    
    % Parameters below this point OFTEN NEEDS TO BE CHANGED
    param.season_name = '2018_Greenland_P3';
    base_dir = '/process/mcords/';
    %   base_dir = '/net/field1/landing/mcords/';
    
    param.adc_folder_name = 'board%b';
    file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
    
    try
      create_segment_raw_file_list_v2;
    catch ME
      warning(ME.getReport())
    end
  end
  
  %% RDS: MCoRDS 5 UWB
  if uwb_en
    param = [];
    
    param.radar_name = 'mcords5';
    param.clk = 1.6e9/8;
    adcs = 1; % Choose a smaller subset of channels to check if data rate an issue
    raw_file_suffix = '.bin';
    reuse_tmp_files = true; % Set to false if you want to overwrite current results
    file_prefix_override = ''; % most of the time
    counter_correction_en = true;
    presum_bug_fixed = false; % Seasons from 2015 Greenland Polar6 onward should be set to true, except for 2017 Antarctica Basler which uses the cresis DDS with this bug
    union_time_epri_gaps = true;
    
    % Parameters below this point OFTEN NEEDS TO BE CHANGED
    param.season_name = '2018_Greenland_Polar6';
    base_dir = 'D:\awi\';
    param.adc_folder_name = 'chan%d';
    file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  
    try
      create_segment_raw_file_list_v2;
    catch ME
      warning(ME.getReport())
    end
  end
  
  %% Snow 8 (OIB)
  if snow8_en
    param = [];
    param.radar_name = 'snow8';
    param.clk = 125e6;
    adcs = 1;
    param.file_version = 8; % 3 for 2013 Gr, 5 for after that
    raw_file_suffix = '.bin';
    reuse_tmp_files = true; % Set to false if you want to overwrite current results
    file_prefix_override = '';
    counter_correction_en = true;
    
    % Parameters below this point OFTEN NEEDS TO BE CHANGED
    param.season_name = '2018_Greenland_P3';
    base_dir = '/process/fmcw/snow/';
    param.adc_folder_name = '';
    file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
    
    try
      create_segment_raw_file_list_v2;
    catch ME
      warning(ME.getReport())
    end
  end
  
  %% Copy GPS data
  if gps_en
    gps_dir = sprintf('/scratch/metadata/2018_Greenland_P3/%s/',day_string);
    if ~exist(gps_dir,'dir')
      mkdir(gps_dir);
    end
   try
      copyfile('/net/field1/landing/mcords/GPS*',sprintf('/scratch/metadata/2018_Greenland_P3/%s/',day_string))
    catch ME
      warning(ME.getReport);
    end
  end
  
  if gps_en
    %% Make GPS
    try
      make_gps_2018_greenland_P3;
    catch ME
      warning(ME.getReport);
    end
  end
  
  fprintf('Pausing for %.0f seconds...\n', pause_time);
  pause(pause_time);
end

return;
