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
  online_mode = true;
  day_string = '20161020';
  kuband_en = 1;
  rds_en = 1;
  snow_en = 1;
  gps_en = 1; % Remember to add a new entry into make_gps_*.m file that you are using
  
  
  %% User Settings that should not generally be changed
  % You may have to set to false to read some of the results from this function when it was first written (should always be true)
  tmp_fn_uses_adc_folder_name = true;
  
  MIN_SEG_SIZE = 2;
  MAX_TIME_GAP = 1000/75;
  MIN_PRF = 100;
  
  %% Ku-band 3
  if kuband_en
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
    
    create_segment_raw_file_list_v2;
  end
  
  %% RDS: MCoRDS 3
  if rds_en
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
    param.season_name = '2016_Antarctica_DC8';
    base_dir = '/process/mcords/';
    %   base_dir = '/net/field1/landing/mcords/';
    
    param.adc_folder_name = 'board%b';
    file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
    
    create_segment_raw_file_list_v2;
  end
  
  %% Snow 3 (OIB)
  if snow_en
    param = [];
    param.radar_name = 'snow3';
    param.clk = 125e6;
    adcs = 1;
    param.file_version = 5; % 3 for 2013 Gr, 5 for after that
    raw_file_suffix = '.bin';
    reuse_tmp_files = true; % Set to false if you want to overwrite current results
    file_prefix_override = '';
    counter_correction_en = true;
    
    % Parameters below this point OFTEN NEEDS TO BE CHANGED
    param.season_name = '2016_Antarctica_DC8';
    base_dir = '/process/fmcw/snow/';
    param.adc_folder_name = '';
    file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
    
    create_segment_raw_file_list_v2;
  end
  
  %% Copy GPS data
  if gps_en
    mkdir(sprintf('/scratch/metadata/2016_Antarctica_DC8/%s/',day_string));
    copyfile('/net/field1/landing/mcords/GPS*',sprintf('/scratch/metadata/2016_Antarctica_DC8/%s/',day_string))
  end
  
  if gps_en
    %% Make GPS
    make_gps_2016_antarctica_DC8;
  end
  
  fprintf('Pausing for 30 seconds...\n');
  pause(30);
end

return;
