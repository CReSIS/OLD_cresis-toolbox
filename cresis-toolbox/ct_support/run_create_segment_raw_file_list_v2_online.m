% script run_create_segment_raw_file_list_v2_online.m
%
% Runs script run_create_segment_raw_file_list_v2_online.m
%
% Instructions:
% 1. Add in sections from run_create_segment_raw_file_list_v2
% 2. Enable GPS copying if needed and update GPS paths if necessary
% 3. Enable make GPS if needed, update script name and update make GPS script if necessary
% 4. Run the script
%
% Author: John Paden

while 1
  
  counter_correction_en = false;
  online_mode = true;
  
  %% User Settings that should not generally be changed
  % You may have to set to false to read some of the results from this function when it was first written (should always be true)
  tmp_fn_uses_adc_folder_name = true;
  
  MIN_SEG_SIZE = 2;
  MAX_TIME_GAP = 1000/75;
  MIN_PRF = 100;
  
  %% Ku-band 3
  if 1
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
    day_string = '20161015'; % Only used for stdout print of the vectors worksheet
    
    create_segment_raw_file_list_v2;
  end
  
  %% RDS: MCoRDS 3
  if 1
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
    day_string = '20161015'; % Only used for stdout print of the vectors worksheet
    
    create_segment_raw_file_list_v2;
  end
  
  %% Snow 3 (OIB)
  if 1
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
    file_midfix = '20161015'; % Data files must contain this string in the middle of their name (usually should be empty)
    day_string = '20161015'; % Only used for stdout print of the vectors worksheet
    
    create_segment_raw_file_list_v2;
  end
  
  %% Copy GPS data
  if 1
    mkdir('/scratch/metadata/2016_Antarctica_DC8/20161015/');
    copyfile('/net/field1/landing/mcords/GPS*','/scratch/metadata/2016_Antarctica_DC8/20161015/')
  end
  
  if 1
    %% Make GPS
    make_gps_2016_antarctica_DC8;
  end
  
  fprintf('Pausing for 30 seconds...\n');
  pause(30);
end

return;
