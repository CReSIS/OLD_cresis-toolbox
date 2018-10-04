% script run_preprocess.m
%
% Runs script preprocess.m
%
% Instructions:
% 1. Set your default radar parameters file
% 2. Set the input directory of the data files
% 3. Set the input directory of the config files if applicable
% 4. Run the script
%
% Author: John Paden

%% User Setup
% =========================================================================
param_override = [];

% Set param.radar_name and param.season_name and get radar default
% parameters.
param = default_radar_params_2018_Greenland_P3_snow;
% param = default_radar_params_2018_Antarctica_TObas_accum3;
% param = default_radar_params_2018_Antarctica_Ground_mcords6;

if ispc
  param.preprocess.base_dir = 'E:\tmp\2018_Antarctica_TObas\';
else
  param.preprocess.base_dir = '/N/dcwan/projects/cresis/2018_Greenland_P3/';
  % param.preprocess.base_dir = '/data/';
end
param.preprocess.config_folder_names = {'20180405/fmcw/snow/'};
param.preprocess.board_folder_names = {'20180405/fmcw/snow/'};
param.preprocess.date_strs = {'20180405'};
% param.preprocess.config_folder_names = {'20180929/'};
% param.preprocess.board_folder_names = {'20180929/%b'};


dbstop if error;
param_override.cluster.type = 'debug';
% param_override.cluster.type = 'matlab';
% param_override.cluster.type = 'debug';
%param_override.cluster.rerun_only = true;
param_override.cluster.desired_time_per_job  = 0*60;
% param_override.cluster.cpu_time_mult  = 2;
% param_override.cluster.mem_mult  = 2;

%% Automated Section
% =====================================================================

% Input checking
global gRadar;
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

ctrl_chain = preprocess(param,param_override);

cluster_print_chain(ctrl_chain);

[chain_fn,chain_id] = cluster_save_chain(ctrl_chain);

% Potentially stop and inspect cluster_print_chain output to adjust
% cluster control parameters before running or to run the next lines on a
% different computer (the save/load functions are for this purpose).

return
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.desired_time_per_job',5*60);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.cpu_time_mult',2);
%ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.mem_mult',2);

[ctrl_chain,chain_fn] = cluster_load_chain([],chain_id);
ctrl_chain = cluster_run(ctrl_chain);






return


% User Settings
% =========================================================================

param = [];
counter_correction_en = false;
online_mode = false;

% Enable Just One Radar Setup
radar_setup = 'snow8';

%% Accum 1
if strcmpi(radar_setup,'ACCUM')
  param.radar_name = 'accum';
  adcs = 1;
  param.file_version = 5;
  raw_file_suffix = '.dat';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2010_Greenland_P3';
  base_dir = '/cresis/snfs1/data/Accum_Data/2010_Greenland_P3/';
  param.adc_folder_name = '20100513B';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20100513'; % Only used for stdout print of the vectors worksheet
  % file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'data'; % (pre 2011, and beginning of 2011)
end

%% Ka-band 3
if strcmpi(radar_setup,'KABAND3')
  param.radar_name = 'kaband3';
  param.clk = 125e6;
  adcs = 1;
  param.file_version = 5;
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = '';
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2015_Greenland_C130';
  base_dir = '/cresis/snfs1/data/kaband/';
  param.adc_folder_name = '20150328';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20150328'; % Only used for stdout print of the vectors worksheet
end

%% Kuband 1
if strcmpi(radar_setup,'KUBAND1')
  param.radar_name = 'kuband';
  adcs = 1;
  param.file_version = 1; % 2 for 2012
  raw_file_suffix = '.dat'; % kuband3, snow3
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2011_Greenland_P3';
  base_dir = '/cresis/snfs1/data/Ku-Band/2011_Greenland_P3/';
  param.adc_folder_name = '20110316';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20110316'; % Only used for stdout print of the vectors worksheet  
  % file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'data'; % (pre 2011, and beginning of 2011)
end

%% Kuband 2
if strcmpi(radar_setup,'KUBAND2')
  param.radar_name = 'kuband2';
  param.clk = 125e6;
  adcs = 1;
  param.file_version = 2; % 2 for 2012
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = 'kuband';
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2012_Greenland_P3';
  base_dir = '/cresis/snfs1/data/Ku-Band/';
  param.adc_folder_name = '20120316';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20120316'; % Only used for stdout print of the vectors worksheet
end

%% Ku-band 3
if strcmpi(radar_setup,'KUBAND3')
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
  base_dir = '/process3/20161115/fmcw/kuband/';
  param.adc_folder_name = '';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20161115'; % Only used for stdout print of the vectors worksheet
end

%% RDS: ACORDS
if strcmpi(radar_setup,'ACORDS')
  param.radar_name = 'acords';
  adcs = 1;
  param.file_version = 406;
  raw_file_suffix = '';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2004_Antarctica_P3chile';
  % base_dir = '/cresis/snfs1/data/ACORDS/airborne2005/';
  base_dir = '/cresis/snfs1/data/ACORDS/Chile_2004/';
  param.adc_folder_name = 'nov21_04';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20041121'; % Only used for stdout print of the vectors worksheet
  file_prefix_override = 'nov21_04'; % most of the time
  file_regexp = '\.[0-9]*$';
end

%% RDS: MCoRDS 3
if strcmpi(radar_setup,'MCORDS3')
  param.radar_name = 'mcords3';
  param.clk = 1e9/9;
  adcs = [1 5 9 13];
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = true;
  presum_bug_fixed = false;
  union_time_epri_gaps = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2018_Greenland_P3';
  base_dir = '/process-archive/';
%   base_dir = '/net/field1/landing/mcords/';
  param.adc_folder_name = '20180414/mcords/board%b';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20180414'; % Only used for stdout print of the vectors worksheet
end

%% RDS: MCoRDS 4
if strcmpi(radar_setup,'MCORDS4')
  base_dir = '/N/dc2/projects/cresis/2013_Antarctica_DC3/20131216/mcords4/';
  param.radar_name = 'mcords4';
  adc_folder_names = {'chan1'};

  param.file_version = 404;

  file_midfix = ''; % Can often be left empty
  day_string = '20131216'; % Only used during printing of the segments
  param.season_name = '2013_Antarctia_Basler';
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  % file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'mcords4'; 
end

%% RDS: MCORDS5
if strcmpi(radar_setup,'MCORDS5')
  param.radar_name = 'mcords5';
  param.clk = 1.6e9/8;
  adcs = 1:8;
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = true;
  presum_bug_fixed = false; % Seasons from 2015 Greenland Polar6 onward should be set to true, except for 2017 Antarctica Basler which uses the cresis DDS with this bug
  union_time_epri_gaps = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2018_Greenland_P3';
  base_dir = '/process-archive/';
  param.adc_folder_name = '20180414/accum/chan%d';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20180414'; % Only used for stdout print of the vectors worksheet
end

%% Snow 1
if strcmpi(radar_setup,'SNOW1')
  param.radar_name = 'snow';
  adcs = 1;
  param.file_version = 1; % 2 for 2012
  raw_file_suffix = '.dat'; % kuband3, snow3
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2011_Greenland_P3';
  base_dir = '/cresis/snfs1/data/SnowRadar/2011_Greenland_P3/';
  param.adc_folder_name = '20110316';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20110316'; % Only used for stdout print of the vectors worksheet  end
  % file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'data'; % (pre 2011, and beginning of 2011)
end

%% Snow 2
if strcmpi(radar_setup,'SNOW2')
  param.radar_name = 'snow2';
  param.clk = 125e6;
  adcs = 1;
  param.file_version = 2; % 2 for 2012, 4 for 2018_Alaska_SO (SIERRA) 
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = 'snow';
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2012_Greenland_P3';
  if strcmpi(param.season_name,'2018_Alaska_SO')
      param.nohack = 1;
  end
  base_dir = '/cresis/snfs1/data/SnowRadar/';
  param.adc_folder_name = '20120316';
  file_midfix = '20120316'; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20120316'; % Only used for stdout print of the vectors worksheet
end

%% Snow 3 (OIB)
if strcmpi(radar_setup,'SNOW3')
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
  day_string = '20161017'; % Only used for stdout print of the vectors worksheet
end

%% Snow 3 (NRL)
if strcmpi(radar_setup,'SNOW3_NRL')
  param.radar_name = 'snow3';
  param.clk = 125e6;
  adcs = 1:10;
  param.file_version = 5; % 3 for 2013 Gr, 5 for after that
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = '';
  counter_correction_en = true;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2015_Alaska_TOnrl';
  base_dir = '/cresis/snfs1/data/SnowRadar/';
  param.adc_folder_name = '20150328/chan%02d';
  file_midfix = '20150328'; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20150328'; % Only used for stdout print of the vectors worksheet
end

%% SNOW5
if strcmpi(radar_setup,'SNOW5')
  param.radar_name = 'snow5';
  param.clk = 125e6;
  adcs = 1:2;
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = false;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2016_Greenland_Polar6';
  base_dir = '/mnt/HDD6/';
  param.adc_folder_name = '1604180601/UWBM/chan%d';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20160418'; % Only used for stdout print of the vectors worksheet
  expected_rec_sizes = [60480      120864      181296];
end

%% SNOW8
if strcmpi(radar_setup,'SNOW8')
  param.radar_name = 'snow8';
  param.clk = 125e6;
  adcs = 1;
  param.file_version = 8;
  raw_file_suffix = '.bin';
  reuse_tmp_files = false; % Set to false if you want to overwrite current results
  file_prefix_override = 'snow8'; % most of the time
  counter_correction_en = false;
  
  % Parameters below this point OFTEN NEEDS TO BE CHANGED
  param.season_name = '2018_Greenland_P3';
  base_dir = '/N/dcwan/projects/cresis/2018_Greenland_P3';
  param.adc_folder_name = '20180414/fmcw/snow/';
  file_midfix = ''; % Data files must contain this string in the middle of their name (usually should be empty)
  day_string = '20180414'; % Only used for stdout print of the vectors worksheet
end

%% User Settings that should not generally be changed
% You may have to set to false to read some of the results from this function when it was first written (should always be true)
tmp_fn_uses_adc_folder_name = true;


MIN_SEG_SIZE = 2;
MAX_TIME_GAP = 1000/75;
MIN_PRF = 100;

%% Automated Section
% =========================================================================

create_segment_raw_file_list_v2;

return;
