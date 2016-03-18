% script run_create_segment_raw_file_list_v2.m
%
% Runs script create_segment_raw_file_list_v2.m
%
% Instructions:
% 1. Find your radar in the "user settings" section below
% 2. Enable it ("if 1") and disable all others ("if 0")
% 3. Change settings in the radar section as required (usually
%    just path changes to match the current date).
% 4. Run the script
%
% Author: John Paden

% =========================================================================
%% User Settings
% =========================================================================

param = [];
counter_correction_en = false;

% Enable Just One Radar Setup

if 0
  base_dir = '/N/dc2/projects/cresis/2013_Antarctica_DC3/20131216/mcords4/';
  param.radar_name = 'mcords4';
  adc_folder_names = {'chan1'};

  param.file_version = 404;

  file_midfix = ''; % Can often be left empty
  day_string = '20131216'; % Only used during printing of the segments
  param.season_name = '2013_Antarctia_Basler';
  raw_file_suffix = '.bin';
  reuse_tmp_files = false; % Set to false if you want to overwrite current results
  % file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'mcords4'; 
end


if 0
  param.radar_name = 'accum';
  base_dir = '/cresis/snfs1/data/Accum_Data/2010_Greenland_P3/';
  adc_folder_names = {'20100513B'}; %27 28 29 31
  param.file_version = 5;

  file_midfix = ''; % Can often be left empty
  day_string = '20110316'; % Only used during printing of the segments
  param.season_name = '2011_Greenland_P3';
  raw_file_suffix = '.dat'; % accum
  reuse_tmp_files = false; % Set to false if you want to overwrite current results
  % file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'data'; % (pre 2011, and beginning of 2011)
end

if 1
  param.radar_name = 'snow';
  base_dir = '/cresis/snfs1/data/SnowRadar/2010_Greenland_P3/';
  adc_folder_names = {'20100508A'};
  param.file_version = 1;

  file_midfix = ''; % Can often be left empty
  day_string = '20100508'; % Only used during printing of the segments
  param.season_name = '2010_Greenland_P3';
  raw_file_suffix = '.dat'; % kuband3, snow3
  reuse_tmp_files = false; % Set to false if you want to overwrite current results
  %file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'data'; % (pre 2011, and beginning of 2011)
end

if 0
  param.radar_name = 'kuband';
  base_dir = '/cresis/snfs1/data/Ku-Band/2011_Greenland_P3/';
  adc_folder_names = {'20110316'};
  param.file_version = 1;

  file_midfix = ''; % Can often be left empty
  day_string = '20110316'; % Only used during printing of the segments
  param.season_name = '2011_Greenland_P3';
  raw_file_suffix = '.dat'; % kuband3, snow3
  reuse_tmp_files = false; % Set to false if you want to overwrite current results
  % file_prefix_override = ''; % most of the time (most of 2011)
  file_prefix_override = 'data'; % (pre 2011, and beginning of 2011)
end

if 0
  data_day = '21';
  param.radar_name = 'acords';
%   base_dir = '/cresis/snfs1/data/ACORDS/airborne2005/';
  base_dir = '/cresis/snfs1/data/ACORDS/Chile_2004/';
  adc_folder_names = {sprintf('nov%s_04',data_day)}; %27 28 29 31
  param.file_version = 406;
  
  
  file_midfix = ''; % Can often be left empty
  day_string = sprintf('200411%s',data_day(1:2)); % Only used during printing of the segmen1
  param.season_name = '2004_Antarctica_P3chile';
  % file_prefix_override = ''; % most of the time
  file_prefix_override = sprintf('nov%s_04',data_day); % for a few older datasets where file prefix was not radar name
  % day_string = '20110328'; % Only used during printing of the segments
  % param.season_name = '2011_Greenland_P3';
  % file_prefix_override = 'data'; % for a few older datasets where file prefix was not radar name
  % raw_file_suffix = '.bin'; % kuband3, snow3
  file_regexp = '\.[0-9]*$';
  raw_file_suffix = ''; % accum
  reuse_tmp_files = false;
end

if 0
  param.radar_name = 'snow2';
  base_dir = '/cresis/snfs1/data/SnowRadar/2012_Greenland_P3/';
  adc_folder_names = {'20120316'};
  param.file_version = 2; % 2 for 2012

  file_midfix = ''; % Can often be left empty
  day_string = '20120316'; % Only used during printing of the segments
  param.season_name = '2012_Greenland_P3';
  raw_file_suffix = '.bin'; % kuband3, snow3
  reuse_tmp_files = false; % Set to false if you want to overwrite current results
  file_prefix_override = 'snow'; % most of the time
end

if 0
  param.radar_name = 'kaband3';
  base_dir = '/process4_SSD/20150324/fmcw/kaband/';
  adc_folder_names = {''};
  param.file_version = 5; % 3 for 2013 Gr, 5 for after that

  file_midfix = ''; % Can often be left empty
  day_string = '20150324'; % Only used during printing of the segments
  param.season_name = '2015_Greenland_LC130';
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = true;
end

if 0
  param.radar_name = 'kuband3';
  base_dir = '/process4_SSD/20150324/fmcw/kuband/';
  adc_folder_names = {''};
  param.file_version = 5; % 3 for 2013 Gr, 5 for after that

  file_midfix = ''; % Can often be left empty
  day_string = '20150324'; % Only used during printing of the segments
  param.season_name = '2015_Greenland_LC130';
  raw_file_suffix = '.bin'; % kuband3, snow3
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = true;
end

if 0
  param.radar_name = 'snow3';
  param.clk = 125e6;
  base_dir = '/cresis/snfs1/data/SnowRadar/20150328/';
  if 0
    % Single channel
    adc_folder_names = {''};
  else
    % Multichannel
    adc_folder_names_chans = [1:10];
    adc_folder_names = {};
    for adc_folder_names_chan_idx = 1:length(adc_folder_names_chans)
      adc_folder_names{adc_folder_names_chan_idx} = sprintf('chan%02d',adc_folder_names_chans(adc_folder_names_chan_idx));
    end
  end
  param.file_version = 5; % 3 for 2013 Gr, 5 for after that

  file_midfix = '20150328'; % Can often be left empty
  day_string = '20150328'; % Only used during printing of the segments
  param.season_name = '2015_Alaska_TOnrl';
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = true;
end

if 0
  param.radar_name = 'mcords5';
  param.clk = 1.6e9/8;
  base_dir = '/cresis/snfs1/data/MCoRDS/2015_Greenland_Polar6/';
  adc_folder_names_chans = [1:24];
  adc_folder_names = {};
  for adc_folder_names_chan_idx = 1:length(adc_folder_names_chans)
    adc_folder_names{adc_folder_names_chan_idx} = sprintf('20150916/chan%d',adc_folder_names_chans(adc_folder_names_chan_idx));
  end

  file_midfix = '20150916'; % Can often be left empty
  day_string = '20150916'; % Only used during printing of the segments
  param.season_name = '2015_Greenland_Polar6';
  raw_file_suffix = '.bin';
  reuse_tmp_files = true; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = true;
  presum_bug_fixed = true; % Seasons from 2015 Greenland Polar6 onward should be set to true
end

if 0
  param.radar_name = 'snow5';
  param.clk = 125e6;
  base_dir = '/cresis/snfs1/data/SnowRadar/2015_Greenland_Polar6/20150917';
  adc_folder_names_chans = [1:2];
  adc_folder_names = {};
  for adc_folder_names_chan_idx = 1:length(adc_folder_names_chans)
    adc_folder_names{adc_folder_names_chan_idx} = sprintf('chan%d',adc_folder_names_chans(adc_folder_names_chan_idx));
  end

  file_midfix = '20150917'; % Can often be left empty
  day_string = '20150917'; % Only used during printing of the segments
  param.season_name = '2015_Greenland_Polar6';
  raw_file_suffix = '.bin';
  reuse_tmp_files = false; % Set to false if you want to overwrite current results
  file_prefix_override = ''; % most of the time
  counter_correction_en = false;
end

%% User Settings that should not generally be changed
% You may have to set to false to read some of the results from this function when it was first written (should always be true)
tmp_fn_uses_adc_folder_name = true;

MIN_SEG_SIZE = 2;
MAX_TIME_GAP = 1000/75;
MIN_PRF = 100;

% =========================================================================
%% Automated Section
% =========================================================================

create_segment_raw_file_list_v2;

return;
