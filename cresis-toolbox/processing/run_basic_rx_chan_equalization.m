
% script run_basic_rx_chan_equalization
%
% Runs basic_rx_chan_equalization
%
% Author: John Paden, Logan Smith

% Enable Just One Radar Setup
radar_setup = 'MCORDS3_DC8';

%% RDS: MCORDS5
if strcmpi(radar_setup,'MCORDS5')
  [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file',
  %   'specific', 'default', or empty to be asked
  param.file_search_mode = '';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'D:\awi\','/mnt/AWI_SSD0/1604261101/UWB/','/mnt/AWI_SSD0/1604261202/UWB/'};
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  wf = 1; adcs = 1:24; ref = 12; % Waveform 1, All ADCs
  %wf = 1; adcs = 1:8; ref = 4; % Waveform 1, Left subarray
  %wf = 1; adcs = 9:16; ref = 4; % Waveform 1, Center subarray
  %wf = 1; adcs = 17:24; ref = 4; % Waveform 1, Right subarray
  param.img = cat(2,wf*ones(length(adcs),1),adcs.'); param.ref_wf_adc = ref;
  
  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 inf];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 20;
  
  % .delay: the method used to calculate delay between the channels
  param.delay = struct('method','xcorr_complex','bin_rng',-20:20,'Mt',10);
end

%% RDS: MCORDS3_P3
if strcmpi(radar_setup,'MCORDS3_P3')
  [param,defaults] = default_radar_params_2016_Greenland_P3_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file', 'specific'
  param.file_search_mode = 'specific';
  param.multiple_files = true;

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'/cresis/snfs1/data/MCoRDS/20140426/'};
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  param.img = cat(2,3*ones(15,1),[1:15].'); param.ref_wf_adc = 3;
  
  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 inf];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 1;
  
  % .delay: the method used to calculate delay between the channels
  param.delay = struct('method','xcorr_complex','bin_rng',-20:20,'Mt',10)
end

%% RDS: MCORDS3_DC8
if strcmpi(radar_setup,'MCORDS3_DC8')
  [param,defaults] = default_radar_params_2016_Antarctica_DC8_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file', 'specific'
  param.file_search_mode = 'segment';
  param.multiple_files = true;

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'/process/20161014/mcords/'};
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  param.img = cat(2,5*ones(6,1),[1:6].'); param.ref_wf_adc = 5;
  
  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 inf];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 1;
  
  % .delay: the method used to calculate delay between the channels
  param.delay = struct('method','xcorr_complex','bin_rng',-20:20,'Mt',10)
end

%% Automated Section

basic_rx_chan_equalization(param,defaults);
