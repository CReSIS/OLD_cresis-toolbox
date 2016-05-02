% script run_basic_rx_chan_equalization
%
% Runs basic_rx_chan_equalization
%
% Author: John Paden, Logan Smith

% Initialization
physical_constants;
param = [];

% Enable Just One Radar Setup
radar_setup = 'MCORDS5';

%% RDS: MCORDS5
if strcmpi(radar_setup,'MCORDS5')
  [param,defaults] = default_radar_params_2016_Greenland_Polar6_mcords;
  
  % .file_search_mode: Specify how to search for a file: 'last_file',
  %   'specific', 'default', or empty to be asked
  param.file_search_mode = '';

  % .base_dir_search: cell vector of paths to search for data files
  param.base_dir_search = {'D:\awi\chan1','/mnt/AWI_SSD0/1604261101/UWB/','/mnt/AWI_SSD0/1604261202/UWB'};
  
  % .img: wf-adc pair list which specifies which waveform-adc pairs to
  %   analyze
  param.img = cat(2,1*ones(24,1),[1:24].'); param.ref_wf_adc = 12;
%   param.img = cat(2,2*ones(24,1),[1:24].'); param.ref_wf_adc = 12;
%   param.img = cat(2,3*ones(24,1),[1:24].'); param.ref_wf_adc = 12;
%   param.img = cat(2,1*ones(8,1),[9:16].'); param.ref_wf_adc = 4; % Center subarray
%   param.img = cat(2,1*ones(8,1),[1:8].'); param.ref_wf_adc = 4; % Left subarray
%   param.img = cat(2,1*ones(8,1),[17:24].'); param.ref_wf_adc = 4; % Right subarray
  
  % .recs: two element vector specifying which records/range-lines to load
  %   [start_record num_records]
  param.recs = [0 inf];
  
  % .presums: Number of additional software presums (coherent averaging) to do
  param.presums = 20;
  
  % .delay: the method used to calculate delay between the channels
  param.delay = struct('method','xcorr_complex','bin_rng',-20:20,'Mt',10);
end

%% RDS: MCORDS3
if strcmpi(radar_setup,'MCORDS3')
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

%% Automated Section

basic_rx_chan_equalization(param,defaults);

return;



