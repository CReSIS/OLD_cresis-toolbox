% script run_basic_rx_chan_equalization
%
% This script is for helping with setting the receiver coefficients
% from raw data. It requires loading one waveform and N receive channels
% and then analyzing these data.
%
% 1. Collect data with so that receiver gain settings are the same as
%    the ordinary data collection, but the air/ice surface response
%    is unsaturated (ideally from high altitude over the ocean)
%    If the TR switch is slow, you will need to collect the data from
%    a high enough altitude where the fast-time gain changes of the TR
%    switch are stable (this is 10 us after the switch change for the
%    original MCoRDS system).
% 2. Pulse duration should not matter, but probably affects results
%    slightly
% 3. Time gate should be large enough to include noise-only data which
%    will be used to determine the SNR of the surface return (SNR
%    threshold is used to exclude low SNR data points from the measurement)
%
% Author: John Paden, Logan Smith

error('Not completed yet');

clear param;
clear pc_param;
physical_constants;
basic_rx_chan_equalization_tstart = tic;

% =======================================================================
% User Settings: Common Settings
% =======================================================================

% use_ginput: use ginput to select signal region and noise region
use_ginput = false;
if ~use_ginput
  % Do not manually select a region with ginput, use these instead
  
  % .rlines = Range lines to process from dataset
  %   These are range lines post presumming
  param.rlines = [1 inf];
  
  % .rbins = Range bins to search for surface in
  %   These range bins are post presumming
  param.rbins = [350:550];
  
  % .noise_rbins,rlines: Use for noise power calculation
  param.noise_rlines = [1 inf];
  param.noise_rbins = [200 300];
end

% Specify pulse compression properties for this waveform
pc_param.f0 = 180e6;
pc_param.f1 = 210e6;
pc_param.f0 = 189.15e6;
pc_param.f1 = 198.65e6;
pc_param.Tpd = 10e-6;
pc_param.tukey = 0.2;
pc_param.window_func = @hanning;

radar_name = 'mcords';
if strcmpi(radar_name,'mcrds')
  % Sampling frequency of radar (required to read data in)
  fs = 120e6;
  
  % .ref_wf_adc_idx = Reference receive channel, index into param.adcs,
  %   (surface location determined from this channel and all phase
  %   measurements made relative to it)
  param.ref_wf_adc_idx = 3;
  
  % .img = which waveform/adc pairs to load
  param.img = [2 1; 2 2; 2 3; 2 4; 2 5; 2 6];
  
  % Correction coefficients
  %   param.td = zeros(size(param.img,1),1);
  %   param.amp = zeros(size(param.img,1),1);
  %   param.phase = zeros(size(param.img,1),1);
  param.td = 1e-9*[0 0 0 0 0 0].';
  param.amp = [0 0 0 0 0 0].';
  param.phase = [0 0 0 0 0 0].';
  
  
  % Base path of data (does not include seg directory)
  base_path = '/cresis/data1/MCRDS/2008_Greenland/20080702A/';
  fn_prefix = 'data.';
  
  % Which segment directory to load from, leave empty for no segment
  % directory
  seg = '';
  
  file_adc_folder_name = 'dataGISMO/';
  
  % File index in filename
  file_nums = [470:474 390:395]; % 90,166, rolls: 330
  
elseif strcmpi(radar_name,'mcords')
  % Sampling frequency of radar (required to read data in)
  fs = 1e9/9;
  
  % .ref_wf_adc_idx = Reference receive channel, index into param.adcs, 
  %   (surface location determined from this channel and all phase
  %   measurements made relative to it)
  param.ref_wf_adc_idx = 3;
  
  % .img = which waveform/adc pairs to load
%   param.img = cat(2,1*ones(6,1),[1 2 3 4 5 6].');
  param.img = [2 1; 2 2; 2 3; 2 4; 2 5];
 
    % Correction coefficients
  param.td = zeros(size(param.img,1),1);
%   param.td = 1e-9*[6.3	0.0	1.1	11.6	30.9	28.5].';
  param.amp = zeros(size(param.img,1),1);
%   param.amp = [-7.4	0.0	0.4	-0.7	-0.9	-1.9].';
  param.phase = zeros(size(param.img,1),1);
%   param.phase = [118.9	-0.0	101.5	87.6	-18.4	152.7].';
  
  % Base path of data (does not include seg directory)
  base_path =     '/N/dcwan/projects/cresis/2011_Chile_DC8/20111014/';

  % Which segment directory to load from, leave empty for no segment
  % directory
  seg = 'seg01';

  % Channel directory naming convention (usually "chan%d" or "ch%d")
  chan = 'ch%d';
  
  % Radar ID (there are 4 mcords systems, IDs 1-4, ID is in filename as rX)
  radar_id = 1;
  
  fn_start = 'mcords';
  
  % File index in filename
  file_nums = 0;
  
  trim = [1 1];
  
elseif strcmpi(radar_name,'mcords2')
  % Sampling frequency of radar (required to read data in)
  fs = 1e9/9;
  
  % .ref_wf_adc_idx = Reference receive channel, index into param.adcs,
  %   (surface location determined from this channel and all phase
  %   measurements made relative to it)
  param.ref_wf_adc_idx = 1;
  
  % .img = which waveform/adc pairs to load
  param.img = cat(2,1*ones(12,1),[2 3 4 5 6 7 9 10 11 12 13 14].');
%   param.img = cat(2,1*ones(6,1),[2 3 4 5 6 7 ].');
  
  % Correction coefficients
  param.td = zeros(14,1);
  param.td = 1e-9*[0.0	0.0	-2.2	-1.2	0.5	0.5	4.7	0.0	-82.4	-82.6	-79.5	-81.3	-80.0	-78.1].';
  param.amp = zeros(14,1);
  param.amp = [0.0	0.0	-4.3	-2.9	-4.4	-4.0	-6.5	0.0	-6.7	-1.3	-0.5	-2.1	0.5	-0.4].';
  param.phase = zeros(14,1);
  param.phase = [0.0	-0.0	-149.0	-60.1	23.0	-127.2	-4.4	0.0	89.7	3.2	-128.3	80.0	-168.2	-81.0].'
  
  % base_path = Base path of data (does not include seg directory)
  % set = Which segment directory to load from, leave empty for no segment directory
  % utc_time_correction
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111202/'; seg = 'seg_00'; utc_time_correction = 86400;% 39 files
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111206/'; seg = 'seg_00'; utc_time_correction = 0;% 29 files
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111206/'; seg = 'seg_08'; utc_time_correction = 86400;% 104 files (ROLLS)
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111207/'; seg = 'seg_00'; utc_time_correction = 0; param.rbins = [750:1250]; % 70 files
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111210/'; seg = 'seg_00'; utc_time_correction = 0;% 64 files
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111213/'; seg = 'seg_00'; utc_time_correction = 0;% 50 files
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111214/'; seg = 'seg_00'; utc_time_correction = 0; param.rbins = [750:1250]; % 80 files
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111215/'; seg = 'seg_00'; utc_time_correction = 0;% 50 files (BROKEN)
%     base_path = '/cresis/data4/MCoRDS/2011_Antarctica_TO/20111217/'; seg = 'seg_00'; utc_time_correction = 0; param.rbins = [750:1250]; % 30 files
%     base_path = '/cresis/data3/MCoRDS/2011_Antarctica_TO/20111219/'; seg = 'seg_00'; utc_time_correction = 0;% 115 files
%     base_path = '/cresis/data3/MCoRDS/2011_Antarctica_TO/20111220/'; seg = 'seg_00'; utc_time_correction = 0;% 105 files
      base_path = '/mnt/source/20120504/mcords/';
      
  
  
  % Optionally restrict search to a particular acquisition number
  % (part of the data files' filenames)
  acquisition_num = [0];
  
  % File index in filename
  file_nums = [1];
 
end

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = '/N/dc/projects/cresis/csarp_support/gps/2011_Antarctica_DC8/gps_20111014.mat';
% utc_time_correction = 0;

% =======================================================================
% User Settings: Usually not changed
% =======================================================================

% param.type sent to motion_comp.m
param.mocomp_type = 4;

% Transmit weights sent to lever_arm_fh
% param.tx_weights = [1 1 1 1 1 1 1];
% param.tx_weights = [1 1 1 1 1 1];
param.tx_weights = [1 1 1 1 1];

% Map of adcs to receivers for each waveform, used with lever_arm_fh
param.rx_paths = {};
param.rx_paths{1} = [1 2 3 4 5];
param.rx_paths{2} = [1 2 3 4 5];


% lever_arm_fh = lever arm
param.lever_arm_fh = @lever_arm_mcords_2011_antarctica_DC8_ATM;

% .plot_en = flag to enable plots
param.plot_en = true;

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR are included in the estimate)
param.snr_threshold = 12;

% presums = Number of presums (coherent averaging) to do
presums = 10;

% .averaging_fh = method to average complex vectors when computing
%   recommended channel compensation (mean is ideal, but median
%   may be necessary to remove outliers)
param.averaging_fh = @mean;

% Specify if cross correlation should be used (required for finding td)
% or if peak finding should be used when comparing channels
param.cross_correlation_flag = 0;

% Combine channels for surface tracker?
param.combine_channels = false;

% .ref_bins: reference bins around surface peak bin to use in correlation
param.ref_bins = [-7 7];
param.ref_bins = [-1 1];
% .search_bins: neighboring bins to search for a peak in other channels
param.search_bins = [-12 12];
param.search_bins = [-1 1];
% .Mt: oversampling factor for correlation output
param.Mt = 100;

basic_rx_chan_equalization;

return;



