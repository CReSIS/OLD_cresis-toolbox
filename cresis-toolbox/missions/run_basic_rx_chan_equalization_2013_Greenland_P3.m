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
  param.rbins = [1000:1150];
  
  % .noise_rbins,rlines: Use for noise power calculation
  param.noise_rlines = [1 inf];
  param.noise_rbins = [1600 1700];
end

% Specify pulse compression properties for this waveform
pc_param.f0 = 180e6;
pc_param.f1 = 210e6;
pc_param.Tpd = 10e-6;
pc_param.tukey = 0.2;
pc_param.window_func = @hanning;

radar_name = 'mcords3';
if strcmpi(radar_name,'mcords3')
  % Sampling frequency of radar (required to read data in)
  fs = 1e9/9;
  
  % .ref_wf_adc_idx = Reference receive channel, index into param.adcs,
  %   (surface location determined from this channel and all phase
  %   measurements made relative to it)
  param.ref_wf_adc_idx = 3;
  
  % .img = which waveform/adc pairs to load
  param.img = cat(2,2*ones(7,1),[1 2 3 4 5 6 7].');
  
  % Correction coefficients
  param.td = zeros(7,1);
%   param.td = 1e-9*[0.0	0.0	-2.2	-1.2	0.5	0.5	4.7	0.0	-82.4	-82.6	-79.5	-81.3	-80.0	-78.1].';
%   param.amp = zeros(7,1);
  param.amp = [1.3	-0.3	0.0	-2.0	0.0	0.7	1.2].';
  % 0.9	-1.1	0.0	-2.5	0.8	1.1	0.5
  % 0.9	-1.1	0.0	-2.1	0.1	0.4	0.8
  % [1.63 0.43 0 -2.73 -0.13 0.3 0.9]
%   param.phase = zeros(7,1);
  param.phase = [-3.1	-44.1	-0.0	-124.4	-6.6	-14.2	7.3].'
  % -4.5	-60.3	0.0	-137.5	-29.7	-40.5	-26.8 <-- low rolls
  % 10.8	-39.5	0.0	-132.2	-28.1	-42.0	-29.7
  % 36.7	-35.9	0.0	-152.9	-56.2	-80.4	-91.3 <-- low rolls
  % [5 -41 0 -134 -19 -28 -15]
  
  % base_path = Base path of data (does not include seg directory)
  % set = Which segment directory to load from, leave empty for no segment directory
  % utc_time_correction
%     base_path = '/cresis/data3/MCoRDS/2011_Antarctica_TO/20111220/'; seg = 'seg_00'; utc_time_correction = 0;% 105 files
      base_path = '/N/dc2/projects/cresis/2013_Greenland_P3/20130426/mcords/'; seg = ''; utc_time_correction = 0;% 105 files
      
  param.season_name = '2013_Greenland_P3';
param.radar_name = 'mcords3';
  
  % Optionally restrict search to a particular acquisition number
  % (part of the data files' filenames)
  acquisition_num = [0];
  
  % File index in filename
  file_nums = [1444:1448];
%   file_nums = [1444:1448]-15; %<-- low roll angles
 
end

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = '/N/dc2/projects/cresis/csarp_support/gps/2013_Greenland_P3/gps_20130426.mat';
% param.gps_fn = '';
% utc_time_correction = 0;

% =======================================================================
% User Settings: Usually not changed
% =======================================================================

% param.type sent to motion_comp.m
param.mocomp_type = 4;

% Transmit weights sent to lever_arm_fh
param.tx_weights = [1 1 1 1 1 1 1];

% Map of adcs to receivers for each waveform, used with lever_arm_fh
param.rx_paths = {};
param.rx_paths{1} = [1 2 3 4 5 6 7];
param.rx_paths{2} = [1 2 3 4 5 6 7];


% lever_arm_fh = lever arm
param.lever_arm_fh = @lever_arm;

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



