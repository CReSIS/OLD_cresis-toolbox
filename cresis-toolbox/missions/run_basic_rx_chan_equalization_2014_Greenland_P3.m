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
  param.rbins = [700:900];
  
  % .noise_rbins,rlines: Use for noise power calculation
  param.noise_rlines = [1 inf];
  param.noise_rbins = [200 300];
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
  param.ref_wf_adc_idx = 4;
  
  % .img = which waveform/adc pairs to load
  param.img = cat(2,param.ref_wf_adc_idx*ones(15,1),[2 3 4 5 6 7 8 9 10 11 12 13 14 15 16].');
  
  % Correction coefficients
  param.td = zeros(15,1);
%   param.td = 1e-9*[-5.03	-5.29	-7.56	0.00	-5.40	-7.50	-5.98	-27.44	-31.98	-32.67	-37.96	-37.71	-33.76	-34.00	-28.37].';
  param.amp = zeros(15,1);
%   param.amp = [4.6	3.5	2.4	0.0	2.3	3.4	4.3	4.8	1.7	3.1	4.1	4.9	3.3	1.7	4.8].';
  param.phase = zeros(15,1);
%   param.phase = [157.5	83.3	-15.2	0.0	103.8	-27.1	91.5	144.6	49.8	-108.5	153.1	-145.5	-166.0	-14.0	152.9].'
  
  % base_path = Base path of data (does not include seg directory)
  % set = Which segment directory to load from, leave empty for no segment directory
  % utc_time_correction
      base_path = 'D:\20140310_26_calibration\'; seg = ''; utc_time_correction = 1;
  param.season_name = '2014_Greenland_P3';
  param.radar_name = 'mcords3';
  
  % Optionally restrict search to a particular acquisition number
  % (part of the data files' filenames)
  acquisition_num = [7];
  
  % File index in filename
  file_nums = [2:8];
%   file_nums = [1444:1448]-15; %<-- low roll angles
 
end

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = 'C:\csarp_support\gps\2014_Greenland_P3\gps_20140310.mat';
% param.gps_fn = '';
% utc_time_correction = 0;

% =======================================================================
% User Settings: Usually not changed
% =======================================================================

% param.type sent to motion_comp.m
param.mocomp_type = 4;

% Transmit weights sent to lever_arm_fh
param.tx_weights = [1 0 0 0 0 0 0];

% Map of adcs to receivers for each waveform, used with lever_arm_fh
param.rx_paths = {};
param.rx_paths{4} = [1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];


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



