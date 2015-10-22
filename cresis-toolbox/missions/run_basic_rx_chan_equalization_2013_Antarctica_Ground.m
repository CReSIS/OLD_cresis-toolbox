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
  param.rbins = [3000:3500];
  
  % .noise_rbins,rlines: Use for noise power calculation
  param.noise_rlines = [1 inf];
  param.noise_rbins = [1000 1500];
end

% Specify pulse compression properties for this waveform
pc_param.f0 = 200e6;
pc_param.f1 = 450e6;
pc_param.Tpd = 3e-6;
pc_param.tukey = 0.2;
pc_param.window_func = @hanning;

param.season_name = '2013_Antarctica_Basler';
param.radar_name = 'mcords4';
radar_name = 'mcords4';
if strcmpi(radar_name,'mcords4')
  % Sampling frequency of radar (required to read data in)
  fs = 1e9/2;
  
  % .ref_wf_adc_idx = Reference receive channel, index into param.adcs,
  %   (surface location determined from this channel and all phase
  %   measurements made relative to it)
  param.ref_wf_adc_idx = 4;
  
  % .img = which waveform/adc pairs to load
  param.img = cat(2,-j*9*ones(8,1),[1 2 3 4 5 6 7 8].');
  
  % Correction coefficients
  param.td = zeros(8,1);
  %   param.td = 1e-9*[0 0 0 0 0 0 0 0].';
  % param.amp = zeros(8,1);
  param.amp = [-0.9	-0.1	-0.1	0.0	-0.0	0.0	-0.3	-1.5].';
  % param.phase = zeros(8,1);
  param.phase = [-7.9	5.6	13.9	0.0	-0.0	10.2	-12.1	-1.9].'

  % base_path = Base path of data (does not include seg directory)
  % set = Which segment directory to load from, leave empty for no segment directory
  % utc_time_correction
  base_path = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/test_flights_20130921/';
  seg = '';
  
  % Optionally restrict search to a particular acquisition number/time
  % (part of the data files' filenames)
  acquisition_num = '20130921_20*03';
  
  % File index in filename
  file_nums = [0 1 2];
  
end

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = '/cresis/projects/dev/cr1/gps/2013_Antarctica_Basler/gps_20130921.mat';
utc_time_correction = 2;

% =======================================================================
% User Settings: Usually not changed
% =======================================================================

% param.type sent to motion_comp.m
param.mocomp_type = 4;

% Transmit weights sent to lever_arm_fh
param.tx_weights = [0 0 0 1 0 0 0 0];

% Map of adcs to receivers for each waveform, used with lever_arm_fh
param.rx_paths = {};
param.rx_paths{9} = [1 2 3 4 5 6 7 8];


% lever_arm_fh = lever arm
param.lever_arm_fh = @lever_arm;

% .plot_en = flag to enable plots
param.plot_en = true;

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR are included in the estimate)
param.snr_threshold = 12;

% presums = Number of presums (coherent averaging) to do
presums = 4;

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

param.multilook = [1];

basic_rx_chan_equalization;

return;



