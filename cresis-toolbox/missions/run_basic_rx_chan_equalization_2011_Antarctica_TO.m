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
  % SET BELOW
  
  % .noise_rbins,rlines: Use for noise power calculation
  param.noise_rlines = [1 inf];
  param.noise_rbins = [250 300];
end

% Specify pulse compression properties for this waveform
pc_param.f0 = 180e6;
pc_param.f1 = 210e6;
pc_param.tukey = 0.0;
pc_param.window_func = @hanning;

radar_name = 'mcords2';
if strcmpi(radar_name,'mcords2')
  % Sampling frequency of radar (required to read data in)
  fs = 1e9/9;
  
  % .ref_wf_adc_idx = Reference receive channel, index into param.adcs,
  %   (surface location determined from this channel and all phase
  %   measurements made relative to it)
  param.ref_wf_adc_idx = 4;
  
  % .img = which waveform/adc pairs to load
  param.img = [1 2; 1 3; 1 4; 1 5; 1 6; 1 7];
  param.ref_wf_adc_idx = 4;
%   param.img = [1 2; 1 3; 1 9; 1 10];
%   param.ref_wf_adc_idx = 1;
%   param.img = [1 9; 1 10; 1 11; 1 12; 1 13; 1 14];
%   param.ref_wf_adc_idx = 4;
%   param.img = [2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 9; 2 10; 2 11; 2 12; 2 13; 2 14];
  
  % Correction coefficients
  %   param.td = zeros(size(param.img,1),1);
  %   param.amp = zeros(size(param.img,1),1);
  %   param.phase = zeros(size(param.img,1),1);
  param.td = 1e-9*[0 0 0 0 0 0 0 0 0 0 0 0].';
  param.td = 1e-9*[-83 -83 -83 -83 -83 -83 0 0 0 0 0 0].';
  param.amp = [0 0 0 0 0 0 0 0 0 0 0 0].';
  param.phase = [0 0 0 0 0 0 0 0 0 0 0 0].';

  if 0
    % fns = get_raw_files(struct('radar_name','mcrds','season_name','2011_Antarctica_TO'),'20111201_04_001')
    % fns = get_raw_files(struct('radar_name','mcrds','season_name','2011_Antarctica_TO'),'20111201_05_006')
    
    % Base path of data (does not include seg directory)
    base_path = '/cresis/snfs1/data/MCoRDS/2011_Antarctica_TO/20111202';

    fn_prefix = 'mcords2'; % 
    
    % Which segment directory to load from, leave empty for no segment
    % directory
    seg = 'seg_04';
    
    acquisition_num = 0;
    
    % Surface low gain
    file_nums = [1:34];
    param.rbins = [70:330];
    pc_param.Tpd = 1e-6;
    
    file_nums = [190:195]; % 214
    param.rbins = [70:330];
    pc_param.Tpd = 1e-6;
    
    %param.rbins = []; % Internal layers low gain
    %param.rbins = []; % Internal layers high gain
    %param.rbins = [280:450]; % Shallow ice bottom high gain
    
    % .gps_fn = Optional GPS file name (leave empty to disable)
    param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2011_Antarctica_TO/gps_20111201.mat';
    utc_time_correction = 0;

  elseif 1
    % fns = get_raw_files(struct('radar_name','mcrds','season_name','2011_Antarctica_TO'),'20111218_02_012')
    
    % Base path of data (does not include seg directory)
    base_path = '/cresis/snfs1/data/MCoRDS/2011_Antarctica_TO/20111219';

    fn_prefix = 'mcords2'; % 
    
    % Which segment directory to load from, leave empty for no segment
    % directory
    seg = 'seg_02';
    
    acquisition_num = 0;
    
    % Surface low gain
    file_nums = [587:645];
    param.rbins = [70:330];
    pc_param.Tpd = 1e-6;
    
    %param.rbins = []; % Internal layers low gain
    %param.rbins = []; % Internal layers high gain
    %param.rbins = [280:450]; % Shallow ice bottom high gain
    
    % .gps_fn = Optional GPS file name (leave empty to disable)
    param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2011_Antarctica_TO/gps_20111218.mat';
    utc_time_correction = 0;
    
  end

  param.season_name = '2011_Antarctica_TO';
  param.radar_name = 'mcords2';
 
end

% param.gps_fn = '';

% =======================================================================
% User Settings: Usually not changed
% =======================================================================

% param.type sent to basic_motion_comp.m
param.mocomp_type = 1;

% Transmit weights sent to lever_arm_fh
param.tx_weights = [1 1 1 1 1 1];

% Map of adcs to receivers for each waveform, used with lever_arm_fh
param.rx_paths = {};
param.rx_paths{1} = [1 7 8 9 10 11 12 1 1 2 3 4 5 6 1 1];
param.rx_paths{2} = [1 7 8 9 10 11 12 1 1 2 3 4 5 6 1 1];


% lever_arm_fh = lever arm
param.lever_arm_fh = @lever_arm;

% .plot_en = flag to enable plots
param.plot_en = true;

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR are included in the estimate)
param.snr_threshold = 12;

% presums = Number of presums (coherent averaging) to do
presums = 25;

% .averaging_fh = method to average complex vectors when computing
%   recommended channel compensation (mean is ideal, but median
%   may be necessary to remove outliers)
param.averaging_fh = @mean;

% Specify if cross correlation should be used (required for finding td)
% or if peak finding should be used when comparing channels
param.cross_correlation_flag = 1;

% Combine channels for surface tracker?
param.combine_channels = true;

% Remove DC/coherent-noise
param.coherent_noise_removal = true;

% .ref_bins: reference bins around surface peak bin to use in correlation
param.ref_bins = [-3 3];
% .search_bins: neighboring bins to search for a peak in other channels
param.search_bins = [-5 5];
% .Mt: oversampling factor for correlation output
param.Mt = 100;

basic_rx_chan_equalization;

return;



