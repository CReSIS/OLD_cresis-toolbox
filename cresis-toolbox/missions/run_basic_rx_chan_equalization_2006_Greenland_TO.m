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
  
  % .noise_rbins,rlines: Use for noise power calculation
  param.noise_rlines = [1 inf];
  param.noise_rbins = [600 700];
end

% Specify pulse compression properties for this waveform
pc_param.tukey = 0.0;
pc_param.window_func = @hanning;

radar_name = 'mcrds';
if strcmpi(radar_name,'mcrds')
  % Sampling frequency of radar (required to read data in)
  fs = 120e6;
  
  % .ref_wf_adc_idx = Reference receive channel, index into param.adcs,
  %   (surface location determined from this channel and all phase
  %   measurements made relative to it)
  param.ref_wf_adc_idx = 3;
  
  
  % Correction coefficients
  %   param.td = zeros(size(param.img,1),1);
  %   param.amp = zeros(size(param.img,1),1);
  %   param.phase = zeros(size(param.img,1),1);
  param.td = 1e-9*[0 0 0 0 0].';
  param.amp = [0 0 0 0 0].';
  param.phase = [0 0 0 0 0].';
  
  
  if 0
    % fns = get_raw_files(struct('radar_name','mcrds','season_name','2006_Greenland_TO'),'20060530_08_001')
    
    % Base path of data (does not include seg directory)
    base_path = '/cresis/snfs1/data/MCRDS/2006_Greenland/20060530b/';
    fn_prefix = 'data.';

    % Which segment directory to load from, leave empty for no segment
    % directory
    seg = '';
    
    file_adc_folder_name = '';
    
    % File index in filename
    file_nums = [63:73];
    file_nums = [379:393];
    
    % .gps_fn = Optional GPS file name (leave empty to disable)
    param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2006_Greenland_TO/gps_20060530.mat';
    
    pc_param.f0 = 140e6;
    pc_param.f1 = 160e6;
    pc_param.Tpd = 3e-6;
    % .img = which waveform/adc pairs to load
    param.img = [1 1; 1 2; 1 3; 1 4; 1 5];
  
    % .rbins = Range bins to search for surface in
    %   These range bins are post presumming
    param.rbins = [50:200];
    
  elseif 0
    % Base path of data (does not include seg directory)
    base_path = '/cresis/snfs1/data/MCRDS/2006_Greenland/20060528/';
    fn_prefix = 'data.';
    
    % Which segment directory to load from, leave empty for no segment
    % directory
    seg = '';
    
    file_adc_folder_name = '';
    
    % File index in filename
    file_nums = [430];
    
    % .gps_fn = Optional GPS file name (leave empty to disable)
    param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2006_Greenland_TO/gps_20060528.mat';
    
    pc_param.f0 = 140e6;
    pc_param.f1 = 160e6;
    pc_param.Tpd = 3e-6;
    % .img = which waveform/adc pairs to load
    param.img = [2 1; 2 2; 2 3; 2 4; 2 5];
  
    % .rbins = Range bins to search for surface in
    %   These range bins are post presumming
    param.rbins = [250:350];
    
  else
    % Base path of data (does not include seg directory)
    base_path = '/cresis/snfs1/data/MCRDS/2006_Greenland/20060528/';
    fn_prefix = 'data.';
    
    % Which segment directory to load from, leave empty for no segment
    % directory
    seg = '';
    
    file_adc_folder_name = '';
    
    % File index in filename
    file_nums = [406:430];
    %file_nums = [428:432]; % Shows roll
    
    % .gps_fn = Optional GPS file name (leave empty to disable)
    param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2006_Greenland_TO/gps_20060528.mat';
    
    pc_param.f0 = 140e6;
    pc_param.f1 = 160e6;
    pc_param.Tpd = 3e-6;
    % .img = which waveform/adc pairs to load
    param.img = [2 1; 2 2; 2 3; 2 4; 2 5];
  
    % .rbins = Range bins to search for surface in
    %   These range bins are post presumming
    param.rbins = [250:350];
    
  end

  param.season_name = '2006_Greenland_TO';
  param.radar_name = 'mcrds';
 
end

% param.gps_fn = '';
utc_time_correction = 0;

% =======================================================================
% User Settings: Usually not changed
% =======================================================================

% param.type sent to motion_comp.m
param.mocomp_type = 4;

% Transmit weights sent to lever_arm_fh
param.tx_weights = [1 1 1 1 1];

% Map of adcs to receivers for each waveform, used with lever_arm_fh
param.rx_paths = {};
param.rx_paths{1} = [1 2 3 4 5];
param.rx_paths{2} = [1 2 3 4 5];


% lever_arm_fh = lever arm
param.lever_arm_fh = @lever_arm;

% .plot_en = flag to enable plots
param.plot_en = true;

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR are included in the estimate)
param.snr_threshold = 12;

% presums = Number of presums (coherent averaging) to do
presums = 20;

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
param.ref_bins = [-3 3];
% .search_bins: neighboring bins to search for a peak in other channels
param.search_bins = [-2 2];
% .Mt: oversampling factor for correlation output
param.Mt = 100;

basic_rx_chan_equalization;

return;



