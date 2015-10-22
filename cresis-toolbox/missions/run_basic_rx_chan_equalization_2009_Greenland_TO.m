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
  param.noise_rbins = [1000 1200];
end

% Specify pulse compression properties for this waveform
pc_param.f0 = 140e6;
pc_param.f1 = 160e6;
pc_param.Tpd = 3e-6;
% pc_param.Tpd = 10e-6;
pc_param.tukey = 0.0;
pc_param.window_func = @hanning;

% Sampling frequency of radar (required to read data in)
fs = 120e6;

% .ref_wf_adc_idx = Reference receive channel, index into param.adcs,
%   (surface location determined from this channel and all phase
%   measurements made relative to it)
param.ref_wf_adc_idx = 4;

% Correction coefficients
%   param.td = zeros(size(param.img,1),1);
%   param.amp = zeros(size(param.img,1),1);
%   param.phase = zeros(size(param.img,1),1);
param.td = 1e-9*[0 0 0 0 0 0].';
param.amp = [0 0 0 0 0 0].';
param.phase = [0 0 0 0 0 0].';

if 1
  % Waveform 1
  param.phase = [-0.033333333	31.56666667	11.03333333	0	21.56666667	-13].';
  param.amp = [-0.133333333	-0.4	-0.233333333	0	-0.1	-0.566666667].';

  % .img = which waveform/adc pairs to load
  param.img = [1 1; 1 2; 1 3; 1 4; 1 5; 1 6]; % Low gain

  
else
  % Waveform 2
%   param.phase = [18 43.4 13.8 0 16.3 -16.3].';
%   param.amp = [0 -0.1 -0.1 0 0.7 -1.1].';
  param.phase = [-0.033333333	31.56666667	11.03333333	0	21.56666667	-13].';
  param.amp = [-0.133333333	-0.4	-0.233333333	0	-0.1	-0.566666667].';
  
  % .img = which waveform/adc pairs to load
  param.img = [2 1; 2 2; 2 3; 2 4; 2 5; 2 6]; % High gain
end

param.season_name = '2009_Greenland_TO';
param.radar_name = 'mcrds';

if 0
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090331_10_008')
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090331_05_010')
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090331_05_001')
  
  % Base path of data (does not include seg directory)
  base_path = '/cresis/data1/MCRDS/2009_Greenland/20090331/dataGISMO/';
  
  fn_prefix = 'data.200903311'; % 20090331_10_008
  fn_prefix = 'data.2009033107'; % 20090331_05_001, 20090331_05_010
  
  % Which segment directory to load from, leave empty for no segment
  % directory
  seg = '';
  
  file_adc_folder_name = '';
  
  % File index in filename
  file_nums = [324:349]; % Surface
  param.rbins = [70:330]; % Surface low gain
  %param.rbins = []; % Internal layers low gain
  
  %file_nums = [333:349]; % Bottom
  %param.rbins = [280:450]; % Shallow ice bottom high gain
  
  %     file_nums = [337 340:341]; % Internal layers high gain
  %     param.rbins = [300:320]; % Internal layers high gain
  
  file_nums = [213:252]; % 20090331_05_010: Surface
  param.rbins = [70:330]; % Surface low gain
  
  %     file_nums = [91:113]; % 20090331_05_001: Surface
  %     param.rbins = [70:330]; % Surface low gain
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2009_Greenland_TO/gps_20090331.mat';
  utc_time_correction = 0;
  
elseif 0
  
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090402_02_028')
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090402_02_032')
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090402_02_030')
  
  % Base path of data (does not include seg directory)
  base_path = '/cresis/snfs1/data/MCRDS/2009_Greenland/20090402/dataGISMO/';
  
  fn_prefix = 'data.';
  
  % Which segment directory to load from, leave empty for no segment
  % directory
  seg = '';
  
  file_adc_folder_name = '';
  
  % File index in filename
  file_nums = [515:530]; % 028: Surface
  file_nums = [515:519]; % 028: Bottom shallow
  %     param.rbins = [90:330]; % Surface low gain
  %     param.rbins = [445:600]; % Shallow ice bottom high gain
  
  %     file_nums = [616:630 655:658]; % 032
  %     file_nums = [616:618]; % 032: Bottom (limited by clutter/SNR)
  file_nums = [657]; % 032: Surface and Bottom
  %     param.rbins = [90:330]; % Surface low gain
  %     param.rbins = [440:600]; % Shallow ice bottom high gain 616-618
  param.rbins = [480:600]; % Shallow ice bottom high gain 655-658
  
  %     file_nums = [566:588 602:608]; % 032
  %     param.rbins = [90:330]; % Surface low gain
  %     file_nums = [566:569]; % 032
  %     param.rbins = [450:600]; % Ice bottom high gain 616-618
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2009_Greenland_TO/gps_20090402.mat';
  utc_time_correction = 3;
  
elseif 0
  % Helheim Glacier
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090409_02_025')
  
  % Base path of data (does not include seg directory)
  base_path = '/cresis/snfs1/data/MCRDS/2009_Greenland/20090409/dataGISMO/';
  
  fn_prefix = 'data.';
  
  % Which segment directory to load from, leave empty for no segment
  % directory
  seg = '';
  
  file_adc_folder_name = '';
  
  % File index in filename
  file_nums = [668:701]; % 002: Surface
  param.rbins = [75:330]; % Surface low gain
  %     param.rbins = [445:600]; % Shallow ice bottom high gain
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2009_Greenland_TO/gps_20090409.mat';
  
  % The radar time has jumps for this segment, so we require special utc time correction
  utc_time_correction = ct_filename_support(setfield(param,'day_seg','20090409_02'),'','records');
    
elseif 1
  % Helheim Glacier
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090411_01_002')
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090411_01_003')
  
  % Base path of data (does not include seg directory)
  base_path = '/cresis/snfs1/data/MCRDS/2009_Greenland/20090411/dataGISMO/';
  
  fn_prefix = 'data.';
  
  % Which segment directory to load from, leave empty for no segment
  % directory
  seg = '';
  
  file_adc_folder_name = '';
  
  if 0
  % File index in filename
    file_nums = [57:83]; % 002: Surface
    param.rbins = [75:330]; % Surface low gain
  elseif 0
    file_nums = [57:70 73:83]; % 002: Bottom (71/72 bad)
    param.rbins = [300:600]; % Ice bottom high gain
  else
    file_nums = [90:120]; % 002: Surface
    param.rbins = [75:330]; % Surface low gain
  end
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2009_Greenland_TO/gps_20090411.mat';
  
  % The radar time has jumps for this segment, so we require special utc time correction
  utc_time_correction = ct_filename_support(setfield(param,'day_seg','20090411_01'),'','records');
  
    
elseif 0
  % Transit
  % fns = get_raw_files(struct('radar_name','mcrds','season_name','2009_Greenland_TO'),'20090407_02_006')
  
  % Base path of data (does not include seg directory)
  base_path = '/cresis/snfs1/data/MCRDS/2009_Greenland/20090407/dataGISMO/';
  
  fn_prefix = 'data.';
  
  % Which segment directory to load from, leave empty for no segment
  % directory
  seg = '';
  
  file_adc_folder_name = '';
  
  if 0
    % Surface low gain
    file_nums = [194:196];
    param.rbins = [75:330];
  elseif 0
    % Low gain internal layers
    file_nums = [194:196];
    param.rbins = [380:440];
  else
    % High gain internal layers
    file_nums = [194:196];
    param.rbins = [380:440];
  end
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  param.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2009_Greenland_TO/gps_20090407.mat';
  
  % The radar time has jumps for this segment, so we require special utc time correction
  utc_time_correction = 0;%ct_filename_support(setfield(param,'day_seg','20090407_02'),'','records');
  
end

radar_name = param.radar_name;

% param.gps_fn = '';

% =======================================================================
% User Settings: Usually not changed
% =======================================================================

% param.type sent to motion_comp.m
param.mocomp_type = 4;

% Transmit weights sent to lever_arm_fh
param.tx_weights = [1 1 1 1 1 1];

% Map of adcs to receivers for each waveform, used with lever_arm_fh
param.rx_paths = {};
param.rx_paths{1} = [1 2 3 4 5 6];
param.rx_paths{2} = [1 2 3 4 5 6];


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
param.search_bins = [-2 2];
% .Mt: oversampling factor for correlation output
param.Mt = 100;

basic_rx_chan_equalization;

return;



