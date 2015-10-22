% CReSIS CSARP-compatible radar simulator run script
%
% Author: Logan Smith
%

clear;% clc;

global gRadar
if isempty(gRadar)
  startup
  %   error('gRadar has been cleared!  Run startup.m to recreate!')
end

param.debug_level = 1;
param.validation_mode_en = 1; % Power received = power transmitted
param.visible_surf = 1;

% =========================================================================
% Output Paths
% =========================================================================
param.day_seg = '20110708_01';
fprintf('day_seg: %s\n',param.day_seg)
% radar data output file path
param.data_path = ...
  '/cresis/scratch1/jliwestc/mdce/mcords/mcords_simulator/simulated_data/motion/20110708_01';
% param.data_path = ...
%   '/cresis/scratch2/mdce/mcords/mcords_simulator/simulated_data/idealcase/20110708_01/';
param.data_fn = 'mcords.rec000.r3-%d.20110708000000.0001.dat';
fprintf('Data Path: %s\n',fullfile(param.data_path,param.data_fn))
% gps file path 
param.gps_path = ...
  '/cresis/scratch1/jliwestc/mdce/csarp_support/gps/mcords_simulator/';
% param.gps_path = ...
%   '/cresis/projects/dev/csarp_support/gps/mcords_simulator/';
param.gps_fn = 'gps_20110708.mat';
fprintf('GPS Path: %s\n',fullfile(param.gps_path,param.gps_fn))

% =========================================================================
% Simulation Parameters
% =========================================================================

% RADAR PARAMETERS
param.radar.fc = 150e6; % center frequency (Hz)
param.radar.fs = 120e6; % sampling frequency (Hz)
param.radar.prf = 10000/64; % pulse repetition frequency (Hz)
param.radar.Tpri  = 1/param.radar.prf; % pulse repetition interval (s)
param.radar.Vpp = 2; % max ADC peak-to-peak input voltage (V)
param.radar.Z0 = 1; % system characteristic impedance (Ohms)
param.radar.beamwidth_deg = 40; % antenna full-beamwidth (deg)
param.radar.tx_amp = .5; 

% WAVEFORM PARAMETERS
param.radar.wfs{1}.presums = 1; % number of hardware presums
param.radar.wfs{1}.f0 = 140e6; % chirp start frequency (Hz)
param.radar.wfs{1}.f1 = 160e6; % chirp stop frequency (Hz)
param.radar.wfs{1}.BW = param.radar.wfs{1}.f1-param.radar.wfs{1}.f0; % chirp bandwidth (Hz)
param.radar.wfs{1}.Tpd = 3e-6; % chirp duration (s)
param.radar.wfs{1}.t0 = -0e-5; % sample delay (s)
param.radar.rx_gain = [1 1 1 1]; % receiver gains (V)
param.radar.tuk_ratio = 0.25;

% ADC PARAMETERS
param.radar.adc_bits = 14; % number of ADC bits
% ADC dynamic range (dB)
param.radar.adc_dynamic_gain = 20*log10(2^param.radar.adc_bits);

param.radar.tx_weights = [1 1 1 1];

% PLATFORM PARAMETERS
param.platform.h = 500; % platform height above air/ice interface (m)
param.platform.v = 60; % platform velocity (m/s)
param.platform.lever_arm_fh = ...
  @lever_arm_mcords_simulator; % antenna lever arm function handle
param.synth_gps = 0;  % 1 for fake gps and ins data; 0 for real gps data
if ~param.synth_gps
    param.real_gps_fn = '/cresis/projects/dev/csarp_support/gps/2011_Greenland_TO/gps_20110428.mat';
end

param.radar.chans = 4; % number of receive elements
param.radar.rx_path = 1:param.radar.chans; % vector of receive elements

% INTRODUCED ERRORS
% param.errors.td = [0 2 4 6]*1e-9; % introduced time delays (s)
% param.errors.phase = [0 10 20 30]; % introduced phase param.errors (deg)
% param.errors.mag = [0 1 2 3]; % introduced amplitude param.errors (dB)
param.errors.td = [0 0 0 0]*1e-9; % introduced time delays (s)
param.errors.phase = [0 0 0 0]; % introduced phase param.errors (deg)
param.errors.mag = [0 0 0 0]; % introduced amplitude param.errors (dB)
% complex channel equalization errors
param.errors.chan_equal = 10.^(param.errors.mag/20).*exp(1j*param.errors.phase/180*pi);
% additional noise power (dB)
param.errors.noise_floor = [0 0 0 0];

% Target locations
% param.target.x = [0 0 0 0 0];%100 -200 -250 300];
% param.target.z = [500 1000 1500 2000 3000];
param.target.x = [0];
param.target.z = [500];
param.ntarget = length(param.target.x);
param.target.Z0 = 1800;           % target area in depth is within [Zc-Z0,Zc+Z0]

param.scene.surface_angle = deg2rad(0);
% =========================================================================

% =========================================================================
% Run Simulator
% =========================================================================
csarp_simulator(param)


