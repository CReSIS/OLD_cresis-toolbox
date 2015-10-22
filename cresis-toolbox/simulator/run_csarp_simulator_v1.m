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
param.validation_mode_en = 0; % Power received = power transmitted
param.visible_surf = 0;       % Do not generate surface signal using same return power as target

% =========================================================================
% Output Paths
% =========================================================================
param.day_seg = '20060530_01';
fprintf('day_seg: %s\n',param.day_seg)

% radar data output file path
param.data_path = ...
  '/cresis/scratch1/petan/mdce/mcords/mcords_simulator/simulated_data/test_run/20060530_01';
%   '/cresis/scratch1/petan/mdce/mcords/mcords_simulator/simulated_data/motion/20060530_01';

param.data_fn = 'mcords.rec000.r3-%d.20110708000000.0001.dat';
fprintf('Data Path: %s\n',fullfile(param.data_path,param.data_fn))

% gps file path 
param.gps_path = ...
  '/cresis/scratch1/petan/mdce/csarp_support/gps/mcords_simulator/';
param.gps_fn = 'gps_20060530.mat';
fprintf('GPS Path: %s\n',fullfile(param.gps_path,param.gps_fn))

% =========================================================================
% Simulation Parameters
% =========================================================================

% RADAR PARAMETERS
param.radar.fc = 150e6;                     % center frequency (Hz)
param.radar.fs = 120e6;                     % sampling frequency (Hz)
param.radar.prf = 10000;%10000/64;          % pulse repetition frequency (Hz)
param.radar.Tpri  = 1/param.radar.prf;      % pulse repetition interval (s)
param.radar.Vpp = 2;                        % max ADC peak-to-peak input voltage (V)
param.radar.Z0 = 1;                         % system characteristic impedance (Ohms)
param.radar.beamwidth_deg = 10;%20;         % antenna full-beamwidth (deg)
param.radar.tx_pwr = 800;                   % antenna total TX power for all channels
param.radar.Ls = 10^(2/10);                 % Radar System loss

% WAVEFORM PARAMETERS
param.radar.wfs{1}.presums = 64;  % number of hardware presums
param.radar.wfs{1}.f0 = 140e6;    % chirp start frequency (Hz)
param.radar.wfs{1}.f1 = 160e6;    % chirp stop frequency (Hz)
param.radar.wfs{1}.BW = param.radar.wfs{1}.f1-param.radar.wfs{1}.f0; % chirp bandwidth (Hz)
param.radar.wfs{1}.Tpd = 3e-6; % chirp duration (s)
param.radar.wfs{1}.t0 = 0e-5; % sample delay (s)
param.radar.rx_gain = [1 1 1 1]; % receiver gains (V)
param.radar.tuk_ratio = 0.20;%0.25;

% ADC PARAMETERS
param.radar.adc_bits = 14; % number of ADC bits
% ADC dynamic range (dB)
param.radar.adc_dynamic_gain = 20*log10(2^param.radar.adc_bits);

param.radar.tx_weights = [1 1 1 1];

% PLATFORM PARAMETERS
param.platform.h = 500; % platform height above air/ice interface (m)
param.platform.v = 60; % platform velocity (m/s)
param.platform.lever_arm_fh = @lever_arm_mcords_simulator; % antenna lever arm function handle
param.synth_gps = 1;  % 1 for fake gps and ins data; 0 for real gps data
if ~param.synth_gps
    param.real_gps_fn = '/cresis/projects/dev/csarp_support/gps/2006_Greenland_TO/gps_20060530.mat';
end

param.radar.chans = 4;%4;                     % number of receive elements
param.radar.rx_path = 1:param.radar.chans;    % vector of receive elements
param.Ant_Gain = 10.0;                        % Antenna Gain is 10 dB at mainlobe

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

% Target Parameters
param.target.x = [0 0 0 0 0];                 
param.target.z = [0 1000 1100 1200 1500];
param.ntarget = length(param.target.x);
param.target.Z0 = 1500;                       % target area in depth is within [Zc-Z0,Zc+Z0]

param.scene.surface_angle = deg2rad(0);

% Clutter Parameters
param.clutter_rcs = 0.01;                     % Average Clutter RCS is -20 dB
param.Gain_clutter = param.Ant_Gain/sqrt(2);  % Antenna Gain at sidelobe for clutter

% Ice parameter
param.ice_atten = 0.02;                       % [dB/m]2-way ice attenuation factor due to absorption/scattering

% =========================================================================

% =========================================================================
% Run Simulator
% =========================================================================
csarp_simulator_v1(param)


