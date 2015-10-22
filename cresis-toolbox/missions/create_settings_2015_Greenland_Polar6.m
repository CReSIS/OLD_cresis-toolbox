% script create_settings_2015_Greenland_Polar6
%
% Creates NI radar depth sounder settings
%
% See link_budgets.xls

% Define waveforms
base_dir = 'D:\waveforms\';

DDC_mode_strings = {'Non','DDC4','DDC8'};

% [sqrt(0.3)*exp(1j*0.5*pi/180)
%             sqrt(0.50)*exp(1j*10*pi/180)
%             sqrt(0.93)*exp(1j*0.2*pi/180)
%             sqrt(1)*exp(1j*0*pi/180)
%             sqrt(1)*exp(1j*0*pi/180)
%             sqrt(0.93)*exp(1j*0.2*pi/180)
%             sqrt(0.5)*exp(1j*10*pi/180)
%             sqrt(0.3)*exp(1j*0.5*pi/180)];

final_DDS_phase = [];
final_DDS_amp = [];
final_DDS_time = [];
if 0
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:4
    final_DDS_phase{idx} = [0 0 0 0 0 0 0 0];
    final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
    final_DDS_amp{idx} = [3000 3000 3000 3000 3000 3000 3000 3000];
    final_DDS_time{idx} =  [0 0 0 0 0 0 0 0];
  end
else
  % After transmit calibration during Sept 11, 2015 test flight
  %   From basic_tx_chan_equalization_SEASON_NAME.m
  % 165-510
  final_DDS_phase{1} = [78.3	45.9	0.7	0.0	20.0	1.1	46.4	94.9];
  final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{1} = [931	2069	2173	2698	3000	2041	1870	1023];
  final_DDS_time{1} =  [-2.49	-2.70	0.05	0.00	0.19	0.02	-2.83	-2.52];
  
  % 170-230
  final_DDS_phase{2} = [86.4	128.0	125.1	0.0	-145.9	174.7	150.5	71.9];
  final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{2} = [796	2058	2507	2953	3000	2552	1952	914];
  final_DDS_time{2} =  [-2.50	-1.88	1.77	0.00	-1.41	1.77	-2.16	-3.30];
    
  % 180-210
  final_DDS_phase{3} = [83.5	126.2	128.8	0.0	-55.1	-173.7	148.4	71.3];
  final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{3} = [743	1901	2480	3000	2853	2493	1900	876];
  final_DDS_time{3} =  [-2.50	-1.88	1.77	0.00	-1.41	1.77	-2.16	-3.30];
    
  % 150-600
  final_DDS_phase{4} = [0 0 0 0 0 0 0 0];
  final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{4} = [3000 3000 3000 3000 3000 3000 3000 3000];
  final_DDS_time{4} =  [0 0 0 0 0 0 0 0];
  
end

f0 = 165e6;
f1 = 510e6;

% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = chebwin(8,30).';

physical_constants;
param = [];
param.season_name = '2015_Greenland_Polar6';
param.radar_name = 'rds';
param.gps_source = 'awi-final';
clear phase_centers;
for tx_chan = 1:8
  % Just enable the current antenna to determine its phase center
  tx_weights = zeros(1,8);
  tx_weights(tx_chan) = 1;
  % Determine phase center for the antenna
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, tx_chan);
end

% Create waveforms directory if it does not exist
if ~exist(base_dir,'dir')
  mkdir(base_dir);
end

%% Survey Mode + loopback, noise, and deconv modes
% <3250 m thick ice, 900 to 1900 ft AGL
f0_list = [f0 170e6 180e6 150e6];
f1_list = [f1 230e6 210e6 600e6];
DDC_select_list = [1 2 2 0];
cal_settings = [1 2 3 1];
presums = [36 36 36 64]
ice_thickness = [3250 3250 3250 2500]
for freq_idx = 1:length(f0_list)
  param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'fs',1600e6,'fs_sync',62.5e6,'fs_dds',1600e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
  param.max_tx = [3000 3000 3000 3000 3000 3000 3000 3000]; param.max_data_rate = 750; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 500*12*2.54/100;
  param.tg.Haltitude = 1400*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = 12000;
  param.presums = [2 2 presums(freq_idx)-4];
  param.wfs(1).atten = 43;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.1;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
  if freq_idx == 1
    % Default Mode
    param.fn = fullfile(base_dir,'default.xml');
    write_cresis_xml(param);
  end
  % Loopback Mode
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 10e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Deconvolution Mode (for over calm lake or sea ice lead)
  param.wfs(1).atten = 43;
  param.wfs(2).atten = 43;
  param.wfs(3).atten = 43;
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 2000*12*2.54/100;
  param.tg.Haltitude = 3000*12*2.54/100;
  param.tg.Hice_thick = 0 * 12*2.54/100/1.78;
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  if freq_idx == 1
    % Noise Mode
    param.tx_weights = [0 0 0 0 0 0 0 0];
    [param.wfs(1:1).tx_mask] = deal([1 1 1 1 1 1 1 1]);
    param.wfs(1).atten = 43;
    param.wfs(2).atten = 0;
    param.wfs(3).atten = 0;
    param.tg.staged_recording = [1 2 3];
    param.tg.altitude_guard = 500*12*2.54/100;
    param.tg.Haltitude = 1400*12*2.54/100;
    param.tg.Hice_thick = 3250;
    param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_NOISE.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
    write_cresis_xml(param);
  end
end

%% Image Mode
% Ice thickness 2000 m +/- 700 m, 6000 +/- 250 ft AGL
param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'fs',1600e6,'fs_sync',62.5e6,'fs_dds',1600e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
param.max_tx = [3000 3000 3000 3000 3000 3000 3000 3000]; param.max_data_rate = 700; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
param.max_data_rate = 755;
param.DDC_select = 1;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
altitude_agl_feet = 6000;
swath_beamwidth_deg = 45;
ice_thickness = 1800;
ice_thickness_guard = 750;
param.tg.staged_recording = false;
param.tg.altitude_guard = 250*12*2.54/100 + ice_thickness_guard*1.78;
param.tg.Haltitude = altitude_agl_feet*12*2.54/100 + ice_thickness;
param.tg.Hice_thick = ice_thickness/cosd(swath_beamwidth_deg) ...
  + ice_thickness*1.78/cosd(asind(sind(swath_beamwidth_deg)/1.78)) - (ice_thickness+ice_thickness*1.78);
param.prf = 12000;
param.presums = [24 24];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{1} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{1};
param.wfs(2).phase = final_DDS_phase{1};
param.wfs(1).delay = final_DDS_time{1} - (0.468 / (c/2) * sind(20)*(0:7))*1e9;
param.wfs(2).delay = final_DDS_time{1} + (0.468 / (c/2) * sind(20)*(0:7))*1e9;
param.f0 = f0;
param.f1 = f1;
param.DDC_freq = (param.f0+param.f1)/2;
[param.wfs(1:2).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,altitude_agl_feet,param.wfs(end).Tpd*1e6,ice_thickness));
write_cresis_xml(param);

% Pattern measurements
param.tg.altitude_guard = 500*12*2.54/100;
param.tg.Haltitude = 3000*12*2.54/100;
param.tg.Hice_thick = 0;
[param.wfs(1:2).atten] = deal(43);
param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_PATTERN.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6));
write_cresis_xml(param);

% Narrow band
altitude_agl_feet = 6000;
swath_beamwidth_deg = 45;
ice_thickness = 1800;
ice_thickness_guard = 1250;
param.presums = [28 28];
param.tg.staged_recording = false;
param.tg.altitude_guard = 250*12*2.54/100 + ice_thickness_guard*1.78;
param.tg.Haltitude = altitude_agl_feet*12*2.54/100 + ice_thickness;
param.tg.Hice_thick = ice_thickness/cosd(swath_beamwidth_deg) ...
  + ice_thickness*1.78/cosd(asind(sind(swath_beamwidth_deg)/1.78)) - (ice_thickness+ice_thickness*1.78);
param.DDC_select = 2;
param.f0 = 180e6;
param.f1 = 210e6;
param.DDC_freq = (param.f0+param.f1)/2;
DDS_amp = final_DDS_amp{3} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.wfs(1).phase = final_DDS_phase{3};
param.wfs(2).phase = final_DDS_phase{3};
param.wfs(1).delay = final_DDS_time{3} - (0.468 / (c/2) * sind(20)*(0:7))*1e9;
param.wfs(2).delay = final_DDS_time{3} + (0.468 / (c/2) * sind(20)*(0:7))*1e9;
param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,altitude_agl_feet,param.wfs(end).Tpd*1e6,ice_thickness));
write_cresis_xml(param);

%% Equalization High Altitude (Using Ocean)
% 7500 +/- 1000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
Haltitude = [1500 3000];
Tpd_list = [1e-6 3e-6];
for Tpd_idx = 1:length(Tpd_list)
  Tpd = Tpd_list(Tpd_idx);
  f0_list = [f0 170e6 180e6 150e6];
  f1_list = [f1 230e6 210e6 600e6];
  DDC_select_list = [1 2 2 0];
  cal_settings = [1 2 3 1];
  for freq_idx = 1:length(f0_list)
    param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'fs',1600e6,'fs_sync',62.5e6,'fs_dds',1600e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
    param.max_tx = [3000 3000 3000 3000 3000 3000 3000 3000]; param.max_data_rate = 700; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
    param.DDC_select = DDC_select_list(freq_idx);
    param.max_duty_cycle = 0.12;
    param.create_IQ = false;
    param.tg.staged_recording = false;
    param.tg.altitude_guard = 500*12*2.54/100;
    param.tg.Haltitude = Haltitude(Tpd_idx)*12*2.54/100;
    param.tg.Hice_thick = 0;
    param.prf = 12000;
    param.presums = [10 10 10 10 10 10 10 10 10];
    [param.wfs(1:8).atten] = deal(35);
    [param.wfs(9:9).atten] = deal(43);
    param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
    param.tukey = 0.1;
    param.Tpd = Tpd;
    for wf=1:9
      param.wfs(wf).phase = final_DDS_phase{cal_settings(freq_idx)};
    end
    param.delay = final_DDS_time{cal_settings(freq_idx)};
    param.f0 = f0_list(freq_idx);
    param.f1 = f1_list(freq_idx);
    param.DDC_freq = (param.f0+param.f1)/2;
    for wf=1:8
      param.wfs(wf).tx_mask = ones(1,8);
      param.wfs(wf).tx_mask(9-wf) = 0;
    end
    for wf=9:9
      param.wfs(wf).tx_mask = [0 0 0 0 0 0 0 0];
    end
    param.fn = fullfile(base_dir,sprintf('txequal_%.0f-%.0fMHz_%.0fft_%.0fus.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.Tpd*1e6));
    write_cresis_xml(param);
  end
end

%% Max power mode with max frequency range (EMI survey)
f0_list = [150e6];
f1_list = [600e6];
DDC_select_list = [0];
cal_settings = [4];
for freq_idx = 1:length(f0_list)
  param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'fs',1600e6,'fs_sync',62.5e6,'fs_dds',1600e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
  param.max_tx = [3000 3000 3000 3000 3000 3000 3000 3000]; param.max_data_rate = 750; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 1400*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.prf = 12000;
  param.presums = [32];
  param.wfs(1).atten = 43;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.1;
  param.wfs(1).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  % RFI 0, 25, 50, 75, 100 % Power Modes
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)} * sqrt(0);
  param.fn = fullfile(base_dir,sprintf('RFI_%.0f-%.0fMHz_%.0fus_000p.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)} * sqrt(0.25);
  param.fn = fullfile(base_dir,sprintf('RFI_%.0f-%.0fMHz_%.0fus_025p.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)} * sqrt(0.5);
  param.fn = fullfile(base_dir,sprintf('RFI_%.0f-%.0fMHz_%.0fus_050p.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)} * sqrt(0.75);
  param.fn = fullfile(base_dir,sprintf('RFI_%.0f-%.0fMHz_%.0fus_075p.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)} * sqrt(1.00);
  param.fn = fullfile(base_dir,sprintf('RFI_%.0f-%.0fMHz_%.0fus_100p.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
end
