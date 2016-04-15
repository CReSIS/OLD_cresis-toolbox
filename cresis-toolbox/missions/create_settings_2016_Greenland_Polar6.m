% script create_settings_2015_Greenland_Polar6
%
% Creates NI radar depth sounder settings
%
% See link_budgets.xls

physical_constants; % c = speed of light

% Define waveforms
base_dir = 'E:\waveforms\';

f0_list = [150e6 180e6];
f1_list = [520e6 210e6];
DDC_select_list = [1 2]; % Which DDC mode to use
cal_settings = [1 2];
prf = 10000;

% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
velocity = 110; % m/s
fudge_factor = [1.5 1]; % Set to 1 for lambda/4, more than 1 to relax the requirement
presums = round(fudge_factor .* c./abs(f1_list+f0_list)/2 / velocity * prf / 2)*2

final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 0
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0 0 0 0 0 0 0 0];
    final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
    final_DDS_amp{idx} = [4000 4000 4000 4000 4000 4000 4000 4000];
    final_DDS_time{idx} =  [0 0 0 0 0 0 0 0];
  end
else
  % After transmit calibration during 20160401 test flight
  %   From basic_tx_chan_equalization_SEASON_NAME.m
  
  % NOTE: These values are valid for when DDS channels 5-8 come up one 1440
  % MHz clock cycle after channels 1-4.
  
  % 150-520 MHz
  final_DDS_phase{end+1} = [63.3	86.1	-16.5	0.0	8.9	-17.4	68.8	51.1];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [1312	2849	2657	3572	4000	2618	2574	1386];
  final_DDS_time{end+1} =  [-2.62	-2.35	-0.13	0.00	-0.56	-0.82	-3.30	-3.50];
    
  % 180-210 MHz
  final_DDS_phase{end+1} = [61.4	85.2	-15.1	0.0	8.7	-1.0	76.8	56.1];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [1172	2550	3026	3650	4000	3106	2361	1223];
  final_DDS_time{end+1} =  [-2.62	-2.35	-0.13	0.00	-0.56	-0.82	-3.30	-3.50];
  
end

% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = chebwin(8,30).';

% =========================================================================
% =========================================================================
% AUTOMATED SECTION
% =========================================================================
% =========================================================================

param = [];
param.season_name = '2016_Greenland_Polar6';
param.radar_name = 'rds';
param.gps_source = 'awi-final';
clear phase_centers;
for tx_chan = 1:8
  % Just enable the current antenna to determine its phase center
  tx_weights = zeros(1,8);
  tx_weights(tx_chan) = 1;
  rxchan = 12; % Fix the receiver (it should not matter which one you choose)
  % Determine phase center for the antenna
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, rxchan);
end
% Adjust phase centers to the mean phase center position
phase_centers = bsxfun(@minus,phase_centers,mean(phase_centers,2));

% Create waveforms directory if it does not exist
if ~exist(base_dir,'dir')
  mkdir(base_dir);
end

%% Survey Mode + loopback, noise, and deconv modes
% <3250 m thick ice, 900 to 1900 ft AGL
ice_thickness = [3250 3250];
for freq_idx = 1:length(f0_list)
  param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 500*12*2.54/100;
  param.tg.Haltitude = 1400*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 2 presums(freq_idx)-4];
  param.wfs(1).atten = 43;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
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
  % Loopback Mode without delay line
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 2000*12*2.54/100;
  param.tg.Haltitude = 0e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK_NO_DELAY.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Loopback Mode (10e-6 delay line)
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 2000*12*2.54/100;
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
  param.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
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


%% Survey Mode for thin ice
% <2500 m thick ice, 900 to 1900 ft AGL
ice_thickness = [2500 2500];
for freq_idx = 1:length(f0_list)
  param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2];
  param.tg.altitude_guard = 500*12*2.54/100;
  param.tg.Haltitude = 1400*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 presums(freq_idx)-2];
  param.wfs(1).atten = 43;
  param.wfs(2).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:2).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('thinice_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Image Mode
% Ice thickness 2000 m +/- 700 m, 6000 +/- 250 ft AGL
param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
param.max_data_rate = 755;
param.DDC_select = 1;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
altitude_agl_feet = 6000;
swath_beamwidth_deg = 45;
ice_thickness = 1800;
ice_thickness_guard = 500;
param.tg.staged_recording = false;
param.tg.altitude_guard = 250*12*2.54/100 + ice_thickness_guard*sqrt(er_ice);
param.tg.Haltitude = altitude_agl_feet*12*2.54/100 + ice_thickness*sqrt(er_ice);
param.tg.Hice_thick = altitude_agl_feet*12*2.54/100/cosd(swath_beamwidth_deg) ...
  + ice_thickness*sqrt(er_ice)/cosd(asind(sind(swath_beamwidth_deg)/sqrt(er_ice))) - (altitude_agl_feet*12*2.54/100+ice_thickness*sqrt(er_ice));
param.prf = prf;
param.presums = [ceil(presums(freq_idx)/4)*2 ceil(presums(freq_idx)/4)*2];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{1} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 10e-6;
param.wfs(2).Tpd = 10e-6;
param.wfs(1).phase = final_DDS_phase{1};
param.wfs(2).phase = final_DDS_phase{1};
beam_angle_deg = 20; % Positive to the left
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
param.wfs(1).delay = final_DDS_time{1} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(2).delay = final_DDS_time{1} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.f0 = f0_list(1);
param.f1 = f1_list(1);
param.DDC_freq = (param.f0+param.f1)/2;
[param.wfs(1:2).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,altitude_agl_feet,param.wfs(end).Tpd*1e6,ice_thickness));
write_cresis_xml(param);

% Pattern measurements
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.tg.altitude_guard = 500*12*2.54/100;
param.tg.Haltitude = 3000*12*2.54/100;
param.tg.Hice_thick = 0;
[param.wfs(1:2).atten] = deal(43);
param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_PATTERN.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6));
write_cresis_xml(param);

% Narrow band
param.wfs(1).Tpd = 10e-6;
param.wfs(2).Tpd = 10e-6;
altitude_agl_feet = 6000;
swath_beamwidth_deg = 45;
ice_thickness = 1800;
ice_thickness_guard = 1250;
param.presums = [ceil(presums(freq_idx)/4)*2 ceil(presums(freq_idx)/4)*2];
param.tg.staged_recording = false;
param.tg.altitude_guard = 250*12*2.54/100 + ice_thickness_guard*sqrt(er_ice);
param.tg.Haltitude = altitude_agl_feet*12*2.54/100 + ice_thickness;
param.tg.Hice_thick = ice_thickness/cosd(swath_beamwidth_deg) ...
  + ice_thickness*sqrt(er_ice)/cosd(asind(sind(swath_beamwidth_deg)/sqrt(er_ice))) - (ice_thickness+ice_thickness*sqrt(er_ice));
param.DDC_select = 2;
param.f0 = f0_list(2);
param.f1 = f1_list(2);
param.DDC_freq = (param.f0+param.f1)/2;
DDS_amp = final_DDS_amp{2} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.wfs(1).phase = final_DDS_phase{2};
param.wfs(2).phase = final_DDS_phase{2};
beam_angle_deg = 20; % Positive to the left
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
param.wfs(1).delay = final_DDS_time{2} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(2).delay = final_DDS_time{2} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,altitude_agl_feet,param.wfs(end).Tpd*1e6,ice_thickness));
write_cresis_xml(param);

% Pattern measurements
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.tg.altitude_guard = 500*12*2.54/100;
param.tg.Haltitude = 3000*12*2.54/100;
param.tg.Hice_thick = 0;
[param.wfs(1:2).atten] = deal(43);
param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_PATTERN.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6));
write_cresis_xml(param);

%% Equalization (Using Ocean)
% Haltitude +/- 1000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
Haltitude = [0 2000 0 3000];
Tpd_list = [1e-6 1e-6 3e-6 3e-6];
for Tpd_idx = 1:length(Tpd_list)
  Tpd = Tpd_list(Tpd_idx);
  for freq_idx = 1:length(f0_list)
    param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
    param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
    param.DDC_select = DDC_select_list(freq_idx);
    param.max_duty_cycle = 0.12;
    param.create_IQ = false;
    param.tg.staged_recording = false;
    param.tg.altitude_guard = 1000*12*2.54/100;
    param.tg.Haltitude = Haltitude(Tpd_idx)*12*2.54/100;
    param.tg.Hice_thick = 0;
    param.prf = prf;
    param.presums = [10 10 10 10 10 10 10 10 10];
    [param.wfs(1:8).atten] = deal(35);
    [param.wfs(9:9).atten] = deal(43);
    param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
    param.tukey = 0.08;
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

%% Max power mode with max frequency range (only for EMI survey)
if 0
  for freq_idx = [1] % 1:length(f0_list)
    param = struct('radar_name','mcords5','num_chan',24,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
    param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 4.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true;
    param.DDC_select = DDC_select_list(freq_idx);
    param.max_duty_cycle = 0.12;
    param.create_IQ = false;
    param.tg.staged_recording = false;
    param.tg.altitude_guard = 1000*12*2.54/100;
    param.tg.Haltitude = 1400*12*2.54/100;
    param.tg.Hice_thick = 0;
    param.prf = prf;
    param.presums = presums(freq_idx);
    param.wfs(1).atten = 43;
    DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
    param.tx_weights = DDS_amp;
    param.tukey = 0;
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
end