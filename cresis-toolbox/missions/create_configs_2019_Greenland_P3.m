% script create_configs_2019_Greenland_P3
%
% Creates NI radar depth sounder settings

physical_constants; % c = speed of light

% Define waveforms
if ispc
  base_dir = 'C:\waveforms\';
else
  base_dir = '/scratch/waveforms/';
end

f0_list = [180e6];
f1_list = [210e6];
DDC_select_list = [1]; % Which DDC mode to use
cal_settings = [1];
prf = 12000;

%% Calculate the desired presums
% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
physical_constants; % c = speed of light
velocity = [140]; % m/s
quarter_wavelength = c./max(abs(f0_list),abs(f1_list)) / 4;
presums = floor(quarter_wavelength ./ velocity * prf / 2) * 2;
fprintf('For prf %g Hz, max presums for quarter wavelength (%g m) sampling is: %g presums\n', prf, quarter_wavelength, presums);

final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 1
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0 0 0 0 0 0 0 0];
    final_DDS_phase_no_time{idx} = [0 0 0 0 0 0 0 0]; % not used usually
    final_DDS_amp{idx} = [40000 40000 40000 40000 40000 40000 40000 0];
    final_DDS_time{idx} =  [0 0 0 0 0 0 0 0];
  end
  % final_tx_mask: Listed backwards antenna 8 to 1
  final_tx_mask = [1 0 0 0 0 0 0 0]; 
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % Tx chan equalization txequal_mcords3_20170307_163925_00.xml
  idx = 1;
  final_DDS_phase{idx} = [131.5	-35.1	0.0	-158.4	113.8	-73.0	139.6	0.0];
  final_DDS_phase_no_time{idx} = [18.6	-8.6	0.0	97.1	-70.1	13.7	56.8	0.0]; % not used usually
  final_DDS_amp{idx} = [21000   25000   32500   40000   30000   27000   25000  0];
  final_DDS_time{idx} =  [-8.42	-0.16	0.00	1.38	-3.50	-0.75	-3.72	0.00];
  
  % final_tx_mask: Listed backwards antenna 8 to 1
  final_tx_mask = [1 0 0 0 0 0 0 0]; 
end

Hwindow_orig = [1 1 1 1 1 1 1 0]; % Desired window created during transmit calibration

physical_constants;
param = [];
param.season_name = '2019_Greenland_P3';
param.radar_name = 'rds';
param.gps_source = 'atm-final20160301';
clear phase_centers;
for tx_chan = 1:7
  tx_weights = zeros(1,7);
  tx_weights(tx_chan) = 1;
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, tx_chan);
end

% Create waveforms directories if they do not exist
if ~exist(base_dir,'dir')
  mkdir(base_dir);
end
calval_dir = fullfile(base_dir, 'calval');
if ~exist(calval_dir,'dir')
  mkdir(calval_dir);
end

%% Survey Mode + loopback, noise, and deconv modes
% <4000 m thick ice, 1250 +/- 750 ft AGL
ice_thickness = [3500];
freq_idx = 1;
param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = ice_thickness(freq_idx);
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_3500mthick.xml');
param.prf = 12000;
param.presums = [3 3 29];
param.wfs(1).atten = 33;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.delay = final_DDS_time{cal_settings(freq_idx)};
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
[param.wfs(:).tx_mask] = deal(final_tx_mask);
write_cresis_xml(param);
% Default Mode
param.fn = fullfile(base_dir,'default.xml');
write_cresis_xml(param);
% TTL Settings Check
param.tx_weights(2:end) = 0;
param.tukey = 0.01;
param.fn = fullfile(calval_dir,'survey_mode_10us_3wf_3500mthick_TTL_CHECK.xml');
write_cresis_xml(param);
% Noise Mode
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_3500mthick_NOISE.xml');
param.tx_weights = DDS_amp;
[param.wfs(:).tx_mask] = deal([1 1 1 1 1 1 1 1]);
write_cresis_xml(param);
% Loopback Mode
param.tg.staged_recording = false;
[param.wfs(:).tx_mask] = deal(final_tx_mask);
param.tg.altitude_guard = 0*12*2.54/100;
param.tg.Haltitude = 0*12*2.54/100;
param.tg.Hice_thick = 20e-6 * 3e8/2/sqrt(3.15); % Long enough for delay line
param.fn = fullfile(calval_dir,'survey_mode_10us_3wf_3500mthick_LOOPBACK.xml');
write_cresis_xml(param);
% Deconvolution Mode (for Sea Ice)
param.prf = 8500;
param.tg.staged_recording = false;
param.tg.altitude_guard = 7500*12*2.54/100;
param.wfs(1).atten = 20;
param.wfs(2).atten = 20;
param.wfs(3).atten = 20;
param.tg.Haltitude = 12500*12*2.54/100;
param.tg.Hice_thick = 0;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_3500mthick_DECONVOLUTION.xml');
write_cresis_xml(param);

%% Survey Mode High Altitude
% <3200 m thick ice, 5000-20000 ft AGL
freq_idx = 1;
param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 7500*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 12500*12*2.54/100;
param.tg.Hice_thick = 3200;
param.fn = fullfile(base_dir,'survey_mode_10us_2wf_3200mthick_high_altitude.xml');
param.prf = 9000;
param.presums = [5 25];
param.wfs(1).atten = 20;
param.wfs(2).atten = 0;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 10e-6;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.delay = final_DDS_time{cal_settings(freq_idx)};
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
[param.wfs(:).tx_mask] = deal(final_tx_mask);
write_cresis_xml(param);

%% Survey Mode Thin Ice
% <2500 m thick ice, 1250 +/- 750 ft AGL
freq_idx = 1;
param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,sprintf('survey_mode_3us_2wf_%.0fmthick_thin_ice.xml',param.tg.Hice_thick));
param.prf = 12000;
param.presums = [5 29];
param.wfs(1).atten = 33;
param.wfs(2).atten = 0;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.delay = final_DDS_time{cal_settings(freq_idx)};
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
[param.wfs(:).tx_mask] = deal(final_tx_mask);
write_cresis_xml(param);

%% Image Mode Thick Ice
% All Ice <3500 m, 500 +/- 250 m AGL
param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = [1 1 2 2 3 3];
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = 3500;
param.fn = fullfile(base_dir,'image_mode_10us_6wf_3500mthick.xml');
param.prf = 12000;
param.presums = [3 3 3 3 13 13];
param.wfs(1).atten = 26;
param.wfs(2).atten = 26;
param.wfs(3).atten = 0;
param.wfs(4).atten = 0;
param.wfs(5).atten = 0;
param.wfs(6).atten = 0;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)} .* [hanning(7).', 0] ./ Hwindow_orig;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)} .* [0 0 hanning(5).', 0] ./ Hwindow_orig; % Antenna 2 bad
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 1e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(4).Tpd = 3e-6;
param.wfs(5).Tpd = 10e-6;
param.wfs(6).Tpd = 10e-6;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.tg.look_angle_deg = [-15 15 -15 15 -15 15];
for wf = 1:length(param.tg.look_angle_deg)
  nadir_vec = [0 0 1];
  beam_vec = [0 sind(param.tg.look_angle_deg(wf)) cosd(param.tg.look_angle_deg(wf))];
  nadir_delay = -nadir_vec * phase_centers;
  beam_delay = -beam_vec * phase_centers;
  nadir_delay = nadir_delay - nadir_delay(4);
  beam_delay = beam_delay - beam_delay(4);
  nadir_delay = nadir_delay / c;
  beam_delay = beam_delay / c - nadir_delay; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay) 0];
  param.wfs(wf).delay = (final_DDS_time{cal_settings(freq_idx)}/1e9 - beam_delay)*1e9;
end
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
[param.wfs(:).tx_mask] = deal(final_tx_mask);
write_cresis_xml(param);

%% Image Mode Thin Ice
% All Ice <2500 m, 500 +/- 250 m AGL
param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 750*12*2.54/100;
param.tg.staged_recording = [1 1 2 2];
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = 2500;
param.fn = fullfile(base_dir,'image_mode_10us_4wf_2000mthin.xml');
param.prf = 12000;
param.presums = [3 3 15 15];
param.wfs(1).atten = 26;
param.wfs(2).atten = 26;
param.wfs(3).atten = 0;
param.wfs(4).atten = 0;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)} .* [hanning(7).', 0] ./ Hwindow_orig;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)} .* [0 0 hanning(5).', 0] ./ Hwindow_orig; % Antenna 2 bad
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 1e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(4).Tpd = 3e-6;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.tg.look_angle_deg = [-15 15 -15 15];
for wf = 1:length(param.tg.look_angle_deg)
  nadir_vec = [0 0 1];
  beam_vec = [0 sind(param.tg.look_angle_deg(wf)) cosd(param.tg.look_angle_deg(wf))];
  nadir_delay = -nadir_vec * phase_centers;
  beam_delay = -beam_vec * phase_centers;
  nadir_delay = nadir_delay - nadir_delay(4);
  beam_delay = beam_delay - beam_delay(4);
  nadir_delay = nadir_delay / c;
  beam_delay = beam_delay / c - nadir_delay; % Only include delay relative to nadir
  beam_delay = [beam_delay - mean(beam_delay) 0];
  param.wfs(wf).delay = (final_DDS_time{cal_settings(freq_idx)}/1e9 - beam_delay)*1e9;
end
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
[param.wfs(:).tx_mask] = deal(final_tx_mask);
write_cresis_xml(param);


%% Ping pong Thin Ice
% All Ice <2500 m, 500 +/- 250 m AGL
param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 2500*12*2.54/100;
param.tg.staged_recording = [0 0];
param.tg.Haltitude = 10000*12*2.54/100;
param.tg.Hice_thick = 2000;
param.fn = fullfile(base_dir,'pingpong_mode_10us_4wf_2000mthin.xml');
param.prf = 10000;
param.presums = [15 15];
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
DDS_amp = [40000 0 0 0 0 40000 40000 0];
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.tg.look_angle_deg = [30 30];
param.delay = final_DDS_time{cal_settings(freq_idx)};
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
[param.wfs(1).tx_mask] = deal([1 0 0 1 1 1 1 1]);
[param.wfs(2).tx_mask] = deal([1 1 1 1 1 1 1 0]);
write_cresis_xml(param);

%% Transmit waveform test mode
% These files are for measuring the transmit waveforms with a DSO
% Attenuators should be connected to each power amp output
% MAX_DATA_RATE DOES NOT MATTER IN THIS MODE, NOT USED TO CAPTURE
% freq_idx = 1;
% param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
% param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 500; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
% param.max_duty_cycle = 0.12;
% param.create_IQ = false;
% param.tg.staged_recording = false;
% param.tg.altitude_guard = 0*12*2.54/100;
% param.tg.Haltitude = 0*12*2.54/100;
% param.tg.Hice_thick = 2e-6 * 3e8/2/sqrt(3.15); % We don't capture in this mode, so does not matter
% param.prf = 12000;
% param.presums = [1];
% param.wfs(1).atten = 0;
% DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
% param.tx_weights = DDS_amp;
% param.phase = final_DDS_phase{cal_settings(freq_idx)};
% param.delay = final_DDS_time{cal_settings(freq_idx)};
% [param.wfs(:).tx_mask] = deal(final_tx_mask);
% % 1e-6 chirp, 0.1 tukey
% param.tukey = 0.2;
% param.Tpd_base_duration = 1e-6;
% param.wfs(1).Tpd = 1e-6;
% param.f0 = f0_list(freq_idx);
% param.f1 = f1_list(freq_idx);
% param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
% write_cresis_xml(param);
% % 3e-6 chirp, 0.1 tukey
% param.tukey = 0.2;
% param.Tpd_base_duration = 1e-6;
% param.wfs(1).Tpd = 3e-6;
% param.f0 = f0_list(freq_idx);
% param.f1 = f1_list(freq_idx);
% param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
% write_cresis_xml(param);
% % 10e-6 chirp, 0.1 tukey
% param.tukey = 0.2;
% param.Tpd_base_duration = 1e-6;
% param.wfs(1).Tpd = 10e-6;
% param.f0 = f0_list(freq_idx);
% param.f1 = f1_list(freq_idx);
% param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
% write_cresis_xml(param);
% % 30 ns CW, 0 tukey
% param.tukey = 0.3;
% param.Tpd_base_duration = 30e-9;
% param.wfs(1).Tpd = 30e-9;
% param.f0 = 195e6;
% param.f1 = 195e6;
% param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
% write_cresis_xml(param);

%% Equalization High Altitude (Using Ocean)
% 5000-20000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
freq_idx = 1;
for Tpd = [10e-6]
  param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);  
  param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.altitude_guard = 7500*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 12500*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('equalization_%.0fus.xml', Tpd*1e6));
  param.prf = 8500;
  param.presums = [15 15 15 15 15 15 15 15];
  [param.wfs(1:8).atten] = deal(15);
  param.wfs(8).atten = deal(30);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.tukey = 0.2;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  for wf=1:7
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(8).tx_mask = zeros(1,8);
  write_cresis_xml(param);
end

%% Equalization Low Altitude
% 1250 +/- 750 ft AGL
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
freq_idx = 1;
for Tpd = [1e-6]
  param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);  
  param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.altitude_guard = 6000*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 6000*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('equalization_%.0fus_%.0fft.xml', Tpd*1e6, param.tg.Haltitude*100/2.54/12));
  param.prf = 12000;
  param.presums = [15 15 15 15 15 15 15 15];
  [param.wfs(1:8).atten] = deal(30);
  param.wfs(8).atten = deal(45);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.tukey = 0.2;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  for wf=1:7
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(8).tx_mask = zeros(1,8);
  write_cresis_xml(param);
end

%% Equalization Lab 
% Use these settings in the lab to demonstrate transmit equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
freq_idx = 1;
for Tpd = 10e-6
  param = struct('radar_name','mcords3','num_chan',15,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1e9/9,'fs_sync',1e9/18,'fs_dds',1e9,'TTL_clock',1e9/18,'TTL_mode',[3e-6 290e-9 -1060e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);  
  param.max_tx = [40000 40000 40000 40000 40000 40000 40000 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 12.18e-6; param.final_tx_mask = final_tx_mask;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.altitude_guard = 5000*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 5000*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(calval_dir,sprintf('equalization_%.0fus_LOOPBACK.xml', Tpd*1e6));
  param.prf = 8500;
  param.presums = [15 15 15 15 15 15 15 15];
  [param.wfs(1:8).atten] = deal(5);
  param.wfs(8).atten = deal(20);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.tukey = 0.2;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  for wf=1:7
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(8).tx_mask = zeros(1,8);
  write_cresis_xml(param);
end
