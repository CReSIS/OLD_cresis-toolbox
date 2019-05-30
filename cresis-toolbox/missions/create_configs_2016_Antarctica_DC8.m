% script create_configs_2016_Antarctica_DC8
%
% Creates NI radar depth sounder settings

physical_constants; % c = speed of light

% Define waveforms
if ispc
  base_dir = 'e:\waveforms\';
else
  base_dir = '/scratch/waveforms/';
end

f0_list = [165e6];
f1_list = [215e6];
DDC_select_list = [1]; % Which DDC mode to use
cal_settings = [1];
prf = 10000;

% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
velocity = [130]; % m/s
presums = round(c./abs(f1_list+f0_list)/2 ./ velocity * prf / 2)*2;

final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 0
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0 0 0 0 0 0 0 0];
    final_DDS_phase_no_time{idx} = [0 0 0 0 0 0 0 0]; % not used usually
    final_DDS_amp{idx} = [57750 57750 65450 61292 61600 54478 0 0];
    final_DDS_amp{idx} = [57750 0 65450 61292 61600 54478 0 0];
    final_DDS_time{idx} =  [0 0 0 0 0 0 0 0];
  end
  final_tx_mask = [1 1 0 0 0 0 0 0];
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % After transmit calibration during Oct 4, 2016 test flight (3 us)
  % final_DDS_phase{end+1} = [17.9	0.0	98.8	-57.9	-92.9	63.2	0.0	0.0];
  % final_DDS_amp{end+1} = [37573	44712	43831	37682	61600	36778	0	0];
  % final_DDS_time{end+1} =  [-4.78	0.00	-3.71	-4.80	-3.78	-8.57	0.00	0.00];
  % final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0  0 0];
  
  % After transmit calibration during Oct 4, 2016 test flight (10 us) <-- Use 10 us!
  final_DDS_phase{end+1} = [-10.8	0.0	94.7	-53.1	-85.5	65.7	0.0	0.0];
  final_DDS_amp{end+1} = [41142	44690	44087	43212	61600	38806	0	0];
  final_DDS_time{end+1} =  [-4.78	0.00	-3.71	-4.80	-3.78	-8.57	0.00	0.00];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0  0 0];
  final_tx_mask = [1 1 0 0 0 0 0 0];
end

Hwindow_orig = [1 1 1 1 1 1 0 0]; % Desired window created during transmit calibration

physical_constants;
param = [];
param.season_name = '2016_Antarctica_DC8';
param.radar_name = 'rds';
param.gps_source = 'atm-final20160301';
clear phase_centers;
for tx_chan = 1:6
  tx_weights = zeros(1,6);
  tx_weights(tx_chan) = 1;
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, tx_chan);
end

%% Survey Mode + loopback, noise, and deconv modes
% <4000 m thick ice, 1250 +/- 1250 ft AGL
ice_thickness = [4500 4500];
freq_idx = 1;
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 8.9e-6;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 1250*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = 4500;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick.xml');
param.prf = 12000;
param.presums = [3 3 25];
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
% Noise Mode
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_NOISE.xml');
[param.wfs(:).tx_mask] = deal([1 1 1 1 1 1 1 1]);
write_cresis_xml(param);
% Loopback Mode
param.tg.staged_recording = false;
[param.wfs(:).tx_mask] = deal(final_tx_mask);
param.tg.altitude_guard = 0*12*2.54/100;
param.tg.Haltitude = 0*12*2.54/100;
param.tg.Hice_thick = 20e-6 * 3e8/2/sqrt(3.15); % Long enough for delay line
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_LOOPBACK.xml');
write_cresis_xml(param);
% Deconvolution Mode (for Sea Ice)
param.prf = 8500;
param.tg.staged_recording = false;
param.tg.altitude_guard = 17500*12*2.54/100;
param.wfs(1).atten = 20;
param.wfs(2).atten = 20;
param.wfs(3).atten = 20;
param.tg.Haltitude = 27500*12*2.54/100;
param.tg.Hice_thick = 0;
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_DECONVOLUTION.xml');
write_cresis_xml(param);

%% Survey Mode High Altitude
% <3200 m thick ice, 20000-35000 ft AGL
freq_idx = 1;
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 8.9e-6;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 7500*12*2.54/100;
param.tg.staged_recording = true;
param.tg.Haltitude = 27500*12*2.54/100;
param.tg.Hice_thick = 3200;
param.fn = fullfile(base_dir,'survey_mode_10us_2wf_3200mthick_high_altitude.xml');
param.prf = 7500;
param.presums = [3 15];
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

%% Transmit waveform test mode
% These files are for measuring the transmit waveforms with a DSO
% Attenuators should be connected to each power amp output
% MAX_DATA_RATE DOES NOT MATTER IN THIS MODE, NOT USED TO CAPTURE
freq_idx = 1;
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 500; param.flight_hours = 7; param.sys_delay = 8.9e-6;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = false;
param.tg.altitude_guard = 0*12*2.54/100;
param.tg.Haltitude = 0*12*2.54/100;
param.tg.Hice_thick = 2e-6 * 3e8/2/sqrt(3.15); % We don't capture in this mode, so does not matter
param.prf = 12000;
param.presums = [1];
param.wfs(1).atten = 0;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
param.tx_weights = DDS_amp;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.delay = final_DDS_time{cal_settings(freq_idx)};
[param.wfs(:).tx_mask] = deal(final_tx_mask);
% 1e-6 chirp, 0.1 tukey
param.tukey = 0.2;
param.Tpd_base_duration = 1e-6;
param.wfs(1).Tpd = 1e-6;
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
write_cresis_xml(param);
% 3e-6 chirp, 0.1 tukey
param.tukey = 0.2;
param.Tpd_base_duration = 1e-6;
param.wfs(1).Tpd = 3e-6;
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
write_cresis_xml(param);
% 10e-6 chirp, 0.1 tukey
param.tukey = 0.2;
param.Tpd_base_duration = 1e-6;
param.wfs(1).Tpd = 10e-6;
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
write_cresis_xml(param);
% 30 ns CW, 0 tukey
param.tukey = 0.3;
param.Tpd_base_duration = 30e-9;
param.wfs(1).Tpd = 30e-9;
param.f0 = 195e6;
param.f1 = 195e6;
param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
write_cresis_xml(param);

%% Image Mode
% Ping pong transmit mode, 0-2500 ft AGL
freq_idx = 1;
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);
param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 8.9e-6;
param.max_duty_cycle = 0.12;
param.create_IQ = false;

param.tg.staged_recording = [1 2 2 3 3];
param.tg.altitude_guard = 1250*12*2.54/100;
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = 3500;

param.fn = fullfile(base_dir,'pingpong_mode_10us_3500mthick_3wf.xml');
param.prf = 12000;
param.presums = [3 3 3 13 13];
param.wfs(1).atten = 23;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
param.wfs(4).atten = 0;
param.wfs(5).atten = 0;
DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
param.tx_weights = DDS_amp;
param.tukey = 0.2;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(4).Tpd = 10e-6;
param.wfs(5).Tpd = 10e-6;
param.phase = final_DDS_phase{cal_settings(freq_idx)};
param.delay = final_DDS_time{cal_settings(freq_idx)};
param.wfs(1).f0 = f0_list(freq_idx);
param.wfs(1).f1 = f1_list(freq_idx);
param.wfs(2).f0 = f0_list(freq_idx);
param.wfs(2).f1 = f1_list(freq_idx);
param.wfs(3).f0 = f0_list(freq_idx);
param.wfs(3).f1 = f1_list(freq_idx);
param.wfs(4).f0 = f0_list(freq_idx);
param.wfs(4).f1 = f1_list(freq_idx);
param.wfs(5).f0 = f0_list(freq_idx);
param.wfs(5).f1 = f1_list(freq_idx);
[param.wfs(1).tx_mask] = [1 1 0 1 1 0 1 1];
[param.wfs(2).tx_mask] = [1 1 0 1 1 0 1 1];
[param.wfs(3).tx_mask] = [1 1 1 1 0 1 1 0];
[param.wfs(4).tx_mask] = [1 1 0 1 1 0 1 1];
[param.wfs(5).tx_mask] = [1 1 1 1 0 1 1 0];
write_cresis_xml(param);

%% Equalization High Altitude (Using Ocean)
% 10000-45000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
freq_idx = 1;
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);  
  param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 8.9e-6;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.altitude_guard = 17500*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 27500*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('equalization_%.0fus.xml', Tpd*1e6));
  param.prf = 8500;
  param.presums = [15 15 15 15 15 15 15];
  [param.wfs(1:7).atten] = deal(15);
  param.wfs(7).atten = deal(30);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.tukey = 0.2;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  for wf=1:6
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(7).tx_mask = zeros(1,8);
  write_cresis_xml(param);
end

%% Ping Pong Paired Transmitter High Altitude (Using Ocean)
% 10000-45000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
freq_idx = 1;
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);  
  param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 8.9e-6;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.altitude_guard = 17500*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 27500*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('equalization_pairedtx_%.0fus.xml', Tpd*1e6));
  param.prf = 8500;
  param.presums = [15 15 15];
  [param.wfs(1:3).atten] = deal(21);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.tukey = 0.2;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.wfs(1).tx_mask = [1 1 1 1 0 1 1 0];
  param.wfs(2).tx_mask = [1 1 1 0 1 1 0 1];
  param.wfs(3).tx_mask = [1 1 0 1 1 0 1 1];
  write_cresis_xml(param);
end

%% Equalization High Altitude (Using Ocean) Without Combined
% 10000-45000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization
% Creates one waveform for each of N DDS-transmitters
%
% NOTE: THESE SETTING DO NOT HAVE COMBINED WAVEFORM... THESE ARE USEFUL
% WHEN ACTIVE S11 FROM BADLY MATCHED CHANNELS COULD DAMAGE PA's.
freq_idx = 1;
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);  
  param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 8.9e-6;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.altitude_guard = 17500*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 27500*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('equalization_%.0fus_nocombinedwf.xml', Tpd*1e6));
  param.prf = 8500;
  param.presums = [15 15 15 15 15 15];
  [param.wfs(1:6).atten] = deal(5);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.tukey = 0.2;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  for wf=1:6
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  write_cresis_xml(param);
end

%% Equalization Lab 
% Use these settings in the lab to demonstrate transmit equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
freq_idx = 1;
for Tpd = 10e-6
  param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',900e6/16,'fs_dds',900e6,'TTL_clock',150e6/2,'TTL_mode',[3e-6 0.35e-6 -880e-9],'xml_version',2.0,'DDC_freq',0,'DDC_select',0);  
  param.max_tx = [57750 57750 65450 61292 61600 54478 0 0]; param.max_data_rate = 150; param.flight_hours = 7; param.sys_delay = 8.9e-6;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.altitude_guard = 5000*12*2.54/100;
  param.tg.staged_recording = false;
  param.tg.Haltitude = 5000*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('equalization_%.0fus_LOOPBACK.xml', Tpd*1e6));
  param.prf = 8500;
  param.presums = [15 15 15 15 15 15 15];
  [param.wfs(1:7).atten] = deal(5);
  param.wfs(7).atten = deal(20);
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.tukey = 0.2;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  for wf=1:6
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(7).tx_mask = zeros(1,8);
  write_cresis_xml(param);
end
