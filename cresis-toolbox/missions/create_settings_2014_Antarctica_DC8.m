% script create_settings_2014_Antarctica_DC8
%
% Creates NI radar depth sounder settings

% Define waveforms
base_dir = 'D:\waveforms';
if 0
  % Initial conditions
  final_DDS_phase = [0 0 0 0 0 0  0 0];
  final_DDS_phase_no_time = [0 0 0 0 0 0  0 0];
  final_DDS_amp = [15000 15000 17000 15920 15000 14150 0 0];
  final_DDS_time =  [0 0 0 0 0 0  0 0];
else
  % After transmit calibration during Oct 7 test flight
  final_DDS_phase = [-77.1	-144.3	-25.2	0.0	-125.1	59.4	0.0	0.0];
  final_DDS_amp = [9312	13066	16610	11214	15000	10915	0	0];
  final_DDS_time =  [-1.06	2.37	-1.30	0.00	-0.77	-4.58	0.00	0.00];
  final_DDS_phase_no_time = [-1.7	45.7	52.5	0.0	-85.1	6.9	0.0	0.0];
end

f0 = 165e6;
f1 = 215e6;

Hwindow_orig = [1 1 1 1 1 1 1 1]; % Desired window created during transmit calibration

physical_constants;
param = [];
param.season_name = '2014_Antarctica_DC8';
param.radar_name = 'rds';
param.gps_source = 'atm-final20140301';
clear phase_centers;
for tx_chan = 1:6
  tx_weights = zeros(1,6);
  tx_weights(tx_chan) = 1;
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, tx_chan);
end

%% Survey Mode + noise mode + loopback mode
% <4500 m thick ice, 0-2500 ft AGL
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);
param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 8.9e-6;
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
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
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
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.tg.altitude_guard = 0*12*2.54/100;
param.tg.Haltitude = 0*12*2.54/100;
param.tg.Hice_thick = 20e-6 * 3e8/2/sqrt(3.15); % Long enough for delay line
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_LOOPBACK.xml');
write_cresis_xml(param);
% Deconvolution Mode (for Sea Ice)
param.wfs(1).atten = 50;
param.wfs(2).atten = 50;
param.wfs(3).atten = 50;
param.tg.Hice_thick = 35e-6 * 3e8/2/sqrt(3.15); % Record length for up to 12000 ft in air
param.fn = fullfile(base_dir,'survey_mode_10us_3wf_4500mthick_DECONVOLUTION.xml');
write_cresis_xml(param);

%% Survey Mode High Altitude
% <3200 m, 25000 ft AGL
% <3200 m thick ice, 20000-35000 ft AGL
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);
param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 8.9e-6;
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
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
write_cresis_xml(param);

%% Transmit waveform test mode
% These files are for measuring the transmit waveforms with a DSO
% Attenuators should be connected to each power amp output
% MAX_DATA_RATE DOES NOT MATTER IN THIS MODE, NOT USED TO CAPTURE
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);
param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 500; param.flight_hours = 7; param.sys_delay = 8.9e-6;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = false;
param.tg.altitude_guard = 0*12*2.54/100;
param.tg.Haltitude = 0*12*2.54/100;
param.tg.Hice_thick = 2e-6 * 3e8/2/sqrt(3.15); % We don't capture in this mode, so does not matter
param.prf = 12000;
param.presums = [1];
param.wfs(1).atten = 0;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
[param.wfs(:).tx_mask] = deal([0 0 0 0 0 0 0 0]);
% 1e-6 chirp, 0.1 tukey
param.tukey = 0.1;
param.Tpd_base_duration = 1e-6;
param.wfs(1).Tpd = 1e-6;
param.f0 = f0;
param.f1 = f1;
param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
write_cresis_xml(param);
% 3e-6 chirp, 0.1 tukey
param.tukey = 0.1;
param.Tpd_base_duration = 1e-6;
param.wfs(1).Tpd = 3e-6;
param.f0 = f0;
param.f1 = f1;
param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
write_cresis_xml(param);
% 10e-6 chirp, 0.1 tukey
param.tukey = 0.1;
param.Tpd_base_duration = 1e-6;
param.wfs(1).Tpd = 10e-6;
param.f0 = f0;
param.f1 = f1;
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
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);
param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 8.9e-6;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 1250*12*2.54/100;
param.tg.staged_recording = false;
param.wfs(1).tg.Haltitude = 1250*12*2.54/100;
param.wfs(2).tg.Haltitude = 1250*12*2.54/100 + 10e-6*3e8/2;
param.wfs(3).tg.Haltitude = 1250*12*2.54/100 + 10e-6*3e8/2;
param.wfs(1).tg.Hice_thick = 1000;
param.wfs(2).tg.Hice_thick = 3500;
param.wfs(3).tg.Hice_thick = 3500;
param.fn = fullfile(base_dir,'pingpong_mode_10us_3500mthick_3wf.xml');
param.prf = 12000;
param.presums = [3 15 15];
param.wfs(1).atten = 30;
param.wfs(2).atten = 15;
param.wfs(3).atten = 15;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 10e-6;
param.wfs(3).Tpd = 10e-6;
param.phase = final_DDS_phase;
param.delay = final_DDS_time;
param.wfs(1).f0 = f0;
param.wfs(1).f1 = f1;
param.wfs(2).f0 = f0;
param.wfs(2).f1 = f1;
param.wfs(3).f0 = f1;
param.wfs(3).f1 = f0;
[param.wfs(1).tx_mask] = [0 0 0 0 0 0 0 0];
[param.wfs(2).tx_mask] = [1 1 0 1 1 0 1 1];
[param.wfs(3).tx_mask] = [1 1 1 1 0 1 1 0];
write_cresis_xml(param);

%% Image Mode (Desert Mode)
% Ping pong transmit mode, narrow pulse, restricted bandwidth for FCC, 0-2500 ft AGL
param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);
param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 8.9e-6; param.Tpd_base_duration = 30e-9;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 1250*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 1250*12*2.54/100;
param.tg.Hice_thick = 1500;
param.fn = fullfile(base_dir,'desert_mode.xml');
param.prf = 12000;
param.presums = [17 17];
param.wfs(1).atten = 25;
param.wfs(2).atten = 25;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.3;
param.wfs(1).Tpd = 30e-9;
param.wfs(2).Tpd = 30e-9;
param.phase = final_DDS_phase_no_time;
param.delay = zeros(size(final_DDS_time));
param.wfs(1).f0 = 195e6;
param.wfs(1).f1 = 195e6;
param.wfs(2).f0 = 195e6;
param.wfs(2).f1 = 195e6;
[param.wfs(1).tx_mask] = [1 1 0 1 1 0 1 1];
[param.wfs(2).tx_mask] = [1 1 1 1 0 1 1 0];
write_cresis_xml(param);

%% Equalization High Altitude (Using Ocean)
% 10000-45000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);  
  param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 8.9e-6;
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
  param.tx_weights = final_DDS_amp;
  param.tukey = 0.1;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase;
  param.delay = final_DDS_time;
  param.f0 = f0;
  param.f1 = f1;
  for wf=1:6
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(7).tx_mask = zeros(1,8);
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
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);  
  param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 8.9e-6;
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
  param.tx_weights = final_DDS_amp;
  param.tukey = 0.1;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase;
  param.delay = final_DDS_time;
  param.f0 = f0;
  param.f1 = f1;
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
for Tpd = 10e-6
  param = struct('radar_name','mcords3','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',150e6,'fs_sync',56.25e6,'fs_dds',900e6,'TTL_mode',[1.3e-6 0.3e-6 -880e-9]);  
  param.max_tx = [15000 15000 17000 15920 15000 14150 0 0]; param.max_data_rate = 135; param.flight_hours = 7; param.sys_delay = 8.9e-6;
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
  param.tx_weights = final_DDS_amp;
  param.tukey = 0.1;
  param.Tpd = Tpd;
  param.phase = final_DDS_phase;
  param.delay = final_DDS_time;
  param.f0 = f0;
  param.f1 = f1;
  for wf=1:6
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  param.wfs(7).tx_mask = zeros(1,8);
  write_cresis_xml(param);
end
