% script create_settings_2015_Greenland_LC130
%
% Creates NI radar depth sounder settings

% Define waveforms
base_dir = '/scratch/waveforms';
if 0
  % Initial conditions
  final_DDS_phase = [0 0 0 0 0 0 0 0];
  final_DDS_phase_no_time = [0 0 0 0 0 0 0 0];
  final_DDS_amp = [30000 30000 0 0 0 0 0 0];
  final_DDS_time =  [0 0 0 0 0 0 0 0];
else
  % After transmit calibration during Mar ?, 2015 test flight
  final_DDS_phase = [0 -85.9 0 0 0 0 0 0];
  final_DDS_phase_no_time = [0 0 0 0 0 0 0 0];
  final_DDS_amp = [30000 26643 0 0 0 0 0 0];
  final_DDS_time =  [0 3.26 0 0 0 0 0 0];
end

f0 = 180e6;
f1 = 450e6;

Hwindow_orig = [1 1 0 0 0 0 0 0]; % Desired window created during transmit calibration

physical_constants;
param = [];
param.season_name = '2015_Greenland_LC130';1
param.radar_name = 'rds';
param.gps_source = 'atm-final';
clear phase_centers;
for tx_chan = 1:2
  % Just enable the current antenna to determine its phase center
  tx_weights = zeros(1,2);
  tx_weights(tx_chan) = 1;
  % Determine phase center for the antenna
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, tx_chan);
end

%% Survey Mode + loopback mode
% <3500 m thick ice, 0-2500 ft AGL
param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 315e6; param.DDC_select = 1;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = [1 2 3 3];
param.tg.Haltitude = 1500*12*2.54/100;
param.tg.Hice_thick = 3500;
param.fn = fullfile(base_dir,'survey_mode_10us_3500mthick.xml');
param.prf = 12000;
param.presums = [3 3 17];
param.wfs(1).atten = 33;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.wfs(1).phase = final_DDS_phase;
param.wfs(2).phase = final_DDS_phase;
param.wfs(3).phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(1:3).tx_mask] = deal([1 1 1 1 1 1 0 0]);
write_cresis_xml(param);
% Default Mode
param.fn = fullfile(base_dir,'default.xml');
write_cresis_xml(param);
% Loopback Mode
param.DDC_select = 0;
param.tg.staged_recording = false;
param.tg.altitude_guard = 0*12*2.54/100;
param.tg.Haltitude = 0*12*2.54/100;
param.tg.Hice_thick = 12e-6 * 3e8/2/sqrt(3.15); % Long enough for 10 us delay line
param.fn = fullfile(base_dir,'survey_mode_10us_3500mthick_LOOPBACK.xml');
write_cresis_xml(param);
param.DDC_select = 1;
% Loopback Mode
param.DDC_select = 0;
param.f0 = 180e6;
param.f1 = 230e6;
param.tg.staged_recording = false;
param.tg.altitude_guard = 0*12*2.54/100;
param.tg.Haltitude = 0*12*2.54/100;
param.tg.Hice_thick = 12e-6 * 3e8/2/sqrt(3.15); % Long enough for 10 us delay line
param.fn = fullfile(base_dir,'glacier_mode_10us_3500mthick_LOOPBACK.xml');
write_cresis_xml(param);
param.f0 = f0;
param.f1 = f1;
param.DDC_select = 1;
% Deconvolution Mode (for over Sea Ice)
param.wfs(1).atten = 30;
param.wfs(2).atten = 30;
param.wfs(3).atten = 30;
param.tg.staged_recording = false;
param.tg.altitude_guard = 7500*12*2.54/100;
param.tg.Haltitude = 12500*12*2.54/100;
param.tg.Hice_thick = 0 * 12*2.54/100/1.78;
param.fn = fullfile(base_dir,'survey_mode_10us_3500mthick_DECONVOLUTION.xml');
write_cresis_xml(param);

%% Noise Mode
param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 0e6; param.DDC_select = 0;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 2000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 2000*12*2.54/100;
param.tg.Hice_thick = 0;
param.fn = fullfile(base_dir,'noise_mode.xml');
param.prf = 12000;
param.presums = [9];
param.wfs(1).atten = 0;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(1).phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(1:1).tx_mask] = deal([1 1 1 1 1 1 1 1]);
write_cresis_xml(param);


% %% Survey Mode high altitude
% % <3000 m thick ice, 11000-17000 ft AGL
% param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
% param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 315e6; param.DDC_select = 1;
% param.max_duty_cycle = 0.12;
% param.create_IQ = false;
% param.tg.altitude_guard = 3000*12*2.54/100;
% param.tg.staged_recording = false;
% param.tg.Haltitude = 14000*12*2.54/100;
% param.tg.Hice_thick = 3000;
% param.fn = fullfile(base_dir,'survey_mode_10us_highaltitude.xml');
% param.prf = 10000;
% param.presums = [13 13];
% param.wfs(1).atten = 0;
% param.wfs(2).atten = 0;
% DDS_amp = final_DDS_amp;
% param.tx_weights = DDS_amp;
% param.tukey = 0.1;
% param.wfs(1).Tpd = 10e-6;
% param.wfs(2).Tpd = 10e-6;
% param.wfs(1).phase = final_DDS_phase;
% param.wfs(2).phase = final_DDS_phase;
% param.delay = final_DDS_time;
% param.f0 = f0;
% param.f1 = f1;
% [param.wfs(1).tx_mask] = deal([1 1 1 1 1 1 0 0]);
% [param.wfs(2).tx_mask] = deal([1 1 1 1 1 1 0 1]);
% write_cresis_xml(param);

%% Glacier Mode
% Narrow bandwidth <3500 m thick ice, 0-2500 ft AGL
param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 205e6; param.DDC_select = 2;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = [1 2 3 3];
param.tg.Haltitude = 1500*12*2.54/100;
param.tg.Hice_thick = 3500;
param.fn = fullfile(base_dir,'glacier_mode_10us_3500mthick.xml');
param.prf = 12000;
param.presums = [3 3 11 11];
param.wfs(1).atten = 33;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
param.wfs(4).atten = 0;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 10e-6;
param.wfs(4).Tpd = 10e-6;
param.wfs(1).phase = final_DDS_phase;
param.wfs(2).phase = final_DDS_phase;
param.wfs(3).phase = final_DDS_phase;
param.wfs(4).phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = 180e6;
param.f1 = 230e6;
[param.wfs(1:3).tx_mask] = deal([1 1 1 1 1 1 0 0]);
[param.wfs(4).tx_mask] = deal([1 1 1 1 1 1 0 1]);
write_cresis_xml(param);
% Deconvolution Mode (for over Sea Ice)
param.wfs(1).atten = 30;
param.wfs(2).atten = 30;
param.wfs(3).atten = 30;
param.wfs(4).atten = 30;
param.tg.staged_recording = false;
param.tg.altitude_guard = 7500*12*2.54/100;
param.tg.Haltitude = 12500*12*2.54/100;
param.tg.Hice_thick = 0 * 12*2.54/100/1.78;
param.fn = fullfile(base_dir,'glacier_mode_10us_3500mthick_DECONVOLUTION.xml');
write_cresis_xml(param);

%% Shallow Ice Mode
% <1500 m thick ice, 0-2500 ft AGL
param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 315e6; param.DDC_select = 1;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = [1 2 2];
param.tg.Haltitude = 1500*12*2.54/100;
param.tg.Hice_thick = 2000;
param.fn = fullfile(base_dir,'shallow_mode_3us_1500mthick.xml');
param.prf = 12000;
param.presums = [3 11 11];
param.wfs(1).atten = 33;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase;
param.wfs(2).phase = final_DDS_phase;
param.wfs(3).phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
[param.wfs(1:2).tx_mask] = deal([1 1 1 1 1 1 0 0]);
[param.wfs(3).tx_mask] = deal([1 1 1 1 1 1 0 1]);
write_cresis_xml(param);

%% Sea ice mode
% <100 m thick ice, 0-2500 ft AGL
param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 14.4e-6;  param.DDC_freq = 0e6; param.DDC_select = 0;
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.staged_recording = false;
param.tg.Haltitude = 1500*12*2.54/100;
param.tg.Hice_thick = 100;
param.fn = fullfile(base_dir,'seaice_mode_1us.xml');
param.prf = 12000;
param.presums = [11 11];
param.wfs(1).atten = 33;
param.wfs(2).atten = 33;
DDS_amp = final_DDS_amp;
param.tx_weights = DDS_amp/2;
param.tukey = 0.1;
param.wfs(1).Tpd = 1e-6;
param.wfs(2).Tpd = 1e-6;
param.wfs(1).phase = final_DDS_phase;
param.wfs(2).phase = final_DDS_phase;
param.delay = final_DDS_time;
param.f0 = f0;
param.f1 = f1;
param.wfs(1).tx_mask = [1 1 1 1 1 1 1 0];
param.wfs(2).tx_mask = [1 1 1 1 1 1 0 1];
write_cresis_xml(param);

% %% Transmit waveform test mode
% % These files are for measuring the transmit waveforms with a DSO
% % Attenuators should be connected to each power amp output
% % MAX_DATA_RATE DOES NOT MATTER IN THIS MODE, NOT USED TO CAPTURE
% param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
% param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = inf; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 315e6; param.DDC_select = 1;
% param.max_duty_cycle = 0.12;
% param.create_IQ = false;
% param.tg.staged_recording = false;
% param.tg.altitude_guard = 0*12*2.54/100;
% param.tg.Haltitude = 0*12*2.54/100;
% param.tg.Hice_thick = 2e-6 * 3e8/2/sqrt(3.15); % We don't capture in this mode, so does not matter
% param.prf = 12000;
% param.presums = [1];
% param.wfs(1).atten = 0;
% DDS_amp = final_DDS_amp;
% param.tx_weights = DDS_amp;
% param.phase = final_DDS_phase;
% param.delay = final_DDS_time;
% [param.wfs(:).tx_mask] = deal([1 1 1 1 1 1 0 0]);
% % 1e-6 chirp, 0.1 tukey
% param.tukey = 0.1;
% param.Tpd_base_duration = 1e-6;
% param.wfs(1).Tpd = 1e-6;
% param.f0 = f0;
% param.f1 = f1;
% param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
% write_cresis_xml(param);
% % 3e-6 chirp, 0.1 tukey
% param.tukey = 0.1;
% param.Tpd_base_duration = 1e-6;
% param.wfs(1).Tpd = 3e-6;
% param.f0 = f0;
% param.f1 = f1;
% param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
% write_cresis_xml(param);
% % 10e-6 chirp, 0.1 tukey
% param.tukey = 0.1;
% param.Tpd_base_duration = 1e-6;
% param.wfs(1).Tpd = 10e-6;
% param.f0 = f0;
% param.f1 = f1;
% param.fn = fullfile(base_dir,sprintf('transmit_wf_capture_DSO_%dns.xml', round(param.wfs(1).Tpd*1e9)));
% write_cresis_xml(param);

%% Equalization High Altitude (Using Ocean)
% 10000-45000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
  param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 315e6; param.DDC_select = 1;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = false;
  if Tpd == 1e-6
    param.tg.altitude_guard = 2000*12*2.54/100;
    param.tg.Haltitude = 3000*12*2.54/100;
    [param.wfs(1:2).atten] = deal(33);
    [param.wfs(3:3).atten] = deal(39);
  else
    param.tg.altitude_guard = 7500*12*2.54/100;
    param.tg.Haltitude = 12500*12*2.54/100;
    [param.wfs(1:2).atten] = deal(17);
    [param.wfs(3:3).atten] = deal(23);
  end
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('tx_equalization_%.0fus_%.0fAGL.xml', Tpd*1e6, param.tg.Haltitude*100/2.54/12));
  param.prf = 8500;
  param.presums = [9 9 9];
  param.tx_weights = final_DDS_amp;
  param.tukey = 0.1;
  param.Tpd = Tpd;
  for wf=1:3
    param.wfs(wf).phase = final_DDS_phase;
  end
%   for wf=4
%     param.wfs(wf).phase = final_DDS_phase + [0 180 0 0 0 0 0 0];
%   end
  param.delay = final_DDS_time;
  param.f0 = f0;
  param.f1 = f1;
  for wf=1:2
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  for wf=3:3
    param.wfs(wf).tx_mask = [1 1 1 1 1 1 0 0];
  end
  write_cresis_xml(param);
end

%% Equalization High Altitude (Using Ocean), Narrowband
% 10000-45000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
for Tpd = [1e-6 3e-6 10e-6]
  param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
  param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 205; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 315e6; param.DDC_select = 2;
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = false;
  if Tpd == 1e-6
    param.tg.altitude_guard = 2000*12*2.54/100;
    param.tg.Haltitude = 3000*12*2.54/100;
    [param.wfs(1:2).atten] = deal(33);
    [param.wfs(3:3).atten] = deal(39);
  else
    param.tg.altitude_guard = 7500*12*2.54/100;
    param.tg.Haltitude = 12500*12*2.54/100;
    [param.wfs(1:2).atten] = deal(17);
    [param.wfs(3:3).atten] = deal(23);
  end
  param.tg.Hice_thick = 0;
  param.fn = fullfile(base_dir,sprintf('tx_equalization_%.0fus_%.0fAGL_narrowbw.xml', Tpd*1e6, param.tg.Haltitude*100/2.54/12));
  param.prf = 8500;
  param.presums = [9 9 9];
  param.tx_weights = final_DDS_amp;
  param.tukey = 0.1;
  param.Tpd = Tpd;
  for wf=1:3
    param.wfs(wf).phase = final_DDS_phase;
  end
  param.delay = final_DDS_time;
  param.f0 = 180e6;
  param.f1 = 230e6;
  for wf=1:2
    param.wfs(wf).tx_mask = ones(1,8);
    param.wfs(wf).tx_mask(9-wf) = 0;
  end
  for wf=3:3
    param.wfs(wf).tx_mask = [1 1 1 1 1 1 0 0];
  end
  write_cresis_xml(param);
end

% %% Equalization Lab 
% % Use these settings in the lab to demonstrate transmit equalization.
% % Creates one waveform for each of N DDS-transmitters plus a combined
% % waveform with all transmitters going.
% for Tpd = 10e-6
%   param = struct('radar_name','mcords5','num_chan',2,'aux_dac',[255 255 255 255 255 255 255 255],'version','10.0','TTL_prog_delay',650,'fs',1600e6,'fs_sync',50e6,'fs_dds',1000e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9]);
%   param.max_tx = [30000 30000 0 0 0 0 0 0]; param.max_data_rate = 320; param.flight_hours = 7.5; param.sys_delay = 12.9e-6;  param.DDC_freq = 315e6; param.DDC_select = 1;
%   param.max_duty_cycle = 0.12;
%   param.create_IQ = false;
%   param.tg.altitude_guard = 5000*12*2.54/100;
%   param.tg.staged_recording = false;
%   param.tg.Haltitude = 5000*12*2.54/100;
%   param.tg.Hice_thick = 0;
%   param.fn = fullfile(base_dir,sprintf('equalization_%.0fus_LOOPBACK.xml', Tpd*1e6));
%   param.prf = 8500;
%   param.presums = [9 9 9];
%   [param.wfs(1:2).atten] = deal(7);
%   [param.wfs(3:3).atten] = deal(13);
%   param.tx_weights = final_DDS_amp;
%   param.tukey = 0.1;
%   param.Tpd = Tpd;
%   for wf=1:3
%     param.wfs(wf).phase = final_DDS_phase;
%   end
% %   for wf=4
% %     param.wfs(wf).phase = final_DDS_phase + [0 180 0 0 0 0 0 0];
% %   end
%   param.delay = final_DDS_time;
%   param.f0 = f0;
%   param.f1 = f1;
%   for wf=1:2
%     param.wfs(wf).tx_mask = ones(1,8);
%     param.wfs(wf).tx_mask(9-wf) = 0;
%   end
%   for wf=3:3
%     param.wfs(wf).tx_mask = [1 1 1 1 1 1 0 0];
%   end
%   write_cresis_xml(param);
% end
