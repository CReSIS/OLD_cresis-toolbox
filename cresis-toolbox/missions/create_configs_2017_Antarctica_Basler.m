% script create_configs_2017_Antarctica_Basler
%
% Creates NI radar depth sounder settings
%
% See link_budgets.xls

physical_constants; % c = speed of light

% Define waveforms
if ispc
  base_dir = 'c:\waveforms\';
  rss_base_dir = 'c:\Temp\';
else
  base_dir = '~/waveforms/';
  rss_base_dir = 'c:\Temp\';
end

f0_list = [150e6 180e6];
f1_list = [450e6 210e6];
DDC_select_list = [1 2]; % Which DDC mode to use
cal_settings = [1 2];
prf = 12000;

% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
velocity = [65 90]; % m/s
presums = round(c./abs(f1_list+f0_list)/2 ./ velocity * prf / 2)*2

final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 1
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0 0 0 0 0 0 0 0];
    final_DDS_phase_no_time = [0 0 0 0 0 0 0 0]; % not used usually
    final_DDS_amp{idx} = [40000 40000 40000 40000 40000 40000 40000 40000];
    final_DDS_time{idx} =  [0 0 0 0 0 0 0 0];
  end
else
  % COPY AND PASTE RESULTS FROM run_basic_tx_chan_equalization.m
  % HERE:
  
  % These are from transmit calibration during 2017? test flight  
  
  % 150-450 MHz
  final_DDS_phase{end+1} = [0 0 0 0 0 0 0 0];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [0 0 0 0 0 0 0 0];
  final_DDS_time{end+1} =  [0 0 0 0 0 0 0 0];
  
  % 180-210 MHz
  final_DDS_phase{end+1} = [-114.5	-152.8	-165.3	0.0	-170.0	0.0	-177.9	-164.6];
  final_DDS_phase_no_time{end+1} = [0 0 0 0 0 0 0 0]; % not used usually
  final_DDS_amp{end+1} = [34461	31170	32120	39658	34905	0	40000	38067];
  final_DDS_time{end+1} =  [-1.86	-2.43	-2.47	0.00	-2.51	0.00	-2.64	-2.50]*1e-9;
end

% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = chebwin(8,30).';

%% SETUP
% =========================================================================

param = [];
param.season_name = '2017_Antarctica_Basler';
param.radar_name = 'mcords5';
param.gps_source = 'cresis-field';
clear phase_centers;
for tx_chan = 1:8
  % Just enable the current antenna to determine its phase center
  tx_weights = zeros(1,8);
  tx_weights(tx_chan) = 1;
  rxchan = 4; % Fix the receiver (it should not matter which one you choose)
  % Determine phase center for the antenna
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, rxchan);
end
% Adjust phase centers to the mean phase center position
phase_centers = bsxfun(@minus,phase_centers,mean(phase_centers,2));

% Create waveforms directories if they do not exist
if ~exist(base_dir,'dir')
  mkdir(base_dir);
end
calval_dir = fullfile(base_dir, 'calval');
if ~exist(calval_dir,'dir')
  mkdir(calval_dir);
end

%% AWI MCoRDS Arena Parameters
arena.awg = [];
arena.awg(end+1).awg = 0;
arena.awg(end).dacs = [0 1];
arena.awg(end).dacClk = [1600e6 1600e6];
arena.awg(end).desiredAlignMin = [0 20];
arena.awg(end).desiredAlignMax = [10 30];
arena.awg(end+1).awg = 1;
arena.awg(end).dacs = [2 3];
arena.awg(end).dacClk = [1600e6 1600e6];
arena.awg(end).desiredAlignMin = [-5 15];
arena.awg(end).desiredAlignMax = [10 30];
arena.awg(end+1).awg = 2;
arena.awg(end).dacs = [4 5];
arena.awg(end).dacClk = [1600e6 1600e6];
arena.awg(end).desiredAlignMin = [0 20];
arena.awg(end).desiredAlignMax = [10 30];
arena.awg(end+1).awg = 3;
arena.awg(end).dacs = [6 7];
arena.awg(end).dacClk = [1600e6 1600e6];
arena.awg(end).desiredAlignMin = [0 20];
arena.awg(end).desiredAlignMax = [10 30];
arena.dacs = [0 1 2 3 4 5 6 7];
arena.dacs_sampFreq = [1600e6 1600e6 1600e6 1600e6 1600e6 1600e6 1600e6 1600e6];
arena.max_tx = [0.63 0.63 0.63 0.63 0.63 0.63 0.63 0.63];
arena.zeropimods = [0 180];
arena.TTL_time = [0.1 0.2 2.2];

arena.TTL_names = {};
for PA = 1:8
  arena.TTL_names{end+1} = sprintf('PA ENA %d',PA);
end
arena.TTL_names{end+1} = 'T/R';
arena.TTL_names{end+1} = 'ISO';
arena.TTL_names{end+1} = 'EPRI';
arena.TTL_names{end+1} = 'PRI';
arena.TTL_names{end+1} = 'EPRI';
arena.TTL_names{end+1} = 'PRI';
arena.TTL_names{end+1} = 'EPRI';
arena.TTL_names{end+1} = 'PRI';
arena.TTL_states{1} = [
  0 1 1 0 % PA ENA 1
  0 1 1 0 % PA ENA 2
  0 1 1 0 % PA ENA 3
  0 1 1 0 % PA ENA 4
  0 1 1 0 % PA ENA 5
  0 1 1 0 % PA ENA 6
  0 1 1 0 % PA ENA 7
  0 1 1 0 % PA ENA 8
  0 1 1 0 % T/R
  0 1 1 0 % ISO
  1 0 0 0 % EPRI
  0 1 0 0 % PRI
  1 0 0 0 % EPRI
  0 1 0 0 % PRI
  1 0 0 0 % EPRI
  0 1 0 0 % PRI
  ];
arena.TTL_states{2} = [
  0 1 1 0 % PA ENA 1
  0 1 1 0 % PA ENA 2
  0 1 1 0 % PA ENA 3
  0 1 1 0 % PA ENA 4
  0 1 1 0 % PA ENA 5
  0 1 1 0 % PA ENA 6
  0 1 1 0 % PA ENA 7
  0 1 1 0 % PA ENA 8
  0 1 1 0 % T/R
  0 1 1 0 % ISO
  0 0 0 0 % EPRI
  0 1 0 0 % PRI
  0 0 0 0 % EPRI
  0 1 0 0 % PRI
  0 0 0 0 % EPRI
  0 1 0 0 % PRI
  ];

%% Survey Mode + loopback, noise, and deconv modes
% <4000 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [4000 4000];
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
  param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 755; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = 1 + [2 4 presums(freq_idx)-6];
  param.wfs(1).atten = 37;
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
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 0e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK_NO_DELAY.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Loopback Mode (10e-6 delay line)
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 1000*12*2.54/100;
  param.tg.Haltitude = 10e-6 * c/2;
  param.tg.Hice_thick = 0; % Long enough for 10 us delay line
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_LOOPBACK.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  % Deconvolution Mode (for over calm lake or sea ice lead)
  param.wfs(1).atten = 43;
  param.wfs(2).atten = 43;
  param.wfs(3).atten = 43;
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 3000*12*2.54/100;
  param.tg.Haltitude = 4000*12*2.54/100;
  param.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
  param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
  if freq_idx == 1
    % Noise Mode
    param.tx_weights = [0 0 0 0 0 0 0 0];
    [param.wfs(1:3).tx_mask] = deal([1 1 1 1 1 1 1 1]);
    param.wfs(1).atten = 37;
    param.wfs(2).atten = 0;
    param.wfs(3).atten = 0;
    param.tg.staged_recording = [1 2 3];
    param.tg.altitude_guard = 500*12*2.54/100;
    param.tg.Haltitude = 1400*12*2.54/100;
    param.tg.Hice_thick = 3250;
    param.fn = fullfile(calval_dir,sprintf('survey_%.0f-%.0fMHz_%.0fus_NOISE.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
    write_cresis_xml(param);
  end
end


%% Survey Mode for thin ice
% <2500 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [2500 2500];
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
  param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.rg_stop_offset = [500*sqrt(3.15) 0 0]; % Keep waveform 1 on for 500 m of ice
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = 1 + [8 2 presums(freq_idx)-10];
  param.wfs(1).atten = 35;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(3).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(3).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.fn = fullfile(base_dir,sprintf('thinice_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end


%% Sea Ice
% 1200 +/- 1200 ft AGL
for freq_idx = [1]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
  param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [0 0];
  param.tg.rg_stop_offset = [0 0]; % Keep waveform 1 on for 500 m of ice
  param.tg.altitude_guard = 1200*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.prf = prf;
  param.presums = 1 + [presums(freq_idx)/2 presums(freq_idx)/2];
  param.wfs(1).atten = 13;
  param.wfs(2).atten = 13;
  DDS_amp = final_DDS_amp{cal_settings(freq_idx)};
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  param.wfs(1).tx_mask = deal([0 1 1 1 1 1 1 1]);
  param.wfs(2).tx_mask = deal([1 1 1 1 1 1 1 0]);
  param.fn = fullfile(base_dir,sprintf('seaice_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Image Mode (Low Altitude, Thick Ice)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
  param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
  param.max_data_rate = 755;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [0 0 0];
  param.tg.start_ref = {'bottom','surface','bottom'};
  param.tg.stop_ref = {'bottom','surface','bottom'};
  param.tg.altitude_guard = 700 * 12*2.54/100;
  param.tg.Haltitude = 1200 * 12*2.54/100;
  param.tg.Hice_thick_min = 1500;
  param.tg.Hice_thick = 3000;
  param.tg.look_angle_deg = [40 0 40];
  param.prf = prf;
  param.presums = 1 + [ceil(presums(freq_idx)/4)*2 2 ceil(presums(freq_idx)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_DDS_amp);
  DDS_amp = round(DDS_amp .* param.max_DDS_amp(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 10e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  param.wfs(3).phase = final_DDS_phase{freq_idx};
  % Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
  beam_angle_deg = 20; % Positive to the left
  param.wfs(1).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = 0; % Nadir
  param.wfs(2).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = -20; % Negative to the right
  param.wfs(3).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.wfs(1).atten = 0;
  param.wfs(2).atten = 37;
  param.wfs(3).atten = 0;
  param.fn = fullfile(base_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Image Mode Pattern Measurements
% 3500 ft +/- 1000 ft AGL
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
  param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
  param.max_data_rate = 755;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [0 0];
  param.tg.start_ref = {'surface','surface'};
  param.tg.altitude_guard = 1000 * 12*2.54/100;
  param.tg.Haltitude = 3500 * 12*2.54/100;
  param.tg.Hice_thick_min = 0;
  param.tg.Hice_thick = 0;
  param.tg.look_angle_deg = [0 0 0];
  param.prf = prf;
  param.presums = 1 + [ceil(presums(freq_idx)/4)*2 ceil(presums(freq_idx)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_DDS_amp);
  DDS_amp = round(DDS_amp .* param.max_DDS_amp(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 3e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  % Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
  beam_angle_deg = 20; % Positive to the left
  param.wfs(1).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = -20; % Negative to the right
  param.wfs(2).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:2).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  [param.wfs(1:2).atten] = deal(43);
  param.fn = fullfile(calval_dir,sprintf('image_%.0f-%.0fMHz_%.0fft_%.0fus_PATTERN.xml', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6));
  write_cresis_xml(param);
end

%% Image Mode (High Altitude, Thin Ice)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
freq_idx_WB = 1;
freq_idx_NB = 2;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
param.max_data_rate = 755;
param.DDC_select = DDC_select_list(freq_idx_WB);
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = [0 0 0];
param.tg.start_ref = {'bottom','surface','bottom'};
param.tg.stop_ref = {'bottom','bottom','bottom'};
param.tg.altitude_guard = 1000 * 12*2.54/100;
param.tg.Haltitude = 6000 * 12*2.54/100;
param.tg.Hice_thick_min = 0;
param.tg.Hice_thick = 1000;
param.tg.look_angle_deg = [40 0 40];
param.prf = prf;
param.presums = 1 + [ceil(presums(freq_idx_WB)/4)*2 4 ceil(presums(freq_idx_WB)/4)*2];
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{freq_idx_WB} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_DDS_amp);
DDS_amp = round(DDS_amp .* param.max_DDS_amp(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 1e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{freq_idx_NB};
param.wfs(2).phase = final_DDS_phase{freq_idx_WB};
param.wfs(3).phase = final_DDS_phase{freq_idx_NB};
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
beam_angle_deg = 20; % Positive to the left
param.wfs(1).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = 0; % Nadir
param.wfs(2).delay = final_DDS_time{freq_idx_WB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(3).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.wfs(1).f0 = f0_list(freq_idx_NB);
param.wfs(2).f0 = f0_list(freq_idx_WB);
param.wfs(3).f0 = f0_list(freq_idx_NB);
param.wfs(1).f1 = f1_list(freq_idx_NB);
param.wfs(2).f1 = f1_list(freq_idx_WB);
param.wfs(3).f1 = f1_list(freq_idx_NB);
param.DDC_freq = (param.wfs(2).f0+param.wfs(2).f1)/2;
[param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.wfs(1).atten = 0;
param.wfs(2).atten = 23;
param.wfs(3).atten = 0;
param.fn = fullfile(base_dir,sprintf('imagehighthin_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
  param.wfs(2).f0/1e6,param.wfs(2).f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
write_cresis_xml(param);

%% Equalization (Using Ocean)
% Haltitude +/- 1000 ft
% For lower altitude, increase attenuation
% Use these settings over ocean or sea ice for fast-time equalization,
% transmit equalization, and receiver equalization.
% Creates one waveform for each of N DDS-transmitters plus a combined
% waveform with all transmitters going.
Haltitude = [1500 1500 0 3000 6000];
Tpd_list = [1e-6 1e-6 3e-6 3e-6 3e-6];
attenuation = [43 39 0 43 43];
fn_hint = {'WATER','ICE','NO_DELAY','WATER','WATER'};
for Tpd_idx = 1:length(Tpd_list)
  Tpd = Tpd_list(Tpd_idx);
  for freq_idx = [1 2]
    param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
    param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
    param.DDC_select = DDC_select_list(freq_idx);
    param.max_duty_cycle = 0.12;
    param.create_IQ = false;
    param.tg.staged_recording = false;
    param.tg.altitude_guard = 2000*12*2.54/100;
    param.tg.Haltitude = Haltitude(Tpd_idx)*12*2.54/100;
    param.tg.Hice_thick = 0;
    param.prf = prf;
    param.presums = 1 + [10 10 10 10 10 10 10 10];
    [param.wfs(1:8).atten] = deal(max(0,attenuation(Tpd_idx)-12));
%     [param.wfs(9:9).atten] = deal(attenuation(Tpd_idx));
    param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
    param.tukey = 0.08;
    param.Tpd = Tpd;
    for wf=1:length(param.wfs)
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
%     for wf=9:9
%       param.wfs(wf).tx_mask = [0 0 0 0 0 0 0 0];
%     end
    param.fn = fullfile(calval_dir,sprintf('txequal_%.0f-%.0fMHz_%.0fft_%.0fus_%s.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.Tpd*1e6,fn_hint{Tpd_idx}));
    write_cresis_xml(param);
  end
end

%% Max power mode with max frequency range (only for EMI survey)
freq_idx = 1;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',50.0e6,'fs_dds',1000e6,'TTL_clock',50e6,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'rss_base_dir',rss_base_dir);
param.max_DDS_amp = [40000 40000 40000 40000 40000 40000 40000 40000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 13.25e-6; param.use_mcords4_names = true;
param.DDC_select = DDC_select_list(freq_idx);
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = false;
param.tg.altitude_guard = 1000*12*2.54/100;
param.tg.Haltitude = 1400*12*2.54/100;
param.tg.Hice_thick = 0;
param.prf = prf;
param.presums = 1 + presums(freq_idx);
param.wfs(1).atten = 43;
param.tukey = 0;
param.wfs(1).Tpd = 10e-6;
param.wfs(1).phase = [0 0 0 0 0 0 0 0];
param.delay = [0 0 0 0 0 0 0 0];
param.f0 = f0_list(freq_idx);
param.f1 = f1_list(freq_idx);
param.DDC_freq = (param.f0+param.f1)/2;
[param.wfs(1:1).tx_mask] = deal([1 1 1 1 1 1 1 1]);
param.tx_weights = [40000 40000 40000 40000 40000 40000 40000 40000] * sqrt(0);
param.fn = fullfile(calval_dir,sprintf('singlewf_%.0f-%.0fMHz_%.0fus_TX_OFF.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
write_cresis_xml(param);
[param.wfs(1:1).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.tx_weights = [40000 40000 40000 40000 40000 40000 40000 40000] * sqrt(1.00);
param.fn = fullfile(base_dir,sprintf('singlewf_%.0f-%.0fMHz_%.0fus_DDS_CHECK.xml',param.f0/1e6,param.f1/1e6,param.wfs(end).Tpd*1e6));
write_cresis_xml(param);

