% script create_settings_2018_Antarctica_TObas
%
% Creates Arena accumulation radar settings
%
% See link_budgets.xls

physical_constants; % c = speed of light

% Define waveforms
if ispc
  arena_base_dir = 'C:\waveforms_arena\';
else
  arena_base_dir = '/arena/waveforms/';
end

f0_list = [600e6];
f1_list = [900e6];
DDC_select_list = [1]; % Which DDC mode to use
cal_settings = [1];
prf = 10000;

% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
velocity = [70]; % m/s
presums = round(c./max(abs(f0_list),abs(f1_list)) ./ velocity * prf / 4)*4

final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 0
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0];
    final_DDS_phase_no_time{idx} = [0]; % not used usually
    final_DDS_amp{idx} = [1];
    final_DDS_time{idx} =  [0];
  end
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % 600-900 MHz
  idx = 1;
  final_DDS_phase{idx} = [0];
  final_DDS_phase_no_time{idx} = [0]; % not used usually
  final_DDS_amp{idx} = [1];
  final_DDS_time{idx} =  [0];
end

% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = chebwin(1).';

%% SETUP
% =========================================================================

param = [];
param.season_name = '2018_Antarctica_TObas';
param.radar_name = 'accum3';
param.gps_source = 'bas-final';
clear phase_centers;
for tx_chan = 1:1
  % Just enable the current antenna to determine its phase center
  tx_weights = zeros(1,1);
  tx_weights(tx_chan) = 1;
  rxchan = 1; % Fix the receiver (it should not matter which one you choose)
  % Determine phase center for the antenna
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, rxchan);
end
% Adjust phase centers to the mean phase center position
phase_centers = bsxfun(@minus,phase_centers,mean(phase_centers,2));

%% Survey Mode + loopback, noise, and deconv modes
% <1700 m thick ice, 1200 +/- 500 ft AGL
ice_thickness = [1700];
for freq_idx = [1]
  [~,defaults] = default_radar_params_2018_Antarctica_TObas;
  param = defaults{1};
  param.flight_hours = 4;
  param.prf = prf;
  
  param.zeropimods = [0 180]; % 180 causes a special mode in some ADC/DAC
  
  param.tg.staged_recording = [1];
  param.tg.altitude_guard = 500*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  
  param.create_IQ = false;
  
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = 250e6;
  param.zeropimods = [0 180];
  
  % wf == 1 is the EPRI always
  param.interleave = 7; % Explicitly programs 2 times as many sequence slots (first set is for interleave loop, second set is for calibration loop)
  param.wfs(1).presums = 4; % Creates 1+Z modes (2 would create 2 modes)
  param.wfs(2).presums = 14; % Creates Z modes
%   param.cal.wfs(1).presums = 7*(4+14)-2; % Create 1+Z modes (epri + pri)
%   param.cal.wfs(2).presums = 1; % Create 1
%   param.cal.wfs(3).presums = 1; % Create 1
  param.tukey = 0.08;
  param.wfs(1).Tpd = 3e-6;
  param.wfs(2).Tpd = 10e-6;
  param.wfs(1).tukey = 0.1;
  param.wfs(2).tukey = 0.1;
  
  param.wfs(1).name = 'low_gain';
  param.wfs(2).name = 'high_gain';
  
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(2).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(1).delay = final_DDS_time{cal_settings(freq_idx)};
  param.wfs(2).delay = final_DDS_time{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  
  idx = find(strcmpi('AttenFirst18dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[1 1],[0 0]};
  param.arena.ctu.out.bit_group(idx).pri = {[1 1],[0 0]};
  
  idx = find(strcmpi('AttenSecond7dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[0 0],[0 0]};
  param.arena.ctu.out.bit_group(idx).pri = {[0 0],[0 0]};
  
  param.arena.fn = fullfile(arena_base_dir,sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_radar_config(param);

  return
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

%% Survey Mode
% <3250 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [3250 3250];
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [2 4 presums(freq_idx)-6];
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
end

%% Survey Mode for thin ice
% <2500 m thick ice, 1200 +/- 700 ft AGL
ice_thickness = [2500 2500];
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3];
  param.tg.rg_stop_offset = [500*sqrt(3.15) 0 0]; % Keep waveform 1 on for 500 m of ice
  param.tg.altitude_guard = 700*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  param.prf = prf;
  param.presums = [8 2 presums(freq_idx)-10];
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
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 750; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [0 0];
  param.tg.rg_stop_offset = [0 0]; % Keep waveform 1 on for 500 m of ice
  param.tg.altitude_guard = 1200*12*2.54/100;
  param.tg.Haltitude = 1200*12*2.54/100;
  param.tg.Hice_thick = 0;
  param.prf = prf;
  param.presums = [round(presums(freq_idx)/2/2)*2 round(presums(freq_idx)/2/2)*2];
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

%% Image Mode (EGRIP Low Altitude, Thick Ice)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
for freq_idx = 2
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 4; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
  param.max_data_rate = 755;
  param.DDC_select = DDC_select_list(freq_idx);
  param.max_duty_cycle = 0.12;
  param.create_IQ = false;
  param.tg.staged_recording = [1 2 3 3];
  param.tg.altitude_guard = 700 * 12*2.54/100;
  param.tg.Haltitude = 1200 * 12*2.54/100;
  param.tg.Hice_thick = 3250;
  param.tg.look_angle_deg = [0 0 40 40];
  param.prf = prf;
  param.presums = [2 4 ceil((presums(freq_idx)-6)/4)*2 ceil((presums(freq_idx)-6)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_tx);
  DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 1e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(4).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  param.wfs(3).phase = final_DDS_phase{freq_idx};
  param.wfs(4).phase = final_DDS_phase{freq_idx};
  param.wfs(1).name = 'nadir1';
  param.wfs(2).name = 'nadir3';
  param.wfs(3).name = 'left';
  param.wfs(4).name = 'right';
  % Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
  beam_angle_deg = 0; % Nadir
  param.wfs(1).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = 0; % Nadir
  param.wfs(2).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = 20; % Positive to the left
  param.wfs(3).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  beam_angle_deg = -20; % Negative to the right
  param.wfs(4).delay = final_DDS_time{freq_idx} ...
    - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
    - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
    + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = (param.f0+param.f1)/2;
  [param.wfs(1:4).tx_mask] = deal([0 0 0 0 0 0 0 0]);
  param.wfs(1).atten = 37;
  param.wfs(2).atten = 0;
  param.wfs(3).atten = 0;
  param.wfs(4).atten = 0;
  param.fn = fullfile(base_dir,sprintf('egrip_image_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
  write_cresis_xml(param);
end

%% Image Mode (Low Altitude, Thick Ice)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
for freq_idx = [1 2]
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
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
  param.presums = [ceil(presums(freq_idx)/4)*2 2 ceil(presums(freq_idx)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_tx);
  DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 10e-6;
  param.wfs(2).Tpd = 1e-6;
  param.wfs(3).Tpd = 10e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  param.wfs(3).phase = final_DDS_phase{freq_idx};
  param.wfs(1).name = 'left';
  param.wfs(2).name = 'nadir';
  param.wfs(3).name = 'right';
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
  param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
  param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
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
  param.presums = [ceil(presums(freq_idx)/4)*2 ceil(presums(freq_idx)/4)*2];
  % Switch from tx calibration window to hanning window to broaden beam
  DDS_amp = final_DDS_amp{freq_idx} .* hanning(8).' ./ Hwindow_orig;
  % Renormalize the amplitudes
  [~,relative_max_idx] = max(DDS_amp./param.max_tx);
  DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
  param.tx_weights = DDS_amp;
  param.tukey = 0.08;
  param.wfs(1).Tpd = 3e-6;
  param.wfs(2).Tpd = 3e-6;
  param.wfs(1).phase = final_DDS_phase{freq_idx};
  param.wfs(2).phase = final_DDS_phase{freq_idx};
  param.wfs(1).name = 'left';
  param.wfs(2).name = 'right';
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
freq_idx_NB = 3;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
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
param.presums = [ceil(presums(freq_idx_WB)/4)*2 4 ceil(presums(freq_idx_WB)/4)*2];
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{freq_idx_WB} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 1e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{freq_idx_NB};
param.wfs(2).phase = final_DDS_phase{freq_idx_WB};
param.wfs(3).phase = final_DDS_phase{freq_idx_NB};
param.wfs(1).name = 'left';
param.wfs(2).name = 'nadir';
param.wfs(3).name = 'right';
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

%% Image Mode (Low Altitude, Ice <2000 m thick)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
freq_idx_WB = 1;
freq_idx_NB = 2;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
param.max_data_rate = 755;
param.DDC_select = DDC_select_list(freq_idx_NB);
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = [0 0 0];
param.tg.start_ref = {'bottom','surface','bottom'};
param.tg.stop_ref = {'bottom','bottom','bottom'};
param.tg.altitude_guard = 500 * 12*2.54/100;
param.tg.Haltitude = 1200 * 12*2.54/100;
param.tg.Hice_thick_min = 0;
param.tg.Hice_thick = 2000;
param.tg.look_angle_deg = [40 0 40];
param.prf = prf;
param.presums = [8 8 8]; %[ceil(presums(freq_idx_WB)/4)*2 4 ceil(presums(freq_idx_WB)/4)*2];
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{freq_idx_NB} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{freq_idx_NB};
param.wfs(2).phase = final_DDS_phase{freq_idx_NB};
param.wfs(3).phase = final_DDS_phase{freq_idx_NB};
param.wfs(1).name = 'left';
param.wfs(2).name = 'nadir';
param.wfs(3).name = 'right';
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
beam_angle_deg = 20; % Positive to the left
param.wfs(1).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = 0; % Nadir
param.wfs(2).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(3).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.wfs(1).f0 = f0_list(freq_idx_NB);
param.wfs(2).f0 = f0_list(freq_idx_NB);
param.wfs(3).f0 = f0_list(freq_idx_NB);
param.wfs(1).f1 = f1_list(freq_idx_NB);
param.wfs(2).f1 = f1_list(freq_idx_NB);
param.wfs(3).f1 = f1_list(freq_idx_NB);
param.DDC_freq = (param.wfs(2).f0+param.wfs(2).f1)/2;
[param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
param.fn = fullfile(base_dir,sprintf('image_angelika_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
  param.wfs(2).f0/1e6,param.wfs(2).f1/1e6,param.tg.Haltitude*100/2.54/12,param.wfs(end).Tpd*1e6,param.tg.Hice_thick));
write_cresis_xml(param);

%% Image Mode (Low Altitude, Ice <1400 m thick)
% Ice thickness "param.tg.Hice_thick_min" m to "param.tg.Hice_thick" m, "param.tg.Haltitude" +/- "param.tg.altitude_guard" ft AGL
freq_idx_WB = 1;
freq_idx_NB = 2;
param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
param.max_data_rate = 755;
param.DDC_select = DDC_select_list(freq_idx_NB);
param.max_duty_cycle = 0.12;
param.create_IQ = false;
param.tg.staged_recording = [0 0 0];
param.tg.start_ref = {'bottom','surface','bottom'};
param.tg.stop_ref = {'bottom','bottom','bottom'};
param.tg.altitude_guard = 500 * 12*2.54/100;
param.tg.Haltitude = 1200 * 12*2.54/100;
param.tg.Hice_thick_min = 0;
param.tg.Hice_thick = 1400;
param.tg.look_angle_deg = [40 0 40];
param.prf = prf;
param.presums = [8 8 8]; %[ceil(presums(freq_idx_WB)/4)*2 4 ceil(presums(freq_idx_WB)/4)*2];
% Switch from tx calibration window to hanning window to broaden beam
DDS_amp = final_DDS_amp{freq_idx_NB} .* hanning(8).' ./ Hwindow_orig;
% Renormalize the amplitudes
[~,relative_max_idx] = max(DDS_amp./param.max_tx);
DDS_amp = round(DDS_amp .* param.max_tx(relative_max_idx) / DDS_amp(relative_max_idx));
param.tx_weights = DDS_amp;
param.tukey = 0.08;
param.wfs(1).Tpd = 3e-6;
param.wfs(2).Tpd = 3e-6;
param.wfs(3).Tpd = 3e-6;
param.wfs(1).phase = final_DDS_phase{freq_idx_NB};
param.wfs(2).phase = final_DDS_phase{freq_idx_NB};
param.wfs(3).phase = final_DDS_phase{freq_idx_NB};
param.wfs(1).name = 'left';
param.wfs(2).name = 'nadir';
param.wfs(3).name = 'right';
% Add in time delays to each position, subtract out the nadir time delays since tx_equalization already took care of those
beam_angle_deg = 20; % Positive to the left
param.wfs(1).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = 0; % Nadir
param.wfs(2).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
beam_angle_deg = -20; % Negative to the right
param.wfs(3).delay = final_DDS_time{freq_idx_NB} ...
  - (phase_centers(2,:) / (c/2) * sind(beam_angle_deg))*1e9 ...
  - (phase_centers(3,:) / (c/2) * cosd(beam_angle_deg))*1e9 ...
  + (phase_centers(3,:) / (c/2) * cosd(0))*1e9;
param.wfs(1).f0 = f0_list(freq_idx_NB);
param.wfs(2).f0 = f0_list(freq_idx_NB);
param.wfs(3).f0 = f0_list(freq_idx_NB);
param.wfs(1).f1 = f1_list(freq_idx_NB);
param.wfs(2).f1 = f1_list(freq_idx_NB);
param.wfs(3).f1 = f1_list(freq_idx_NB);
param.DDC_freq = (param.wfs(2).f0+param.wfs(2).f1)/2;
[param.wfs(1:3).tx_mask] = deal([0 0 0 0 0 0 0 0]);
param.wfs(1).atten = 0;
param.wfs(2).atten = 0;
param.wfs(3).atten = 0;
param.fn = fullfile(base_dir,sprintf('image_angelika_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick.xml', ...
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
attenuation = [43 39 43 43 43];
fn_hint = {'WATER','ICE','NO_DELAY','WATER','WATER'};
for Tpd_idx = 1:length(Tpd_list)
  Tpd = Tpd_list(Tpd_idx);
  for freq_idx = [1 2]
    param = struct('radar_name','mcords5','num_chan',8,'aux_dac',[255 255 255 255 255 255 255 255],'version','14.0f1','TTL_prog_delay',650,'xml_version',2.0,'fs',1600e6,'fs_sync',90.0e6,'fs_dds',1440e6,'TTL_clock',1440e6/16,'TTL_mode',[2.5e-6 260e-9 -1100e-9],'arena_base_dir',arena_base_dir);
    param.max_tx = [4000 4000 4000 4000 4000 4000 4000 4000]; param.max_data_rate = 700; param.flight_hours = 3.5; param.sys_delay = 0.75e-6; param.use_mcords4_names = true; param.arena = arena;
    param.DDC_select = DDC_select_list(freq_idx);
    param.max_duty_cycle = 0.12;
    param.create_IQ = false;
    param.tg.staged_recording = false;
    param.tg.altitude_guard = 1000*12*2.54/100;
    param.tg.Haltitude = Haltitude(Tpd_idx)*12*2.54/100;
    param.tg.Hice_thick = 0;
    param.prf = prf;
    param.presums = [10 10 10 10 10 10 10 10 10];
    [param.wfs(1:8).atten] = deal(attenuation(Tpd_idx)-12);
    [param.wfs(9:9).atten] = deal(attenuation(Tpd_idx));
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
    param.fn = fullfile(calval_dir,sprintf('txequal_%.0f-%.0fMHz_%.0fft_%.0fus_%s.xml',param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.Tpd*1e6,fn_hint{Tpd_idx}));
    write_cresis_xml(param);
  end
end
