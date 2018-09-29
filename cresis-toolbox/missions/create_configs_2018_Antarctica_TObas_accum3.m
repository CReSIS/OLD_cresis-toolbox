% script create_configs_2018_Antarctica_TObas_accum3
%
% Creates Arena accumulation radar settings
%
% See link_budgets.xls

%% Output Directories
if ispc
  arena_base_dir = 'C:\waveforms_arena\';
else
  arena_base_dir = '/arena/waveforms/';
end
calval_dir = 'cal_val';

%% Frequency ranges, DDC settings
f0_list = [600e6];
f1_list = [900e6];
cal_settings_list = [1];
prf = 1/40e-6;

%% Calculate the desired presums
% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
physical_constants; % c = speed of light
velocity = [70]; % m/s
quarter_wavelength = c./max(abs(f0_list),abs(f1_list)) / 4;
presums = round(quarter_wavelength ./ velocity * prf);
fprintf('For prf %g Hz, max presums for quarter wavelength (%g m) sampling is: %g presums\n', prf, quarter_wavelength, presums);

%% Calibration Settings for AWG/DAC
final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 0
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0];
    final_DDS_phase_no_time{idx} = [0]; % not used usually
    final_DDS_amp{idx} = [0.7];
    final_DDS_time{idx} =  [0];
  end
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % 600-900 MHz
  idx = 1;
  final_DDS_phase{idx} = [0];
  final_DDS_phase_no_time{idx} = [0]; % not used usually
  final_DDS_amp{idx} = [0.7];
  final_DDS_time{idx} =  [0];
end

% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = chebwin(1).';

%% Survey Mode + survey-test, loopback, noise, and deconv modes
% <1700 m thick ice, 1200 +/- 500 ft AGL
ice_thickness = [1600];
for freq_idx = 1
  [~,defaults] = default_radar_params_2018_Antarctica_TObas_accum3;
  param = defaults{1};
  param.flight_hours = 4;
  param.prf = 1/40e-6;
  
  param.tg.staged_recording = {[1],[2]};
  param.tg.altitude_guard = 300*12*2.54/100;
  param.tg.Haltitude = 1400*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  
  param.create_IQ = false;
  
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = 250e6;
  param.zeropimods = [0 180];
  
  param.wfs(1).presums = 48;
  param.wfs(1).Tpd = 2e-6;
  param.wfs(1).tukey = 0.1;
  
  param.wfs(1).name = 'single';
  
  param.tx_weights = final_DDS_amp{cal_settings_list(freq_idx)};
  param.wfs(1).phase = final_DDS_phase{cal_settings_list(freq_idx)};
  param.phase = final_DDS_phase{cal_settings_list(freq_idx)};
  param.wfs(1).delay = final_DDS_time{cal_settings_list(freq_idx)};
  param.delay = final_DDS_time{cal_settings_list(freq_idx)};
  
  idx = find(strcmpi('AttenFirst18dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  param.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  
  idx = find(strcmpi('AttenSecond7dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[1 1]};
  param.arena.ctu.out.bit_group(idx).pri = {[1 1]};
  
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick);
  param.arena.fn = fullfile(arena_base_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  
  % Survey Test Mode (begins recording at 0 ft altitude)
  old_param = param;
  param.tg.Haltitude = param.tg.altitude_guard; % Start time 0 sec
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick_TEST', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick);
  param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  param = old_param;
  
  % Loopback Mode (0e-6 delay line)
  old_param = param;
  param.tg.Haltitude = 3e8/2 * 0e-6; % Set the delay line here
  param.tg.altitude_guard = 3e8/2 * 0e-6;
  param.tg.Hice_thick = 800;
  param.tg.staged_recording = {[0],[0]};
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fusDelay_%.0fus_LOOPBACK', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude/(3e8/2)*1e6,param.wfs(end).Tpd*1e6);
  param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  param = old_param;
  
  % Deconvolution Mode (for over calm lake or sea ice lead)
  old_param = param;
  idx = find(strcmpi('AttenFirst18dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  param.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  idx = find(strcmpi('AttenSecond7dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  param.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  param.tg.staged_recording = {[0],[0]};
  param.tg.altitude_guard = 3000*12*2.54/100;
  param.tg.Haltitude = 4000*12*2.54/100;
  param.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6);
  param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  param = old_param;
  
  % Noise Mode
  old_param = param;
  param.tx_weights = [0];
  param.tx_enable = [0];
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_NOISE', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6);
  param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  param = old_param;
end


%% Survey Mode + survey-test, loopback, noise, and deconv modes
% <1700 m thick ice, 1200 +/- 500 ft AGL
ice_thickness = [2500];
for freq_idx = 1
  [~,defaults] = default_radar_params_2018_Antarctica_TObas_accum3;
  param = defaults{1};
  param.flight_hours = 4;
  param.prf = 1/50e-6;
  
  param.tg.staged_recording = {[1],[2]};
  param.tg.altitude_guard = 300*12*2.54/100;
  param.tg.Haltitude = 1400*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  
  param.create_IQ = false;
  
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = 250e6;
  param.zeropimods = [0 180];
  
  param.wfs(1).presums = 60;
  param.wfs(1).Tpd = 2e-6;
  param.wfs(1).tukey = 0.1;
  
  param.wfs(1).name = 'single';
  
  param.tx_weights = final_DDS_amp{cal_settings_list(freq_idx)};
  param.wfs(1).phase = final_DDS_phase{cal_settings_list(freq_idx)};
  param.phase = final_DDS_phase{cal_settings_list(freq_idx)};
  param.wfs(1).delay = final_DDS_time{cal_settings_list(freq_idx)};
  param.delay = final_DDS_time{cal_settings_list(freq_idx)};
  
  idx = find(strcmpi('AttenFirst18dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  param.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  
  idx = find(strcmpi('AttenSecond7dB',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[1 1]};
  param.arena.ctu.out.bit_group(idx).pri = {[1 1]};
  
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick);
  param.arena.fn = fullfile(arena_base_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  
  % Survey Test Mode (begins recording at 0 ft altitude)
  old_param = param;
  param.tg.Haltitude = param.tg.altitude_guard; % Start time 0 sec
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick_TEST', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick);
  param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  param = old_param;
end
