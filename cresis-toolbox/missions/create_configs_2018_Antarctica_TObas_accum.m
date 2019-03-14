% script create_configs_2018_Antarctica_TObas_accum
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

%% Survey Mode + survey-test, loopback, noise, and deconv modes
% <1700 m thick ice, 1400 +/- 300 ft AGL
ice_thickness = [1600];
for freq_idx = 1
  param = default_radar_params_2018_Antarctica_TObas_accum;
  config = param.config;
  
  config.flight_hours = 4;
  config.prf = 1/40e-6;
  
  config.tg.staged_recording = {[1],[2]};
  config.tg.altitude_guard = 300*12*2.54/100;
  config.tg.Haltitude = 1400*12*2.54/100;
  config.tg.Hice_thick = ice_thickness(freq_idx);
  
  config.create_IQ = false;
  
  config.f0 = f0_list(freq_idx);
  config.f1 = f1_list(freq_idx);
  config.DDC_freq = 250e6;
  config.zeropimods = [0 180];
  
  config.wfs(1).presums = 48;
  config.wfs(1).Tpd = 2e-6;
  config.wfs(1).tukey = 0.1;
  
  config.wfs(1).name = 'single';
  
  config.tx_weights = final_DDS_amp{cal_settings_list(freq_idx)};
  config.wfs(1).phase = final_DDS_phase{cal_settings_list(freq_idx)};
  config.phase = final_DDS_phase{cal_settings_list(freq_idx)};
  config.wfs(1).delay = final_DDS_time{cal_settings_list(freq_idx)};
  config.delay = final_DDS_time{cal_settings_list(freq_idx)};
  
  idx = find(strcmpi('AttenFirst18dB',{config.arena.ctu.out.bit_group.name}));
  config.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  config.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  
  idx = find(strcmpi('AttenSecond7dB',{config.arena.ctu.out.bit_group.name}));
  config.arena.ctu.out.bit_group(idx).epri = {[1 1]};
  config.arena.ctu.out.bit_group(idx).pri = {[1 1]};
  
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6,config.tg.Hice_thick);
  config.arena.fn = fullfile(arena_base_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  
  % Survey Test Mode (begins recording at 0 ft altitude)
  old_config = config;
  config.tg.Haltitude = config.tg.altitude_guard; % Start time 0 sec
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick_TEST', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6,config.tg.Hice_thick);
  config.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  config = old_config;
  
  % Loopback Mode (0e-6 delay line)
  old_config = config;
  config.tg.Haltitude = 3e8/2 * 0e-6; % Set the delay line here
  config.tg.altitude_guard = 3e8/2 * 0e-6;
  config.tg.Hice_thick = 800;
  config.tg.staged_recording = {[0],[0]};
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fusDelay_%.0fus_LOOPBACK', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude/(3e8/2)*1e6,config.wfs(end).Tpd*1e6);
  config.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  config = old_config;
  
  % Deconvolution Mode (for over calm lake or sea ice lead)
  old_config = config;
  idx = find(strcmpi('AttenFirst18dB',{config.arena.ctu.out.bit_group.name}));
  config.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  config.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  idx = find(strcmpi('AttenSecond7dB',{config.arena.ctu.out.bit_group.name}));
  config.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  config.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  config.tg.staged_recording = {[0],[0]};
  config.tg.altitude_guard = 3000*12*2.54/100;
  config.tg.Haltitude = 4000*12*2.54/100;
  config.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6);
  config.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  config = old_config;
  
  % Noise Mode
  old_config = config;
  config.tx_weights = [0];
  config.tx_enable = [0];
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_NOISE', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6);
  config.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  config = old_config;
end


%% Survey Mode + survey-test, loopback, noise, and deconv modes
% <2500 m thick ice, 1400 +/- 300 ft AGL
ice_thickness = [2500];
for freq_idx = 1
  param = default_radar_params_2018_Antarctica_TObas_accum;
  config = param.config;
  
  config.flight_hours = 4;
  config.prf = 1/50e-6;
  
  config.tg.staged_recording = {[1],[2]};
  config.tg.altitude_guard = 300*12*2.54/100;
  config.tg.Haltitude = 1400*12*2.54/100;
  config.tg.Hice_thick = ice_thickness(freq_idx);
  
  config.create_IQ = false;
  
  config.f0 = f0_list(freq_idx);
  config.f1 = f1_list(freq_idx);
  config.DDC_freq = 250e6;
  config.zeropimods = [0 180];
  
  config.wfs(1).presums = 60;
  config.wfs(1).Tpd = 2e-6;
  config.wfs(1).tukey = 0.1;
  
  config.wfs(1).name = 'single';
  
  config.tx_weights = final_DDS_amp{cal_settings_list(freq_idx)};
  config.wfs(1).phase = final_DDS_phase{cal_settings_list(freq_idx)};
  config.phase = final_DDS_phase{cal_settings_list(freq_idx)};
  config.wfs(1).delay = final_DDS_time{cal_settings_list(freq_idx)};
  config.delay = final_DDS_time{cal_settings_list(freq_idx)};
  
  idx = find(strcmpi('AttenFirst18dB',{config.arena.ctu.out.bit_group.name}));
  config.arena.ctu.out.bit_group(idx).epri = {[0 0]};
  config.arena.ctu.out.bit_group(idx).pri = {[0 0]};
  
  idx = find(strcmpi('AttenSecond7dB',{config.arena.ctu.out.bit_group.name}));
  config.arena.ctu.out.bit_group(idx).epri = {[1 1]};
  config.arena.ctu.out.bit_group(idx).pri = {[1 1]};
  
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6,config.tg.Hice_thick);
  config.arena.fn = fullfile(arena_base_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  
  % Survey Test Mode (begins recording at 0 ft altitude)
  old_config = config;
  config.tg.Haltitude = config.tg.altitude_guard; % Start time 0 sec
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick_TEST', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6,config.tg.Hice_thick);
  config.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  config = old_config;
end
