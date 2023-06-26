% script create_configs_2018_Antarctica_Ground_rds
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
f0_list = [170e6];
f1_list = [230e6];
cal_settings_list = [1];
prf = 1/100e-6;

%% Calculate the desired presums
% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
physical_constants; % c = speed of light
velocity = [4]; % m/s
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
    final_DDS_phase{idx} = [0 0];
    final_DDS_phase_no_time{idx} = [0 0]; % not used usually
    final_DDS_amp{idx} = [1 1];
    final_DDS_time{idx} =  [0 0];
  end
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % 170-230 MHz
  idx = 1;
  final_DDS_phase{idx} = [0 0];
  final_DDS_phase_no_time{idx} = [0 0]; % not used usually
  final_DDS_amp{idx} = [1 1];
  final_DDS_time{idx} =  [0 0];
end

%% Survey Mode + loopback, noise, and deconv modes
% <1700 m thick ice, 1200 +/- 500 ft AGL
ice_thickness = [4500];
for freq_idx = [1]
  param = default_radar_params_2018_Antarctica_GroundLab_rds;
  config = param.config;
  
  config.flight_hours = 10;
  config.prf = prf;
  
  config.zeropimods = [0 180]; % 180 causes a special mode in some ADC/DAC
  
  config.tg.staged_recording = {[1 2]};
  config.tg.altitude_guard = 100*12*2.54/100;
  config.tg.Haltitude = 100*12*2.54/100;
  config.tg.Hice_thick = ice_thickness(freq_idx);
  
  config.create_IQ = false;
  
  config.f0 = f0_list(freq_idx);
  config.f1 = f1_list(freq_idx);
  config.DDC_freq = 160e6;
  config.zeropimods = [0 180];
  
  config.wfs(1).presums = 256;
  config.wfs(2).presums = 256;
  config.tukey = 0.08;
  config.wfs(1).Tpd = 3e-6;
  config.wfs(2).Tpd = 10e-6;
  config.wfs(1).tukey = 0.1;
  config.wfs(2).tukey = 0.1;
  
  config.wfs(1).name = '3us_high';
  config.wfs(2).name = '10us_high';
  
  config.tx_weights = final_DDS_amp{cal_settings_list(freq_idx)};
  config.wfs(1).phase = final_DDS_phase{cal_settings_list(freq_idx)};
  config.wfs(2).phase = final_DDS_phase{cal_settings_list(freq_idx)};
  config.phase = final_DDS_phase{cal_settings_list(freq_idx)};
  config.wfs(1).delay = final_DDS_time{cal_settings_list(freq_idx)};
  config.wfs(2).delay = final_DDS_time{cal_settings_list(freq_idx)};
  config.delay = final_DDS_time{cal_settings_list(freq_idx)};
  
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6,config.tg.Hice_thick);
  config.arena.fn = fullfile(arena_base_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  
  for tr_antenna = 1:2
    % Loopback Mode (0e-6 delay line)
    old_config = config;
    config.tg.Haltitude = 3e8/2 * 0e-6; % Set the delay line here
    config.tg.altitude_guard = 3e8/2 * 0e-6;
    config.tg.Hice_thick = 2000;
    config.tg.staged_recording = {[0],[0]};
    config.tx_weights = [0 0];
    config.tx_enable = [0 0];
    config.tx_weights(tr_antenna) = 1;
    config.tx_enable(tr_antenna) = 1;
    config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fusDelay_%.0fus_LOOPBACK_TRANTENNA%d', ...
      config.f0/1e6,config.f1/1e6,config.tg.Haltitude/(3e8/2)*1e6,config.wfs(end).Tpd*1e6, tr_antenna);
    config.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',config.arena.psc_name));
    write_radar_config(config);
    config = old_config;
  end
  
  % Noise Mode
  old_config = config;
  config.tx_weights = [0 0];
  config.tx_enable = [0 0];
  config.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_NOISE', ...
    config.f0/1e6,config.f1/1e6,config.tg.Haltitude*100/12/2.54,config.wfs(end).Tpd*1e6);
  config.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',config.arena.psc_name));
  write_radar_config(config);
  config = old_config;
end
