% script create_settings_2018_Antarctica_Ground
%
% Creates Arena accumulation radar settings
%
% See link_budgets.xls

physical_constants; % c = speed of light

% Define waveforms
if ispc
  arena_base_dir = 'C:\waveforms_arena\';
  calval_dir = 'cal_val';
else
  arena_base_dir = '/arena/waveforms/';
  calval_dir = 'cal_val';
end

f0_list = [160e6];
f1_list = [230e6];
DDC_select_list = [1]; % Which DDC mode to use
cal_settings = [1];
prf = 12500;

% presums: Ensure lambda/4 sampling (fudge factor allows difference) and an
%   that presums are an even number.
velocity = [4]; % m/s
presums = round(c./max(abs(f0_list),abs(f1_list)) ./ velocity * prf / 4)

final_DDS_phase = [];
final_DDS_phase_no_time = [];
final_DDS_amp = [];
final_DDS_time = [];
if 0
  % Initial conditions (usually all zeros phase/time with max amplitude)
  for idx = 1:length(f0_list)
    final_DDS_phase{idx} = [0 0 0 0];
    final_DDS_phase_no_time{idx} = [0 0 0 0]; % not used usually
    final_DDS_amp{idx} = [1 1 1 1];
    final_DDS_time{idx} =  [0 0 0 0];
  end
else
  % COPY AND PASTE RESULTS FROM basic_tx_chan_equalization_SEASON_NAME.m
  % HERE:
  
  % 160-230 MHz
  idx = 1;
  final_DDS_phase{idx} = [0 0 0 0];
  final_DDS_phase_no_time{idx} = [0 0 0 0]; % not used usually
  final_DDS_amp{idx} = [1 1 1 1];
  final_DDS_time{idx} =  [0 0 0 0];
end

% Hwindow_orig: Desired window created during transmit calibration
%  This is used any time a window that is different from that used
%  during calibration is to be used.
Hwindow_orig = chebwin(4).';

%% SETUP
% =========================================================================

param = [];
param.season_name = '2018_Antarctica_Ground';
param.radar_name = 'mcords6';
param.gps_source = 'arena-field';
clear phase_centers;
num_tx = 4;
for tx_chan = 1:num_tx
  % Just enable the current antenna to determine its phase center
  tx_weights = zeros(1,num_tx);
  tx_weights(tx_chan) = 1;
  rxchan = 1; % Fix the receiver (it should not matter which one you choose)
  % Determine phase center for the antenna
  phase_centers(:,tx_chan) = lever_arm(param, tx_weights, rxchan);
end
% Adjust phase centers to the mean phase center position
phase_centers = bsxfun(@minus,phase_centers,mean(phase_centers,2));

%% Survey Mode + loopback, noise, and deconv modes
% <1700 m thick ice, 1200 +/- 500 ft AGL
ice_thickness = [4500];
for freq_idx = [1]
  [~,defaults] = default_radar_params_2018_Antarctica_Ground_mcords6;
  param = defaults{1};
  param.flight_hours = 4;
  param.prf = prf;
  
  param.zeropimods = [0 180]; % 180 causes a special mode in some ADC/DAC
  
  param.tg.staged_recording = [1 2];
  param.tg.altitude_guard = 0*12*2.54/100;
  param.tg.Haltitude = 0*12*2.54/100;
  param.tg.Hice_thick = ice_thickness(freq_idx);
  
  param.create_IQ = false;
  
  param.f0 = f0_list(freq_idx);
  param.f1 = f1_list(freq_idx);
  param.DDC_freq = 160e6;
  param.zeropimods = [0 180];
  
  param.wfs(1).presums = 64;
  param.wfs(2).presums = 128;
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
  
  idx = find(strcmpi('Atten',{param.arena.ctu.out.bit_group.name}));
  param.arena.ctu.out.bit_group(idx).epri = {[20 20],[0 0]};
  param.arena.ctu.out.bit_group(idx).pri = {[20 20],[0 0]};
  
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_%.0fmthick', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6,param.tg.Hice_thick);
  param.arena.fn = fullfile(arena_base_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);

  % Loopback Mode (0e-6 delay line)
  param.tg.Haltitude = 3e8/2 * 0e-6;
  param.tg.staged_recording = false;
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fusDelay_%.0fus_LOOPBACK', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude/(3e8/2)*1e6,param.wfs(end).Tpd*1e6);
  param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  
  if freq_idx == 1
    % Noise Mode
    param.tg.staged_recording = [1 2];
    param.tx_weights = [0 0 0 0];
    param.tx_eable = [0 0 0 0];
    param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_NOISE', ...
      param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6);
    param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
    write_radar_config(param);
  end
end
