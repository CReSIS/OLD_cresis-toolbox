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
  calval_dir = 'cal_val';
end

f0_list = [600e6];
f1_list = [900e6];
DDC_select_list = [1]; % Which DDC mode to use
cal_settings = [1];
prf = 20000;

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
  final_DDS_amp{idx} = [0.2];
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
  [~,defaults] = default_radar_params_2018_Antarctica_TObas_accum3;
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
  
  param.wfs(1).presums = 128;
  param.wfs(1).Tpd = 2e-6;
  param.wfs(1).tukey = 0.1;
  
  param.wfs(1).name = 'single';
  
  param.tx_weights = final_DDS_amp{cal_settings(freq_idx)};
  param.wfs(1).phase = final_DDS_phase{cal_settings(freq_idx)};
  param.phase = final_DDS_phase{cal_settings(freq_idx)};
  param.wfs(1).delay = final_DDS_time{cal_settings(freq_idx)};
  param.delay = final_DDS_time{cal_settings(freq_idx)};
  
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
  
  % Loopback Mode (10e-6 delay line)
  old_param = param;
  param.tg.Haltitude = 3e8/2 * 10e-6; % Set the delay line here
  param.tg.altitude_guard = 3e8/2 * 3e-6;
  param.tg.Hice_thick = 0;
  param.tg.staged_recording = false;
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
  param.tg.staged_recording = false;
  param.tg.altitude_guard = 3000*12*2.54/100;
  param.tg.Haltitude = 4000*12*2.54/100;
  param.tg.Hice_thick = 0 * 12*2.54/100/sqrt(er_ice);
  param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_DECONV', ...
    param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6);
  param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
  write_radar_config(param);
  param = old_param;
  
  if freq_idx == 1
    % Noise Mode
    old_param = param;
    param.tx_weights = [0];
    param.tx_enable = [0];
    param.tg.staged_recording = false;
    param.arena.psc_name = sprintf('survey_%.0f-%.0fMHz_%.0fft_%.0fus_NOISE', ...
      param.f0/1e6,param.f1/1e6,param.tg.Haltitude*100/12/2.54,param.wfs(end).Tpd*1e6);
    param.arena.fn = fullfile(arena_base_dir,calval_dir,sprintf('%s.xml',param.arena.psc_name));
    write_radar_config(param);
    param = old_param;
  end
end
