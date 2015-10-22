% script basic_tx_chan_equalization
%
% This script is for helping with setting the DDS start phase values
% and amplitude values while on the plane.  It requires loading
% a single receive channel file and then analyzing it with this
% program.
%
% Useful Tip: Search for "NEEDS" and that will help you set all the fields
% that you will need to change.
%
% Steps for N-element antenna array
% 1. Collect data with N+1 waveforms, waveform 1 should transmit from
%    antenna 1, waveform 2 should transmit from antenna 2, etc. with
%    the last waveform being all transmitters together.
%    All other properties should be the same for each waveform
% 2. Only one receiver is required by this routine
% 3. Characterization should be done for each pulse duration since
%    there appears to be some variation with pulse duration
% 4. Time gate should be large enough to include noise-only data which
%    will be used to determine the SNR of the surface return (SNR
%    threshold is used to exclude low SNR data points from the measurement)
% 5. Set ref_wf to correspond to an element in the center of the array
%    because correlation statistics are better for shorter baselines
%
% Sometimes two runs are required:
% 1. One run with xlim,ylim,caxis at their defaults
% 2. A second run with rlines, rbins (and xlim, ylim, caxis) set to
%    a specific region
%
% Author: John Paden

physical_constants;
close all;
tstart = tic;

% =======================================================================
% User Settings
% =======================================================================

fs = 1e9/9;

% .plot_en = flag to enable plots
param.plot_en = true;

% .caxis = Color axis limits (leave empty first time since this causes it to use the
% defaults).
param.caxis = [];
%param.caxis = [50 120];

% .ylim = Leave empty the first time (it just uses the default limits then)
param.ylim = [];
%param.ylim = [1 500];

% .xlim = Leave empty the first time (it just uses the default limits then)
%   These limits are in range lines post presumming
param.xlim = [];
%param.xlim = [200 450];

% .ref_wf = Reference transmit channel (surface location determined from this
%   channel and all phase measurements made relative to it)
param.ref_wf = 3;
% .rlines = Range lines to process
%   These are range lines post presumming
param.rlines = [1 inf];
% .rbins = Range bins to search for surface in
param.rbins = [1 inf];
% .noise_rbins = Range bins to use for noise power calculation (THIS OFTEN NEEDS TO BE SET)
param.noise_rbins = 2000:3000;

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR for every transmit waveform are included in the estimate,
%   if even just one waveform does not meet the threshold the
%   range line is not used)
param.snr_threshold = 10;

% .DDS_start_mag = current DDS waveform attenuation in dB or DDS counts
%   (0 to 65535 DDS counts correspond to 0 to 2 Volts, linear map)
%
% DO NOT EXCEED 35000
% Aux DAC IS ON 80 for all
param.DDS_start_mag = [30000 30000 30000 30000 30000 30000 30000];
param.DDS_start_mag = [14904	18361	20073	30000	19609	16481	18362];
param.DDS_start_mag = [11338	14499	16929	30000	15997	12999	14028];

% param.DDS_start_mag_units: Options are "dB" (NI/ledford) and "DDS"
%   THIS OFTEN NEEDS TO BE SET
param.DDS_start_mag_units = 'DDS';

% .DDS_start_phase = current DDS start phase in deg or DDS counts
%   (0 to 65535 DDS counts correspond to 0 to 2*pi, linear map)
%   (THIS NEEDS TO BE SET EVERYTIME)
param.DDS_start_phase = [0 0 0 0 0 0 0];
param.DDS_start_phase = [29.1	66.8	0	149.4	258.9	340.2	17.9];

% param.DDS_start_phase_units: Options are "deg" (NI) and "DDS" (ledford)
%   THIS OFTEN NEEDS TO BE SET
param.DDS_start_phase_units = 'deg';

param.radar_name = 'mcords3'; % param.radar_name: THIS OFTEN NEEDS TO BE SET
if strcmpi(param.radar_name,'mcords3')
  % .board = board number from 0 to 3
  param.board = 0;
  % .adc = the receive channel to use (relative to the board # so that it
  %    is always contained in [1,4] since there are 4 channels per board)
  param.adc = 3;
  param.acquisition_num = 0;
  
  % Parameters to locate specific file of interest
  %  (THIS NEEDS TO BE SET EVERYTIME)
  param.seg = '';
  param.file_num = 1;
  param.base_path = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_P3/20131108/';
end

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = '/cresis/projects/dev/cr1/gps/2013_Antarctica_P3/gps_20131108.mat';

% .presums = Number of presums (coherent averaging) to do
param.presums = 10;

% .noise_removal_en = mcords noise removal (i.e. should be false for
%    mcrds, mcords2, etc)
param.noise_removal_en = false;

% param.pc_param: pulse compression parameters
%   THIS OFTEN NEEDS TO BE SET
param.pc_param.f0 = 180e6;
param.pc_param.f1 = 210e6;
param.pc_param.Tpd = 10e-6;
param.pc_param.tukey = 0.2;

% Hwindow_desired = The desired window function for the transmit amplitude
%   THIS OFTEN NEEDS TO BE SET
window_func = @boxcar;
Hwindow_desired = window_func(param.DDS_start_mag).';
Hwindow_desired = [1 1 1 1 1 1 1];
Hwindow_desired = [1 1 1 1 1 1 1];

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================
num_wf = length(param.DDS_start_phase);

% =======================================================================
% Load data
% =======================================================================
fprintf('========================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(tstart));

if strcmpi(param.radar_name,'mcords3')
  if isempty(param.seg)
    fn_dir = fullfile(param.base_path, sprintf('board%d',param.board));
  else
    fn_dir = fullfile(param.base_path, sprintf('board%d',param.board), ...
      param.seg);
  end
  file_prefix = sprintf('mcords3_%d_',param.board);
  if isempty(param.acquisition_num)
    file_suffix = sprintf('%04d.bin',param.file_num);
  else
    file_suffix = sprintf('%02d_%04d.bin',param.acquisition_num,param.file_num);
  end
  fprintf('  Path: %s\n', fn_dir);
  fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
  fn = get_filename(fn_dir, file_prefix, '', file_suffix);
  if isempty(fn)
    fprintf('  Could not find any files which match\n');
    return;
  end
  fprintf('  Loading file %s\n', fn);
  [hdr,data] = basic_load_mcords3(fn,struct('clk',fs));
  for wf = 1:length(data)
    data{wf} = data{wf}(:,:,param.adc);
  end
end

[fn_dir fn_name] = fileparts(fn);

if length(data) < num_wf
  fprintf('Number of waveforms required %d does not match data %d\n', num_wf, length(data));
  return;
end
wf = 1;

% =======================================================================
% Load GPS
% =======================================================================
if ~isempty(param.gps_fn)
  fprintf('Loading GPS (%.1f sec)\n', toc(tstart));
  gps = load(param.gps_fn);
  
  [year month day] = datevec(epoch_to_datenum(gps.gps_time(1)));
  hdr.gps_time = datenum_to_epoch(datenum(year,month,day,0,0,hdr.utc_time_sod));
  hdr.gps_time = hdr.gps_time + utc_leap_seconds(hdr.gps_time(1));
  
  roll = interp1(gps.gps_time, gps.roll, hdr.gps_time);
  roll = fir_dec(roll,param.presums);
  
  figure(100); clf;
  plot(roll*180/pi);
  grid on;
  ylabel('Roll (deg)');
  if ~isempty(param.xlim)
    xlim(param.xlim);
  end
end

% =======================================================================
% Time domain burst noise (digital errors) removal
% =======================================================================
if param.noise_removal_en
  fprintf('Noise/digital error removal (%.1f sec)\n', toc(tstart));
  for wf = 1:num_wf
    % Noise Removal
    data{wf} = data{wf} - median(data{wf}(:,1));
    data_pow = abs(data{wf}).^2;
    cfar_threshold = medfilt2(data_pow,[5 3]);
    cfar_threshold(:,1:3) = repmat(cfar_threshold(:,5),[1 3]);
    cfar_threshold(:,end-2:end) = repmat(cfar_threshold(:,end-4),[1 3]);
    cfar_threshold(1:5,:) = repmat(cfar_threshold(7,:),[5 1]);
    cfar_threshold(end-4:end,:) = repmat(cfar_threshold(end-6,:),[5 1]);

    data{wf}(data_pow > cfar_threshold*1000) = 0;

    % For debugging:
    %imagesc(lp(filter2(ones(5,21),data{wf})))
  end
end

% =======================================================================
% Presumming/coherent averaging
% =======================================================================
if param.presums > 1
  fprintf('Coherent averaging (%.1f sec)\n', toc(tstart));
  for wf = 1:num_wf
    for adc_idx = 1:size(data{wf},3)
      data_out{wf}(:,:,adc_idx) = fir_dec(data{wf}(:,:,adc_idx),param.presums);
    end
  end
  data = data_out;
  clear data_out;
end

% =======================================================================
% Pulse compression
% =======================================================================
fprintf('Pulse compression (%.1f sec)\n', toc(tstart));
clear pc_param time;
for wf = 1:num_wf
  pc_param(wf).f0 = param.pc_param.f0;
  pc_param(wf).f1 = param.pc_param.f1;
  pc_param(wf).Tpd = param.pc_param.Tpd;
  pc_param(wf).tukey = param.pc_param.tukey;
  pc_param(wf).time = hdr.wfs(wf).t0 + (0:size(data{wf},1)-1)/fs;
  [data{wf},time{wf}] = pulse_compress(data{wf},pc_param(wf));
end

if param.rbins(2) > size(data{wf},1)
  param.rbins(2) = size(data{wf},1);
end
if param.rlines(2) > size(data{wf},2)
  param.rlines(2) = size(data{wf},2);
end
param.rbins = param.rbins(1):param.rbins(2);
param.rlines = param.rlines(1):param.rlines(2);

% =======================================================================
% Echogram plots
% =======================================================================
if param.plot_en
  for wf = 1:num_wf
    figure(wf); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('E %d',wf));
    imagesc(lp(data{wf}));
    title(sprintf('Wf %d File %s\nTime-Space Relative Power', wf, fn_name),'Interpreter','none');
    grid on;
    colorbar
    if ~isempty(param.caxis)
      caxis(param.caxis);
    end
    if ~isempty(param.ylim)
      ylim(param.ylim);
    end
    if ~isempty(param.xlim)
      xlim(param.xlim);
    end
  end
end

% =======================================================================
% Surface tracker
% =======================================================================
surf_data = filter2(ones(1,10),abs(data{param.ref_wf}.^2));
[surf_vals surf_bins] = max(surf_data(param.rbins,param.rlines));
surf_bins = param.rbins(1)-1 + surf_bins;

if param.plot_en
  for wf = 1:num_wf
    figure(wf);
    hold on;
    plot(param.rlines, surf_bins,'k');
    hold off;
  end
end

% =======================================================================
% Noise power estimate and SNR threshold
% =======================================================================
noise_power = mean(mean(abs(data{param.ref_wf}(param.noise_rbins,param.rlines)).^2));
clear tx_phases tx_powers;
for wf = 1:num_wf
  for rline_idx = 1:length(param.rlines)
    rline = param.rlines(rline_idx);
    tx_phases(wf,rline_idx) = data{wf}(surf_bins(rline_idx),rline);
    tx_powers(wf,rline_idx) = abs(data{wf}(surf_bins(rline_idx),rline)).^2;
  end
end
tx_snr = tx_powers ./ noise_power;
good_meas = lp(tx_snr) > param.snr_threshold;
good_rlines = zeros(size(param.rlines));
good_rlines(sum(good_meas) == size(tx_snr,1)) = 1;
good_rlines = logical(good_rlines);

num_good_rlines = sum(good_rlines);
fprintf('Number of good range lines: %d out of %d\n', num_good_rlines, length(good_rlines));
fprintf('========================================================\n');

% =======================================================================
% Amplitude Settings
% =======================================================================
fprintf('Relative power for each waveform (dB)\n');
clear delta_power;
for wf = 1:num_wf
  ref_power = tx_powers(wf,:)./tx_powers(param.ref_wf,:);
  if param.plot_en
    figure(10+wf); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Pow %d',wf));
    plot(lp(ref_power),1);
  end
  delta_power(wf) = lp(median(ref_power(good_rlines)),1);
  fprintf('%10.2f\n', delta_power(wf));
  %   fprintf('WF %d: relative power: %10.2f dB\n', wf, lp(median(ref_power(good_rlines))));
  %   fprintf('    Std. dev. power: %.2f dB\n', lp(std(ref_power(good_rlines))));
  ref_power(~good_rlines) = NaN;
  if param.plot_en
    hold on;
    plot(lp(ref_power),'ro');
    hold off;
    title(sprintf('Relative Power (%d to ref %d)', wf, param.ref_wf));
  end
end
if strcmpi(param.DDS_start_mag_units,'DDS')
  fprintf('%s WINDOW: Recommended new DDS amplitude settings (DDS counts):\n', ...
    upper(func2str(window_func)));
  new_DDS_amp = param.DDS_start_mag./(10.^(delta_power/20) .* Hwindow_desired);
  fprintf('%.0f\t', new_DDS_amp(1:end-1));
  fprintf('%.0f', new_DDS_amp(end));
  fprintf('\n');
elseif strcmpi(param.DDS_start_mag_units,'dB')
  fprintf('%s WINDOW: Recommended new DDS amplitude settings (dB):\n', ...
    upper(func2str(window_func)));
  new_DDS_amp = param.DDS_start_mag + 20*log10(10.^(delta_power/20) ./ Hwindow_desired);
  fprintf('%.2f\t', new_DDS_amp(1:end-1));
  fprintf('%.2f', new_DDS_amp(end));
  fprintf('\n');
end
fprintf('========================================================\n');

% =======================================================================
% Phase Settings
% =======================================================================
clear ref_phase_median;
for wf = 1:num_wf
  ref_phase = angle(tx_phases(wf,:)./tx_phases(param.ref_wf,:));
  if param.plot_en
    figure(20+wf); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Ang %d',wf));
    plot(ref_phase);
    ylim([-pi pi]);
  end
  ref_phase_median(wf) = median(ref_phase(good_rlines));
  fprintf('WF %d: relative phase: %10.4f rad, %10.1f deg\n', wf, ...
    median(ref_phase(good_rlines)), median(ref_phase(good_rlines))*180/pi);
  fprintf('    Std. dev. phase: %.4f rad, %.1f deg\n', ...
    std(ref_phase(good_rlines)), std(ref_phase(good_rlines))*180/pi);
  ref_phase(~good_rlines) = NaN;
  if param.plot_en
    hold on;
    plot(ref_phase,'ro');
    hold off;
    title(sprintf('Relative Phase (%d to ref %d)', wf, param.ref_wf));
  end
end
DDS_error = ref_phase_median * 65536/(2*pi);
if strcmpi(param.DDS_start_phase_units,'DDS')
  fprintf('Recommended new DDS phase settings (DDS counts):\n');
  new_DDS_phase = mod(param.DDS_start_phase - DDS_error, 65536);
  fprintf('%.0f\t', new_DDS_phase(1:end-1));
  fprintf('%.0f', new_DDS_phase(end));
  fprintf('\n');
elseif strcmpi(param.DDS_start_phase_units,'deg')
  fprintf('Recommended new DDS phase settings (deg):\n');
  new_DDS_phase = mod(param.DDS_start_phase - DDS_error/65536*360, 360);
  fprintf('%.1f\t', new_DDS_phase(1:end-1));
  fprintf('%.1f', new_DDS_phase(end));
  fprintf('\n');
end
fprintf('========================================================\n');

