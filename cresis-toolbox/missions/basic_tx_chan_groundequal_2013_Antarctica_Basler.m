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
% 5. Set ref_tx to correspond to an element in the center of the array
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

fs = 1e9/2;

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

% .ref_tx = Reference transmit channel (surface location determined from this
%   channel and all phase measurements made relative to it)
param.ref_tx = 1;
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
% Aux DAC IS ON FF for all
param.DDS_start_mag = [32500 32500 32500 32500 32500 32500 32500 32500];

% param.DDS_start_mag_units: Options are "dB" (NI/ledford) and "DDS"
%   THIS OFTEN NEEDS TO BE SET
param.DDS_start_mag_units = 'DDS';

% .DDS_start_phase = current DDS start phase in deg or DDS counts
%   (0 to 65535 DDS counts correspond to 0 to 2*pi, linear map)
%   (THIS NEEDS TO BE SET EVERYTIME)
param.DDS_start_phase = [0 0 0 0 0 0 0 0];

param.DDS_start_time = [0 0 0 0 0 0 0 0];

% param.DDS_start_phase_units: Options are "deg" (NI) and "DDS" (ledford)
%   THIS OFTEN NEEDS TO BE SET
param.DDS_start_phase_units = 'deg';

param.adc = 1;

param.radar_name = 'mcords4'; % param.radar_name: THIS OFTEN NEEDS TO BE SET

fn_dir = '/basler/data/mcords4/chan1/';
fns = {};
% //raw data, without compensation
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_231947_00_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_232155_01_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_232334_02_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_232543_03_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_232751_04_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_232926_05_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_233211_06_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131209_233340_07_0000.bin');

% //after delay compensation
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_012148_00_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_012329_01_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_012507_02_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_013442_03_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_013645_04_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_013831_05_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_014426_06_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_014600_07_0000.bin');

% //after phase and delay compensation
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_022625_08_0001.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_022807_09_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_022939_10_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_023933_11_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_024119_12_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_024313_13_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_024527_14_0000.bin');
% fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_024709_15_0000.bin');

% after phase, amp, delay compensation
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_035739_00_0000.bin');
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_035916_01_0000.bin');
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_040055_02_0000.bin');
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_040311_03_0000.bin');
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_040454_04_0000.bin');
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_040629_05_0000.bin');
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_040826_06_0000.bin');
fns{end+1} = fullfile(fn_dir,'mcords4_01_20131210_041004_07_0000.bin');

wf = 1;

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = '';
utc_time_correction = 0;

% .presums = Number of presums (coherent averaging) to do
param.presums = 4;

% .noise_removal_en = mcords noise removal (i.e. should be false for
%    mcrds, mcords2, etc)
param.noise_removal_en = false;

% param.pc_param: pulse compression parameters
%   THIS OFTEN NEEDS TO BE SET
param.pc_param.f0 = 200e6;
param.pc_param.f1 = 450e6;
param.pc_param.Tpd = 10e-6;
param.pc_param.tukey = 0.2;

% Hwindow_desired = The desired window function for the transmit amplitude
%   THIS OFTEN NEEDS TO BE SET
window_func = @boxcar;
Hwindow_desired = window_func(length(param.DDS_start_mag)).';
Hwindow_desired = [1 1 1 1 1 1 1 1];

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

% =======================================================================
% Load data
% =======================================================================
fprintf('========================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(tstart));

if strncmpi(param.radar_name,'mcords4',length('mcords4'))
  data={};
  for fn_idx=1:length(fns)
    fn = fns{fn_idx};
    fprintf('  Loading file %s\n', fn);
    [hdr,tmp_data] = basic_load_mcords4(fn,struct('clk',fs/4));
    data{fn_idx} = tmp_data{wf} -j*tmp_data{wf+1};
  end
end

[fn_dir fn_name] = fileparts(fn);

% =======================================================================
% Load GPS
% =======================================================================
if ~isempty(param.gps_fn)
  fprintf('Loading GPS (%.1f sec)\n', toc(tstart));
  gps = load(param.gps_fn);
  finfo = fname_info_mcords2(fn);
  [year,month,day] = datevec(finfo.datenum);
  hdr.gps_time = datenum_to_epoch(datenum(year,month,day,0,0,hdr.utc_time_sod)); % Still UTC time
  hdr.gps_time = hdr.gps_time + utc_leap_seconds(hdr.gps_time(1)) + utc_time_correction; % Convert from UTC to GPS
  
  roll = interp1(gps.gps_time, gps.roll, hdr.gps_time);
  roll = fir_dec(roll,param.presums);
  
  elev = interp1(gps.gps_time, gps.elev, hdr.gps_time);
  elev = fir_dec(elev,param.presums);
  
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
  for tx_idx = 1:length(data)
    % Noise Removal
    data{tx_idx} = data{tx_idx} - median(data{tx_idx}(:,1));
    data_pow = abs(data{tx_idx}).^2;
    cfar_threshold = medfilt2(data_pow,[5 3]);
    cfar_threshold(:,1:3) = repmat(cfar_threshold(:,5),[1 3]);
    cfar_threshold(:,end-2:end) = repmat(cfar_threshold(:,end-4),[1 3]);
    cfar_threshold(1:5,:) = repmat(cfar_threshold(7,:),[5 1]);
    cfar_threshold(end-4:end,:) = repmat(cfar_threshold(end-6,:),[5 1]);

    data{tx_idx}(data_pow > cfar_threshold*1000) = 0;

    % For debugging:
    %imagesc(lp(filter2(ones(5,21),data{wf})))
  end
end

% =======================================================================
% Presumming/coherent averaging
% =======================================================================
if param.presums > 1
  fprintf('Coherent averaging (%.1f sec)\n', toc(tstart));
  for tx_idx = 1:length(data)
    data_out{tx_idx} = fir_dec(data{tx_idx},param.presums);
  end
  data = data_out;
  clear data_out;
end

% =======================================================================
% Pulse compression
% =======================================================================
fprintf('Pulse compression (%.1f sec)\n', toc(tstart));
clear pc_param time;
for tx_idx = 1:length(data)
  pc_param(tx_idx).f0 = param.pc_param.f0;
  pc_param(tx_idx).f1 = param.pc_param.f1;
  pc_param(tx_idx).Tpd = param.pc_param.Tpd;
  pc_param(tx_idx).tukey = param.pc_param.tukey;
  pc_param(tx_idx).decimate = true;
  pc_param(tx_idx).time = hdr.wfs(wf).t0 + (0:size(data{tx_idx},1)-1)/fs;
  [data{tx_idx},time{tx_idx}] = pulse_compress(data{tx_idx},pc_param(tx_idx));
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
  for tx_idx = 1:length(data)
    figure(tx_idx); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('E %d',tx_idx));
    imagesc(lp(data{tx_idx}));
    title(sprintf('Tx %d File %s\nTime-Space Relative Power', tx_idx, fn_name),'Interpreter','none');
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
% surf_data = filter2(ones(1,10),abs(data{param.ref_tx}.^2));
surf_data = filter2(ones(1,1),abs(data{param.ref_tx}.^2));
[surf_vals surf_bins] = max(surf_data(param.rbins,param.rlines));
surf_bins = param.rbins(1)-1 + surf_bins;

if param.plot_en
  for tx_idx = 1:length(data)
    figure(tx_idx);
    hold on;
    plot(param.rlines, surf_bins,'k');
    hold off;
  end
end

% =======================================================================
% Find peaks values
% =======================================================================
param.ref_bins = [-12 12];
param.search_bins = [-7 7];
param.Mt = 100;
ref_bins = param.ref_bins(1):param.ref_bins(2);
search_bins = param.search_bins(1)+param.ref_bins(1) : param.search_bins(2)+param.ref_bins(2);
zero_padding_offset = length(search_bins) - length(ref_bins);
Hcorr_wind = hanning(length(ref_bins));
clear tx_phases tx_powers peak_val peak_offset;

for tx_idx = 1:length(data)
  for rline_idx = 1:length(param.rlines)
    rline = param.rlines(rline_idx);
    tx_phases(tx_idx,rline_idx) = data{tx_idx}(surf_bins(rline_idx),rline);
    tx_powers(tx_idx,rline_idx) = abs(data{tx_idx}(surf_bins(rline_idx),rline)).^2;
    
    [corr_out,lags] = xcorr(data{tx_idx}(surf_bins(rline_idx)+search_bins,rline), ...
      data{param.ref_tx}(surf_bins(rline_idx)+ref_bins,rline) .* Hcorr_wind);
    corr_int = interpft(corr_out,param.Mt*length(corr_out));
    [peak_val(tx_idx,rline_idx) peak_offset(tx_idx,rline_idx)] = max(corr_int);
    peak_offset(tx_idx,rline_idx) = (peak_offset(tx_idx,rline_idx)-1)/param.Mt+1 ...
      + ref_bins(1) + search_bins(1) - 1 - zero_padding_offset;
  end
end

% All range lines are good
good_rlines = logical(ones(1,size(data{1},2)));

% =======================================================================
% Time Offset Settings
% =======================================================================
clear ref_time_median;
peak_offset_time = peak_offset * (time{1}(2)-time{1}(1));
for tx_idx = 1:length(data)
  ref_time = peak_offset_time(tx_idx,:);
  if param.plot_en
    figure(120+tx_idx); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Time %d',tx_idx));
    plot(ref_time);
%     ylim([min(peak_offset_time(:)) max(peak_offset_time(:))]);
  end

  %ref_time_median(tx_idx) = median(peak_offset(good_rlines));
  ref_time_median(tx_idx) = mean(ref_time(good_rlines));
  fprintf('TX %d: relative time: %10.4f ns\n', tx_idx, ...
    1e9*median(ref_time(good_rlines)));
  fprintf('    Std. dev. time: %10.4f ns\n', ...
    1e9*std(ref_time(good_rlines)));
  ref_time(~good_rlines) = NaN;
  if param.plot_en
    hold on;
    plot(ref_time,'ro');
    hold off;
    title(sprintf('Relative Time (%d to ref %d)', tx_idx, param.ref_tx));
  end
end
fprintf('Recommended new DDS time offset (ns):\n');
new_DDS_time = param.DDS_start_time - ref_time_median;
fprintf('%.4f\t', new_DDS_time(1:end-1)*1e9);
fprintf('%.4f', new_DDS_time(end)*1e9);
fprintf('\n');
fprintf('========================================================\n');

% =======================================================================
% Amplitude Settings
% =======================================================================
fprintf('Relative power for each waveform (dB)\n');
clear delta_power;
for tx_idx = 1:length(data)
  ref_power = tx_powers(tx_idx,:)./tx_powers(param.ref_tx,:);
  if param.plot_en
    figure(10+tx_idx); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Pow %d',tx_idx));
    plot(lp(ref_power),1);
  end
  %delta_power(tx_idx) = lp(median(ref_power(good_rlines)),1);
  delta_power(tx_idx) = lp(mean(ref_power(good_rlines)),1);
  fprintf('%10.2f\n', delta_power(tx_idx));
  %   fprintf('WF %d: relative power: %10.2f dB\n', wf, lp(median(ref_power(good_rlines))));
  %   fprintf('    Std. dev. power: %.2f dB\n', lp(std(ref_power(good_rlines))));
  ref_power(~good_rlines) = NaN;
  if param.plot_en
    hold on;
    plot(lp(ref_power),'ro');
    hold off;
    title(sprintf('Relative Power (%d to ref %d)', tx_idx, param.ref_tx));
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
for tx_idx = 1:length(data)
  ref_phase = angle(tx_phases(tx_idx,:)./tx_phases(param.ref_tx,:));
  if param.plot_en
    figure(20+tx_idx); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Ang %d',tx_idx));
    plot(ref_phase);
    ylim([-pi pi]);
  end
  %ref_phase_median(tx_idx) = median(ref_phase(good_rlines));
  ref_phase_median(tx_idx) = mean(ref_phase(good_rlines));
  fprintf('WF %d: relative phase: %10.4f rad, %10.1f deg\n', tx_idx, ...
    median(ref_phase(good_rlines)), median(ref_phase(good_rlines))*180/pi);
  fprintf('    Std. dev. phase: %.4f rad, %.1f deg\n', ...
    std(ref_phase(good_rlines)), std(ref_phase(good_rlines))*180/pi);
  ref_phase(~good_rlines) = NaN;
  if param.plot_en
    hold on;
    plot(ref_phase,'ro');
    hold off;
    title(sprintf('Relative Phase (%d to ref %d)', tx_idx, param.ref_tx));
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





