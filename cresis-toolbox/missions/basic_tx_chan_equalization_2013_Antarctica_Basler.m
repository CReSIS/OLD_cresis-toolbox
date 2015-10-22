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
% 5. Set ref_chan to correspond to an element in the center of the array
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

% .ref_chan = Reference transmit channel (surface location determined from this
%   channel and all phase measurements made relative to it)
param.ref_chan = 4;
% .rlines = Range lines to process
%   These are range lines post presumming
param.rlines = [1 inf];
% .rbins = Range bins to search for surface in
param.rbins = [5750 6500];
% .noise_rbins = Range bins to use for noise power calculation (THIS OFTEN NEEDS TO BE SET)
param.noise_rbins = 7000:8000;

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR for every transmit waveform are included in the estimate,
%   if even just one waveform does not meet the threshold the
%   range line is not used)
param.snr_threshold = 10;

param.radar_name = 'mcords4'; % param.radar_name: THIS OFTEN NEEDS TO BE SET
if strcmpi(param.radar_name,'mcords4')
  param.base_path = '/mnt/backup-iu/array1/20131216/mcords4/';
  param.base_path = '/basler/fp1/mcords4/';
  xml_fn = '';
%   xml_fn = '/basler/data/20131213/mcords4/mcords4_20131213_061516_00.xml';
  out_xml_fn_dir = '/basler/scratch2/';

  % .adc = the receive channel to use (relative to the board # so that it
  %    is always contained in [1,4] since there are 4 channels per board)
  param.adc = 1;
end

% param.wf_mapping (each entry indicates the waveform that should be used
%   for that tx channel... number of channels )
param.wf_mapping = [1 3 5 7 9 11 13 15];
param.bad_chan_mask = [0 0 0 0 0 0 0 0];

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = '/cresis/projects/dev/cr1/gps/2013_Antarctica_Basler/gps_20130921.mat';
param.gps_fn = '';
utc_time_correction = 3;

% .presums = Number of presums (coherent averaging) to do
param.presums = 20;

% .noise_removal_en = mcords noise removal (i.e. should be false for
%    mcrds, mcords2, etc)
param.noise_removal_en = false;

% Hwindow_desired = The desired window function for the transmit amplitude
%   THIS OFTEN NEEDS TO BE SET
window_func = @boxcar;
Hwindow_desired = [0.5 1 1 1 1 1 1 0.5];
% Hwindow_desired = [1 1 1 1 1 1 1];
max_DDS_amp = 50000;

time_delay_desired = [0 0 0 0 0 0 0 0];

phase_desired = [0 0 0 0 0 0 0 0];

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


% First, get the last XML file made
if isempty(xml_fn)
  xml_fns = get_filenames(param.base_path,'','','.xml');
  for xml_idx = 1:length(xml_fns)
    xml_fn = xml_fns{xml_idx};
    fprintf('=================== xml index %d ================\n', xml_idx);
    fprintf('%s\n', xml_fn);
    [settings,settings_enc] = read_cresis_xml(xml_fn);
    fprintf(' # of waveforms: %d\n', settings.DDS_Setup.Wave);
    settings.DDS_Setup.Waveforms(1)
  end
  xml_idx = input('Input xml index: ');
  xml_fn = xml_fns{xml_idx};
end

[settings,settings_enc] = read_cresis_xml(xml_fn);

xml_fns = get_filenames(param.base_path,'','','.xml');
for xml_idx = 1:length(xml_fns)
  xml_fn = xml_fns{xml_idx};
  [tmp xml_fn_name] = fileparts(xml_fn);
  time_stamp_idx = find(xml_fn_name == '_',1) + 1;
  year = str2double(xml_fn_name(time_stamp_idx + (0:3)));
  month = str2double(xml_fn_name(time_stamp_idx + 4 + (0:1)));
  day = str2double(xml_fn_name(time_stamp_idx + 6 + (0:1)));
  hour = str2double(xml_fn_name(time_stamp_idx + 9 + (0:1)));
  min = str2double(xml_fn_name(time_stamp_idx + 11 + (0:1)));
  sec = str2double(xml_fn_name(time_stamp_idx + 13 + (0:1)));
  %fprintf('  %s: %d\n', new_settings.fn, length(new_settings.Configuration.Waveforms));
  xml_date(xml_idx) = datenum(year,month,day,hour,min,sec);
end

xml_idx = find(strcmp(settings.fn, xml_fns));
date_begin = settings.datenum;
if xml_idx == length(xml_fns)
  date_end = inf;
else
  date_end = xml_date(xml_idx+1);
end

% Find the files that were created after this XML file
param.board = adc_to_board('mcords4',param.adc);
file_prefix = sprintf('mcords4_%02d_',param.board);
fns = get_filenames(fullfile(param.base_path,sprintf('chan%d',param.board)), 'mcords4', '', '.bin');
fn_dir = fullfile(param.base_path, sprintf('chan%d',param.board));
if isempty(fns)
  fprintf('  Could not find any files which match\n');
  return;
end

fn_mask = logical(size(fns));
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  
  fname = fname_info_mcords2(fn);
  fn_date(fn_idx) = fname.datenum;
  
  if fn_date(fn_idx) >= date_begin && fn_date(fn_idx) <= date_end
    fn_mask(fn_idx) = 1;
  end
end

fns = fns(fn_mask);
for fn_idx = 1:length(fns)
  fprintf('%d: %s\n', fn_idx, fns{fn_idx});
end
file_idxs = input('Input file indexes: ');
fns = fns(file_idxs);

% .DDS_start_mag = current DDS waveform attenuation in dB or DDS counts
%   (0 to 65535 DDS counts correspond to 0 to 2 Volts, linear map)
%
% Aux DAC IS ON FF for all
param.DDS_start_mag = double(settings.DDS_Setup.Ram_Amplitude);

% param.DDS_start_mag_units: Options are "dB" (NI/ledford) and "DDS"
param.DDS_start_mag_units = 'DDS';

% .DDS_start_time = current DDS start time in nanoseconds
param.DDS_start_time = settings.DDS_Setup.Waveforms(1).Delay/1e9;

% .DDS_start_phase = current DDS start phase in deg or DDS counts
%   (0 to 65535 DDS counts correspond to 0 to 2*pi, linear map)
param.DDS_start_phase = 180/pi*angle(exp(j*settings.DDS_Setup.Waveforms(1).Phase_Offset/180*pi));

% param.DDS_start_phase_units: Options are "deg" (NI) and "DDS" (ledford)
%   THIS OFTEN NEEDS TO BE SET
param.DDS_start_phase_units = 'deg';

% param.pc_param: pulse compression parameters
param.pc_param.f0 = settings.DDS_Setup.Waveforms(1).Start_Freq(1);
param.pc_param.f1 = settings.DDS_Setup.Waveforms(1).Stop_Freq(1);
param.pc_param.Tpd = settings.DDS_Setup.Base_Len * double(settings.DDS_Setup.Waveforms(1).Len_Mult);
param.pc_param.tukey = settings.DDS_Setup.RAM_Taper;

if ~isempty(param.gps_fn)
  fprintf('Loading GPS (%.1f sec)\n', toc(tstart));
  gps = load(param.gps_fn);
end

results = [];
for file_idx = 1:length(fns)
  fn = fns{file_idx};
  
  fprintf('Loading %s\n', fn);
  
  clear data;
  [hdr,data] = basic_load_mcords4(fn,struct('clk',fs/4));
  data_out = {};
  for chan = 1:length(param.wf_mapping)
    wf = param.wf_mapping(chan);
    data_out{chan} = data{wf}(:,:) - j*data{wf+1}(:,:);
  end
  data = data_out;
  
  [fn_dir fn_name] = fileparts(fn);
  
  % =======================================================================
  % Load GPS
  % =======================================================================
  if ~isempty(param.gps_fn)
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
    for chan = 1:length(param.wf_mapping)
      % Noise Removal
      data_pow = abs(data{chan}).^2;
      cfar_threshold = medfilt2(data_pow,[5 3]);
      cfar_threshold(:,1:3) = repmat(cfar_threshold(:,5),[1 3]);
      cfar_threshold(:,end-2:end) = repmat(cfar_threshold(:,end-4),[1 3]);
      cfar_threshold(1:5,:) = repmat(cfar_threshold(7,:),[5 1]);
      cfar_threshold(end-4:end,:) = repmat(cfar_threshold(end-6,:),[5 1]);
      
      data{chan}(data_pow > cfar_threshold*1000) = 0;
      
      % For debugging:
      %imagesc(lp(filter2(ones(5,21),data{wf})))
    end
  end
  
  % =======================================================================
  % Presumming/coherent averaging
  % =======================================================================
  if param.presums > 1
    fprintf('Coherent averaging (%.1f sec)\n', toc(tstart));
    clear data_out;
    for chan = 1:length(param.wf_mapping)
      for adc_idx = 1:size(data{chan},3)
        data_out{chan}(:,:,adc_idx) = fir_dec(data{chan}(:,:,adc_idx),param.presums);
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
  time = [];
  for chan = 1:length(param.wf_mapping)
    wf = param.wf_mapping(chan);
    pc_param.f0 = param.pc_param.f0;
    pc_param.f1 = param.pc_param.f1;
    pc_param.Tpd = param.pc_param.Tpd;
    pc_param.tukey = param.pc_param.tukey;
    pc_param.time = hdr.wfs(wf).t0 + (0:size(data{chan},1)-1)/fs;
    [data{chan},time{chan}] = pulse_compress(data{chan},pc_param);
  end
  
  chan = 1;
  if param.rbins(2) > size(data{chan},1)
    param.rbins(2) = size(data{chan},1);
  end
  if param.rlines(2) > size(data{chan},2)
    param.rlines(2) = size(data{chan},2);
  end
  rbins = param.rbins(1):param.rbins(2);
  rlines = param.rlines(1):param.rlines(2);
  
  % =======================================================================
  % Echogram plots
  % =======================================================================
  if param.plot_en
    for chan = 1:length(data)
      figure(chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('E %d',chan));
      imagesc(lp(data{chan}));
      title(sprintf('Chan %d File %s\nTime-Space Relative Power', chan, fn_name),'Interpreter','none');
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
  % surf_data = filter2(ones(1,10),abs(data{param.ref_chan}.^2));
  surf_data = filter2(ones(1,1),abs(data{param.ref_chan}.^2));
  [surf_vals surf_bins] = max(surf_data(rbins,rlines));
  surf_bins = rbins(1)-1 + surf_bins;
  
  if param.plot_en
    for chan = 1:length(data)
      figure(chan);
      hold on;
      plot(rlines, surf_bins,'k');
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
  
  % =======================================================================
  % Noise power estimate and SNR threshold
  % =======================================================================
  noise_power = mean(mean(abs(data{param.wf_mapping(param.ref_chan)}(param.noise_rbins,rlines)).^2));
  for chan = 1:length(param.wf_mapping)
    wf = param.wf_mapping(chan);
    for rline_idx = 1:length(rlines)
      rline = rlines(rline_idx);
      tx_phases(chan,rline_idx) = data{chan}(surf_bins(rline_idx),rline);
      tx_powers(chan,rline_idx) = abs(data{chan}(surf_bins(rline_idx),rline)).^2;
      
      [corr_out,lags] = xcorr(data{chan}(surf_bins(rline_idx)+search_bins,rline), ...
        data{param.ref_chan}(surf_bins(rline_idx)+ref_bins,rline) .* Hcorr_wind);
      corr_int = interpft(corr_out,param.Mt*length(corr_out));
      [peak_val(chan,rline_idx) peak_offset(chan,rline_idx)] = max(corr_int);
      peak_offset(chan,rline_idx) = (peak_offset(chan,rline_idx)-1)/param.Mt+1 ...
        + ref_bins(1) + search_bins(1) - 1 - zero_padding_offset;
    end
  end
  tx_snr = tx_powers ./ noise_power;
  good_meas = lp(tx_snr) > param.snr_threshold;
  good_rlines = zeros(size(rlines));
  good_rlines(sum(good_meas(~param.bad_chan_mask,:)) == sum(~param.bad_chan_mask)) = 1;
  good_rlines = logical(good_rlines);
  
  num_good_rlines = sum(good_rlines);
  fprintf('Number of good range lines: %d out of %d\n', num_good_rlines, length(good_rlines));
  fprintf('========================================================\n');
  
  
  % =======================================================================
  % Time Offset Settings
  % =======================================================================
  clear ref_time_median;
  peak_offset_time = peak_offset * (time{1}(2)-time{1}(1));
  for chan = 1:length(param.wf_mapping)
    wf = param.wf_mapping(chan);
    ref_time = peak_offset_time(chan,:);
    if param.plot_en
      figure(120+chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Time %d',chan));
      plot(ref_time);
      %     ylim([min(peak_offset_time(:)) max(peak_offset_time(:))]);
    end
    
    %ref_time_median(chan) = median(peak_offset(good_rlines));
    ref_time_median(chan) = mean(ref_time(good_rlines));
    fprintf('TX %d: relative time: %10.4f ns\n', chan, ...
      1e9*median(ref_time(good_rlines)));
    fprintf('    Std. dev. time: %10.4f ns\n', ...
      1e9*std(ref_time(good_rlines)));
    ref_time(~good_rlines) = NaN;
    if param.plot_en
      hold on;
      plot(ref_time,'ro');
      hold off;
      title(sprintf('Relative Time (%d to ref %d)', chan, param.ref_chan));
    end
  end
  fprintf('Recommended new DDS time offset (ns):\n');
  new_DDS_time = param.DDS_start_time - ref_time_median;
  fprintf('%.4f\t', new_DDS_time(1:end-1)*1e9);
  fprintf('%.4f', new_DDS_time(end)*1e9);
  fprintf('\n');
  fprintf('========================================================\n');
  results.DDS_time_error(file_idx,:) = ref_time_median;
  results.DDS_time(file_idx,:) = new_DDS_time;
  
  % =======================================================================
  % Amplitude Settings
  % =======================================================================
  fprintf('Relative power for each waveform (dB)\n');
  clear delta_power;
  for chan = 1:length(param.wf_mapping)
    wf = param.wf_mapping(chan);
    ref_power = tx_powers(chan,:)./tx_powers(param.ref_chan,:);
    if param.plot_en
      figure(10+chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Pow %d',chan));
      plot(lp(ref_power),1);
    end
    delta_power(chan) = lp(median(ref_power(good_rlines)),1);
    fprintf('%10.2f\n', delta_power(chan));
    %   fprintf('WF %d: relative power: %10.2f dB\n', wf, lp(median(ref_power(good_rlines))));
    %   fprintf('    Std. dev. power: %.2f dB\n', lp(std(ref_power(good_rlines))));
    ref_power(~good_rlines) = NaN;
    if param.plot_en
      hold on;
      plot(lp(ref_power),'ro');
      hold off;
      title(sprintf('Relative Power (%d to ref %d)', chan, param.ref_chan));
    end
  end
  if strcmpi(param.DDS_start_mag_units,'DDS')
    fprintf('%s WINDOW: Recommended new DDS amplitude settings (DDS counts):\n', ...
      upper(func2str(window_func)));
    new_DDS_amp = param.DDS_start_mag./(10.^(delta_power/20) ./ Hwindow_desired);
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
  results.DDS_amp_error(file_idx,:) = delta_power;
  results.DDS_amp(file_idx,:) = new_DDS_amp;

  % =======================================================================
  % Phase Settings
  % =======================================================================
  clear ref_phase_median;
  for chan = 1:length(param.wf_mapping)
    wf = param.wf_mapping(chan);
    ref_phase = angle(tx_phases(chan,:)./tx_phases(param.ref_chan,:));
    if param.plot_en
      figure(20+chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Ang %d',chan));
      plot(ref_phase);
      ylim([-pi pi]);
    end
    ref_phase_median(chan) = angle(median(exp(j*ref_phase(good_rlines))));
    fprintf('WF %d: relative phase: %10.4f rad, %10.1f deg\n', chan, ...
      ref_phase_median(chan), ref_phase_median(chan)*180/pi);
    ref_phase(~good_rlines) = NaN;
    if param.plot_en
      hold on;
      plot(ref_phase,'ro');
      hold off;
      title(sprintf('Relative Phase (%d to ref %d)', chan, param.ref_chan));
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
    new_DDS_phase = 180/pi*angle(exp(j*(param.DDS_start_phase - DDS_error/65536*360)/180*pi));
    fprintf('%.1f\t', new_DDS_phase(1:end-1));
    fprintf('%.1f', new_DDS_phase(end));
    fprintf('\n');
  end
  fprintf('========================================================\n');
  results.DDS_phase_error(file_idx,:) = DDS_error/65536*360;
  results.DDS_phase(file_idx,:) = new_DDS_phase;
    
end

%% Print DDS Time
fprintf('DDS_time_error (ns)\n');
for file_idx = 1:length(fns)
  fprintf('%s', fns{file_idx});
  for wf = 1:size(results.DDS_time_error,2)
    fprintf('\t%.2f', results.DDS_time_error(file_idx,wf)*1e9);
  end
  fprintf('\n');
end
fprintf('DDS_time (ns)\n');
for file_idx = 1:length(fns)
  fprintf('%s', fns{file_idx});
  for wf = 1:size(results.DDS_time,2)
    fprintf('\t%.2f', results.DDS_time(file_idx,wf)*1e9);
  end
  fprintf('\n');
end
fprintf('Mean');
final_DDS_time = mean(results.DDS_time,1)*1e9;
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.2f', final_DDS_time(wf));
end
fprintf('\n');
fprintf('Median');
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.2f', median(results.DDS_time(:,wf))*1e9);
end
fprintf('\n');
fprintf('Stdev');
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.2f', std(results.DDS_time(:,wf))*1e9);
end
fprintf('\n');

%% Print DDS Amplitude
fprintf('DDS_amp_error (dB)\n');
for file_idx = 1:length(fns)
  fprintf('%s', fns{file_idx});
  for wf = 1:size(results.DDS_amp_error,2)
    fprintf('\t%.1f', results.DDS_amp_error(file_idx,wf));
  end
  fprintf('\n');
end
fprintf('DDS_amp (DDS Counts, linear)\n');
for file_idx = 1:length(fns)
  fprintf('%s', fns{file_idx});
  for wf = 1:size(results.DDS_amp,2)
    fprintf('\t%.0f', results.DDS_amp(file_idx,wf));
  end
  fprintf('\n');
end
fprintf('Mean');
final_DDS_amp = mean(results.DDS_amp,1);
final_DDS_amp = final_DDS_amp / max(final_DDS_amp(~param.bad_chan_mask)) * max_DDS_amp;
for wf = 1:size(results.DDS_amp,2)
  fprintf('\t%.0f', final_DDS_amp(wf));
end
fprintf('\n');
fprintf('Median');
for wf = 1:size(results.DDS_amp,2)
  fprintf('\t%.0f', median(results.DDS_amp(:,wf)));
end
fprintf('\n');
fprintf('Stdev');
for wf = 1:size(results.DDS_amp,2)
  fprintf('\t%.0f', std(results.DDS_amp(:,wf)));
end
fprintf('\n');

%% Print DDS Phase
fprintf('DDS_phase_error (deg)\n');
for file_idx = 1:length(fns)
  fprintf('%s', fns{file_idx});
  for wf = 1:size(results.DDS_phase_error,2)
    fprintf('\t%.1f', results.DDS_phase_error(file_idx,wf));
  end
  fprintf('\n');
end
fprintf('DDS_phase (deg)\n');
for file_idx = 1:length(fns)
  fprintf('%s', fns{file_idx});
  for wf = 1:size(results.DDS_phase,2)
    fprintf('\t%.0f', results.DDS_phase(file_idx,wf));
  end
  fprintf('\n');
end
fprintf('Mean');
final_DDS_phase = angle(mean(exp(j*results.DDS_phase/180*pi),1))*180/pi;
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.1f', final_DDS_phase(wf));
end
fprintf('\n');
fprintf('Median');
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.1f', angle(median(exp(j*results.DDS_phase(:,wf)/180*pi)))*180/pi);
end
fprintf('\n');
fprintf('Stdev');
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.1f', std(results.DDS_phase(:,wf)));
end
fprintf('\n');

settings_enc.DDSa20Setup.Ra61ma20Amplitude(~param.bad_chan_mask) = uint16(final_DDS_amp(~param.bad_chan_mask));
for wf = 1:length(settings_enc.DDSa20Setup.Wa61veforms)
  if mod(wf,2)
    settings_enc.DDSa20Setup.Wa61veforms(wf).Pha61sea20Offset(~param.bad_chan_mask) ...
      = double(final_DDS_phase(~param.bad_chan_mask) + phase_desired(~param.bad_chan_mask));
  else
    settings_enc.DDSa20Setup.Wa61veforms(wf).Pha61sea20Offset(~param.bad_chan_mask) = ...
      angle(exp(j*(double(final_DDS_phase(~param.bad_chan_mask) + phase_desired(~param.bad_chan_mask))/180*pi+pi/2)))*180/pi;
  end
  settings_enc.DDSa20Setup.Wa61veforms(wf).Dela61y(~param.bad_chan_mask) ...
    = double(final_DDS_time(~param.bad_chan_mask) + time_delay_desired(~param.bad_chan_mask));
end

[xml_fn_dir xml_fn_name xml_fn_ext] = fileparts(settings.fn);
out_xml_fn = fullfile(out_xml_fn_dir, sprintf('%s_txequal%s', xml_fn_name, xml_fn_ext));

settings_enc = rmfield(settings_enc,'fn');
settings_enc = rmfield(settings_enc,'datenum');
fprintf('Writing %s\n', out_xml_fn);
fid = fopen(out_xml_fn,'w');
fprintf(fid,'<?xml version=''1.0'' standalone=''yes'' ?>\n');
fprintf(fid,'<LVData xmlns="http://www.ni.com/LVData">\n');
write_ni_xml_object(settings_enc,fid,true);
fprintf(fid,'</LVData>');
fclose(fid);

return;








