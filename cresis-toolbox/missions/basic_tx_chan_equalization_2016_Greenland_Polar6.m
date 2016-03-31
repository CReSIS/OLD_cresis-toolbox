% script basic_tx_chan_equalization
%
% This script is for helping with setting the DDS delay, phase,
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
%    the last waveform being all transmitters together (this script
%    does not make use of this last waveform).
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
param = [];

% =======================================================================
% User Settings
% =======================================================================

% fs = sampling frequency (Hz)
fs = 1600e6;

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

% .rlines = 1x2 vector specifying range lines to process (to select all
%   range lines to the end, set second element to inf). These are range
%   lines post presumming
param.rlines = [1 inf];
% .rbins = 1x2 vector specifying range bins to search for surface in (to
%   select range bins to the end, set second element to inf).These are
%   post pulse compression which complex basebands and decimates the data.
%   (THIS OFTEN NEEDS TO BE SET)
param.rbins = [1 inf];
param.rbins = [1000 1300]; % 150-520 MHz 3 us
% param.rbins = [50 150]; % 180-210 MHz 3 us
% param.rbins = [300 500]; % 150-520 MHz 1 us
% param.rbins = [20 60]; % 180-210 MHz 1 us
% .noise_rbins = 1xN vector specifying range bins to use for noise power
%   estimate on each range line (THIS OFTEN NEEDS TO BE SET)
param.noise_rbins = 400:600; % 150-520 MHz 3 us
% param.noise_rbins = 1:40; % 180-210 MHz 3 us
% param.noise_rbins = 100:200; % 150-520 MHz 1 us
% param.noise_rbins = 80:100; % 180-210 MHz 1 us

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR for every transmit waveform are included in the estimate,
%   if even just one waveform does not meet the threshold the
%   range line is not used)
param.snr_threshold = 10;

% param.radar_name: THIS OFTEN NEEDS TO BE SET
param.radar_name = 'mcords5';

% param.base_dir = (THIS OFTEN NEEDS TO BE SET)
param.base_dir = '\\172.18.1.103\d\';
param.base_dir = '\\172.18.1.32\mcords\';
param.base_dir = '\\172.18.1.33\mcords\';
% param.base_dir = 'D:\';
param.base_dir = 'E:\awi\awi\';

% xml_fn = optional string containing xml file to load (normally
%   leave this empty and select during run time)
xml_fn = '';

% xml_version = See cresis_xml_mapping and wiki page on National
%   Instruments XML files to determine which version to use
xml_version = 2.0;

% out_xml_fn_dir = String containg the directory where the new XML file
%   will be placed (THIS OFTEN NEEDS TO BE SET)
out_xml_fn_dir = 'E:\waveforms\';

% param.adc = The ADC channel to load. Generally not critical which one is used.
param.adc = 8;

% param.wf_mapping: each entry indicates the waveform that should be used
%   for that DDS. An entry set to zero indicates that that DDS will not
%   be used.
% param.wf_mapping = [-1j -3j -5j -7j -9j -11j -13j -15j];
param.wf_mapping = [1 2 3 4 5 6 7 8];

% param.bad_chan_mask: logical array equal in size to number of DDS
%   Any DDS with '1' will not be used in calculations or its settings affected
param.bad_chan_mask = [0 0 0 0 0 0 0 0];

% .ref_chan = Reference transmit channel (surface location determined from this
%   channel and all measurements made relative to it)
param.ref_chan = 4;

% .gps_fn = Optional GPS file name (leave empty to disable), not used by
%   the program
param.gps_fn = '';

% utc_time_correction = scalar containing number of seconds to add to the
%   time stored in the raw data files (if radar time was correct, this
%   value would be 0)
utc_time_correction = 1;

% param.presums = Interger containing number of presums (coherent
%   averaging) to do
param.presums = 10;

% param.noise_removal_en = mcords noise removal (i.e. should be false for
%    all other systems)
param.noise_removal_en = false;

% Desired delay, phase, and amplitudes after equalization (THIS SECTION
%   OFTEN NEEDS TO BE SET)
% Hwindow_desired = The desired window function for the transmit
%   amplitudes. Units are relative.
Hwindow_desired = chebwin(8,30).';
% max_DDS_amp = Amplitudes will be normalized to this scalar value
% max_DDS_amp = 15000;
max_DDS_amp = [4000 4000 4000 4000 4000 4000 4000 4000];
% time_delay_desired = 1xNc vector of desired delays. Units are
%   nanoseconds. A positive value means that this channel will be
%   transmitted later relative to the other channels.
time_delay_desired = [0 0 0 0 0 0 0 0];
% phase_desired = 1xNc vector of desired phases. Units are
%   degrees.
phase_desired = [0 0 0 0 0 0 0 0];
% update_delay: logical which causes the DDS delay values to be updated
%   and the phase to be delay-compensated
update_delay = false; % START WITH THIS ONE
% update_amplitude: logical which causes the DDS RAM (amplitude) values to be updated
update_amplitude = true; % AFTER UPDATING DELAY, ENABLE AMP/PHASE UPDATE and DISABLE DELAY
% update_phase: logical which causes the DDS phase values to be updated
update_phase = true; % AFTER UPDATING DELAY, ENABLE AMP/PHASE UPDATE and DISABLE DELAY

% param.ref_bins = 1x2 vector specifying bins relative to peak to use in
%   correlation
param.ref_bins = [-12 12];
% param.search_bins = 1x2 vector specifying max search range to use when
%   looking for the best correlation (this is to ensure that each output
%   correlation value has full support... i.e. no roll off effect)
param.search_bins = [-15 15];
% param.Mt = amount to oversample the correlation
param.Mt = 100;

% =======================================================================
% =======================================================================
%% Automated Section
% =======================================================================
% =======================================================================

cresis_xml_mapping;

% =======================================================================
% Load data
% =======================================================================
fprintf('========================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(tstart));


%% User selects XML file to use
if isempty(xml_fn)
  [settings,settings_enc] = read_ni_xml_directory(param.base_dir,xml_file_prefix,false);
  for xml_idx = 1:length(settings)
    xml_fn = settings(xml_idx).fn;
    fprintf('=================== xml index %d ================\n', xml_idx);
    fprintf('%s\n', xml_fn);
    fprintf(' # of waveforms: %d\n', length(settings(xml_idx).(config_var).Waveforms));
    settings(xml_idx).(config_var).Ram_Amplitude
    settings(xml_idx).(config_var).Waveforms(1)
  end
  xml_idx = input('Input xml index: ');
  xml_fn = settings(xml_idx).fn;
end

%% Determine time range that this XML file corresponds to
time_begin = settings(xml_idx).datenum;
if xml_idx == length(settings)
  time_end = inf;
else
  time_end = settings(xml_idx+1).datenum;
end
settings = settings(xml_idx);
settings_enc = settings_enc(xml_idx);

%% Allow operation with XML files with no Delay variable
for wf = 1:length(settings.(config_var).Waveforms)
  if ~isfield(settings.(config_var).Waveforms(wf), 'Delay')
    settings.(config_var).Waveforms(wf).Delay = zeros(size(settings.(config_var).Waveforms(wf).(phase_var)));
  end
end

%% Find the files that were created during the time range associated with
% this XML file
param.board = adc_to_board(param.radar_name,param.adc);
file_prefix = sprintf('%s_%02d_',param.radar_name,param.board);
if strcmpi(param.radar_name,'mcords3')
  fns = get_filenames(fullfile(param.base_dir,sprintf('board%d',param.board)), 'mcords3', '', '.bin');
  fn_dir = fullfile(param.base_dir, sprintf('board%d',param.board));
elseif strcmpi(param.radar_name,'mcords4')
  fns = get_filenames(fullfile(param.base_dir,sprintf('chan%d',param.board)), 'mcords4', '', '.bin');
  fn_dir = fullfile(param.base_dir, sprintf('chan%d',param.board));
elseif strcmpi(param.radar_name,'mcords5')
  fns = get_filenames(fullfile(param.base_dir,sprintf('chan%d',param.board)), 'mcords5', '', '.bin');
  fn_dir = fullfile(param.base_dir, sprintf('chan%d',param.board));
end
if isempty(fns)
  fprintf('  Could not find any files which match\n');
  return;
end

fn_mask = logical(zeros(size(fns)));
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  
  fname = fname_info_mcords2(fn);
  fn_date(fn_idx) = fname.datenum;
  
  if fn_date(fn_idx) >= time_begin && fn_date(fn_idx) <= time_end
    fn_mask(fn_idx) = 1;
  end
end

%% User selects which files to use in transmit equalization
fns = fns(fn_mask);
for fn_idx = 1:length(fns)
  fprintf('%d: %s\n', fn_idx, fns{fn_idx});
end
file_idxs = input('Input file indexes: ');
fns = fns(file_idxs);

%% Extract important parameters from the XML file
% .DDS_start_mag = current DDS waveform attenuation in dB or DDS counts
%   (0 to 65535 DDS counts correspond to 0 to 2 Volts, linear map)
param.DDS_start_mag = double(settings.(config_var).(ram_var));

% param.DDS_start_mag_units: Options are "dB" (NI/ledford) and "DDS"
param.DDS_start_mag_units = 'DDS';

% .DDS_start_time = current DDS start time in nanoseconds
param.DDS_start_time = settings.(config_var).Waveforms(1).Delay/1e9;

% .DDS_start_phase = current DDS start phase in deg or DDS counts
%   (0 to 65535 DDS counts correspond to 0 to 2*pi, linear map)
param.DDS_start_phase = 180/pi*angle(exp(j*settings.(config_var).Waveforms(1).(phase_var)/180*pi));

% param.DDS_start_phase_units: Options are "deg" (NI) and "DDS" (ledford)
%   THIS OFTEN NEEDS TO BE SET
param.DDS_start_phase_units = 'deg';

% param.pc_param: pulse compression parameters
if strcmpi(param.radar_name,'mcords5') && settings.DDC_Ctrl.DDC_sel.Val >= 1
  %% DDC Enabled
  param.pc_param.f0 = settings.(config_var).Waveforms(1).Start_Freq(1);
  param.pc_param.f1 = settings.(config_var).Waveforms(1).Stop_Freq(1);
else
  %% No DDC
  param.pc_param.f0 = settings.(config_var).Waveforms(1).Start_Freq(1);
  param.pc_param.f1 = settings.(config_var).Waveforms(1).Stop_Freq(1);
end
param.pc_param.Tpd = settings.(config_var).Base_Len * double(settings.(config_var).Waveforms(1).Len_Mult);
param.pc_param.tukey = settings.(config_var).RAM_Taper;

%% Iterate through each file and extract equalization information
results = [];
for file_idx = 1:length(fns)
  fn = fns{file_idx};
  
  fprintf('Loading %s\n', fn);
  
  clear data;
  if strcmp(param.radar_name,'mcords3')
    [hdr,data] = basic_load_mcords3(fn,struct('clk',fs));
  elseif strcmp(param.radar_name,'mcords4')
    [hdr,data] = basic_load_mcords4(fn,struct('clk',fs/4));
  elseif strcmp(param.radar_name,'mcords5')
    [hdr,data] = basic_load_mcords5(fn,struct('clk',fs));
  end
  data_out = cell(size(param.bad_chan_mask));
  for chan = 1:length(param.wf_mapping)
    wf = param.wf_mapping(chan);
    if wf ~= 0
      if strcmp(param.radar_name,'mcords3')
        if isreal(wf)
          data_out{chan} = data{wf}(:,:,mod(param.adc-1,4)+1);
        else
          data_out{chan} = data{abs(wf)}(:,:,param.adc) + sign(wf)*data{abs(wf)+1}(:,:,mod(param.adc-1,4)+1);
        end
      else
        if isreal(wf)
          data_out{chan} = data{wf}(:,:);
        else
          data_out{chan} = data{abs(wf)}(:,:) + sign(wf)*data{abs(wf)+1}(:,:);
        end
      end
    end
  end
  data = data_out;
  
  [fn_dir fn_name] = fileparts(fn);
  
  %% Optional loading of the GPS data (this is not used, but can be useful
  % for debugging
  if ~isempty(param.gps_fn)
    fprintf('Loading GPS (%.1f sec)\n', toc(tstart));
    gps = load(param.gps_fn);

    % A not very robust way to get the UTC date so that we can get the
    % absolute time that the radar data was collected (radars only record
    % seconds of day)
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
  %% Time domain burst noise (digital errors) removal
  % =======================================================================
  if param.noise_removal_en
    fprintf('Noise/digital error removal (%.1f sec)\n', toc(tstart));
    for chan = 1:length(param.wf_mapping)
      if param.wf_mapping(chan) ~= 0
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
  end
  
  % =======================================================================
  %% Presumming/coherent averaging
  % =======================================================================
  if param.presums > 1
    fprintf('Coherent averaging (%.1f sec)\n', toc(tstart));
    data_out = cell(size(param.bad_chan_mask));
    for chan = 1:length(param.wf_mapping)
      if param.wf_mapping(chan) ~= 0
        for adc_idx = 1:size(data{chan},3)
          data_out{chan}(:,:,adc_idx) = fir_dec(data{chan}(:,:,adc_idx),param.presums);
        end
      end
    end
    data = data_out;
    clear data_out;
  end
  
  % =======================================================================
  %% Pulse compression
  % =======================================================================
  fprintf('Pulse compression (%.1f sec)\n', toc(tstart));
  clear pc_param time;
  time = [];
  for chan = 1:length(param.wf_mapping)
    if param.wf_mapping(chan) ~= 0
      wf = abs(param.wf_mapping(chan));
      pc_param.decimate = true;
      pc_param.zero_pad = true;
      pc_param.window_func = @hanning;
      pc_param.f0 = param.pc_param.f0;
      pc_param.f1 = param.pc_param.f1;
      pc_param.Tpd = param.pc_param.Tpd;
      pc_param.tukey = param.pc_param.tukey;
      if settings.DDC_Ctrl.DDC_sel.Val
        pc_param.DDC_mode = true;
        pc_param.DDC_freq = double(settings.DDC_Ctrl.NCO_freq) * 1e6;
        pc_param.time = hdr.wfs(wf).t0 + (0:size(data{chan},1)-1)/fs ...
          * 2^(double(settings.DDC_Ctrl.DDC_sel.Val)+1);
      else
        pc_param.DDC_mode = false;
        pc_param.time = hdr.wfs(wf).t0 + (0:size(data{chan},1)-1)/fs;
      end
      [data{chan},time{chan}] = pulse_compress(data{chan},pc_param);
    end
  end

  %% Determine region of data to use for estimating equalization coefficients
  chan = find(param.wf_mapping,1); % Get first non-zero entry
  if param.rbins(2) > size(data{chan},1)
    param.rbins(2) = size(data{chan},1);
  end
  if param.rlines(2) > size(data{chan},2)
    param.rlines(2) = size(data{chan},2);
  end
  rbins = param.rbins(1):param.rbins(2);
  rlines = param.rlines(1):param.rlines(2);
  
  % =======================================================================
  %% Echogram plots
  % =======================================================================
  if param.plot_en
    for chan = 1:length(data)
      if param.wf_mapping(chan) ~= 0
        figure(chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('E %d',chan));
        imagesc(lp(data{chan}));
        h_echogram_axes(chan) = gca;
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
    linkaxes(h_echogram_axes,'xy');
  end
  
  % =======================================================================
  %% Surface tracker
  % =======================================================================
  % Incoherent along-track filtering
  surf_data = filter2(ones(1,6),abs(data{param.ref_chan}.^2));
  % Simple max search to find surface
  [surf_vals surf_bins] = max(surf_data(rbins,rlines));
  surf_bins = rbins(1)-1 + surf_bins;
  
  if param.plot_en
    for chan = 1:length(data)
      if param.wf_mapping(chan) ~= 0
        figure(chan);
        hold on;
        plot(rlines, surf_bins,'k');
        hold off;
      end
    end
  end
    
  % =======================================================================
  %% Noise power estimate and SNR threshold
  % =======================================================================
  noise_power = mean(mean(abs(data{abs(param.wf_mapping(param.ref_chan))}(param.noise_rbins,rlines)).^2));

  % =======================================================================
  %% Extract delay (using oversampled cross correlation), phase and amplitude
  % differences between channels
  % =======================================================================
  ref_bins = param.ref_bins(1):param.ref_bins(2);
  search_bins = param.search_bins(1)+param.ref_bins(1) : param.search_bins(2)+param.ref_bins(2);
  zero_padding_offset = length(search_bins) - length(ref_bins);
  Hcorr_wind = hanning(length(ref_bins));
  clear tx_phases tx_powers peak_val peak_offset;
  
  tx_powers = zeros(length(data),length(rlines));
  for chan = 1:length(param.wf_mapping)
    if param.wf_mapping(chan) ~= 0
      wf = abs(param.wf_mapping(chan));
      for rline_idx = 1:length(rlines)
        rline = rlines(rline_idx);
        % Get the phases and powers right at the peak of the reference channel
        % surface
        tx_phases(chan,rline_idx) = data{chan}(surf_bins(rline_idx),rline);
        tx_powers(chan,rline_idx) = abs(data{chan}(surf_bins(rline_idx),rline)).^2;

        % Gets the time offset relative to the reference channel (a postive
        % offset means that the channel leads the reference channel)
        [corr_out,lags] = xcorr(data{chan}(surf_bins(rline_idx)+search_bins,rline), ...
          data{param.ref_chan}(surf_bins(rline_idx)+ref_bins,rline) .* Hcorr_wind);
        corr_int = interpft(corr_out,param.Mt*length(corr_out));
        [peak_val(chan,rline_idx) peak_offset(chan,rline_idx)] = max(corr_int);
        peak_offset(chan,rline_idx) = (peak_offset(chan,rline_idx)-1)/param.Mt+1 ...
          + ref_bins(1) + search_bins(1) - 1 - zero_padding_offset;
      end
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
  if num_good_rlines == 0
    error('Cannot continue: no range lines exceeded the snr_threshold');
  end
    
  % =======================================================================
  %% Time Offset Settings
  % =======================================================================
  ref_time_mean = zeros(size(param.DDS_start_time));
  peak_offset_time = peak_offset * (time{1}(2)-time{1}(1));
  for chan = 1:length(param.wf_mapping)
    if param.wf_mapping(chan) ~= 0
      wf = abs(param.wf_mapping(chan));
      ref_time = peak_offset_time(chan,:);
      
      if chan == param.ref_chan
        median_mask = ones(size(good_rlines));
        ref_time_mean(chan) = 0;
      else
        % Remove outliers and take mean
        std_val = std(ref_time(good_rlines));
        median_val = median(ref_time(good_rlines));
        median_mask = ref_time >= median_val - std_val ...
          & ref_time <= median_val + std_val;
        ref_time_mean(chan) = mean(ref_time(good_rlines & median_mask));
      end
      
      fprintf('TX %d: %4.2f ns (%4.2f ns)\n', chan, 1e9*ref_time_mean(chan), ...
        1e9*std(ref_time(good_rlines & median_mask)));
      if param.plot_en
        figure(120+chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Time %d',chan));
        plot(ref_time);
      end
      ref_time(~(good_rlines & median_mask)) = NaN;
      if param.plot_en
        hold on;
        plot(ref_time,'ro');
        hold off;
        title(sprintf('Relative Time (%d to ref %d)', chan, param.ref_chan));
        ylim([min(min(ref_time),-1e-11) max(1e-11,max(ref_time))]);
      end
    end
  end
  if 0
    % Limit time delay precision to 1/20th of a range bin
    ref_time_mean = round(ref_time_mean / (time{1}(2)-time{1}(1)) * 20) ...
      * (time{1}(2)-time{1}(1)) / 20;
  end
  fprintf('Recommended new DDS time offset (ns):\n');
  new_DDS_time = param.DDS_start_time - ref_time_mean;
  fprintf('%.4f\t', new_DDS_time(1:end-1)*1e9);
  fprintf('%.4f', new_DDS_time(end)*1e9);
  fprintf('\n');
  fprintf('========================================================\n');
  results.DDS_time_error(file_idx,:) = ref_time_mean;
  results.DDS_time(file_idx,:) = new_DDS_time;
  
  % =======================================================================
  %% Amplitude Settings
  % =======================================================================
  fprintf('Relative power for each waveform (dB)\n');
  delta_power = zeros(size(param.DDS_start_mag));
  for chan = 1:length(param.wf_mapping)
    if param.wf_mapping(chan) ~= 0
      wf = abs(param.wf_mapping(chan));
      ref_power = tx_powers(chan,:)./tx_powers(param.ref_chan,:);
      
      if chan == param.ref_chan
        median_mask = ones(size(good_rlines));
        delta_power(chan) = 0;
      else
        % Remove outliers and take mean
        std_val = std(ref_power(good_rlines));
        mean_val = mean(ref_power(good_rlines));
        median_mask = ref_power >= mean_val - std_val ...
          & ref_power <= mean_val + std_val;
        delta_power(chan) = lp(mean(ref_power(good_rlines & median_mask)));
      end
      
      fprintf('TX %d: %4.2f dB (%4.2f dB), desired %4.2f\n', chan, delta_power(chan), ...
        lp(std(ref_power(good_rlines & median_mask))), 20*log10(Hwindow_desired(chan)));
      if param.plot_en
        figure(10+chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Pow %d',chan));
        plot(lp(ref_power,1));
      end
      ref_power(~(good_rlines & median_mask)) = NaN;
      if param.plot_en
        hold on;
        plot(lp(ref_power,1),'ro');
        hold off;
        title(sprintf('Relative Power (%d to ref %d)', chan, param.ref_chan));
      end
    end
  end
  if strcmpi(param.DDS_start_mag_units,'DDS')
    fprintf('Recommended new DDS amplitude settings (DDS counts):\n');
    new_DDS_amp = param.DDS_start_mag .* (Hwindow_desired ./ (10.^(delta_power/20)));
    fprintf('%.0f\t', new_DDS_amp(1:end-1));
    fprintf('%.0f', new_DDS_amp(end));
    fprintf('\n');
  elseif strcmpi(param.DDS_start_mag_units,'dB')
    fprintf('Recommended new DDS amplitude settings (dB):\n');
    new_DDS_amp = param.DDS_start_mag + 20*log10(  Hwindow_desired ./ (10.^(delta_power/20))  ) ;
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
  ref_phase_mean = zeros(size(param.DDS_start_phase));
  for chan = 1:length(param.wf_mapping)
    if param.wf_mapping(chan) ~= 0
      wf = abs(param.wf_mapping(chan));
      ref_phase = tx_phases(chan,:).*conj(tx_phases(param.ref_chan,:));
      ref_values_real = real(ref_phase);
      ref_values_imag = imag(ref_phase);
      
      median_mask = ones(size(good_rlines));
      ref_phase_mean(chan) = mean(ref_phase);
%       if chan == param.ref_chan
%         median_mask = ones(size(good_rlines));
%         ref_phase_mean(chan) = 1;
%       else
%         % Remove outliers and take mean (real and imaginary done separately)
%         std_val = std(ref_values_real(good_rlines));
%         mean_val = mean(ref_values_real(good_rlines));
%         median_mask = ref_values_real > mean_val - std_val ...
%           & ref_values_real < mean_val + std_val;
%         
%         std_val = std(ref_values_imag(good_rlines));
%         mean_val = mean(ref_values_imag(good_rlines));
%         median_mask = median_mask & ref_values_imag > mean_val - std_val ...
%           & ref_values_imag < mean_val + std_val;
%         
%         ref_values_real_mean = mean(ref_values_real(good_rlines & median_mask));
%         ref_values_imag_mean = mean(ref_values_imag(good_rlines & median_mask));
%         
%         ref_phase_mean(chan) = ref_values_real_mean + j*ref_values_imag_mean;
%       end
      fprintf('WF %d: relative phase: %1.3f radians, %3.1f deg\n', chan, ...
        angle(ref_phase_mean(chan)), angle(ref_phase_mean(chan))*180/pi);
      
      if param.plot_en
        figure(20+chan); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Ang %d',chan));
        plot(angle(ref_phase));
        ylim([-pi pi]);
      end
      ref_phase(~(good_rlines & median_mask)) = NaN;
      if param.plot_en
        hold on;
        plot(angle(ref_phase),'ro');
        hold off;
        title(sprintf('Relative Phase (%d to ref %d)', chan, param.ref_chan));
      end
    else
      ref_phase_mean(chan) = 1;
    end
  end
  DDS_error = angle(ref_phase_mean) * 65536/(2*pi);
  if strcmpi(param.DDS_start_phase_units,'DDS')
    fprintf('Recommended new DDS phase settings (DDS counts):\n');
    new_DDS_phase = exp(j*param.DDS_start_phase/65536*2*pi - ref_phase_mean);
  elseif strcmpi(param.DDS_start_phase_units,'deg')
    fprintf('Recommended new DDS phase settings (deg):\n');
    new_DDS_phase = param.DDS_start_phase/360*2*pi - angle(ref_phase_mean);
  end
  fprintf('%.1f\t', new_DDS_phase(1:end-1)*180/pi);
  fprintf('%.1f', new_DDS_phase(end)*180/pi);
  fprintf('\n');
  fprintf('========================================================\n');
  results.DDS_phase_error(file_idx,:) = ref_phase_mean;
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
fprintf('Original');
for chan = 1:length(param.DDS_start_time)
  fprintf('\t%.2f', param.DDS_start_time(chan)*1e9);
end
fprintf('\n');
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
fprintf('Original');
for chan = 1:length(param.DDS_start_mag)
  fprintf('\t%.0f', param.DDS_start_mag(chan));
end
fprintf('\n');
fprintf('Mean');
final_DDS_amp = mean(results.DDS_amp,1);
% Normalize DDS amplitude to the maximum relative to the max_DDS_amp vector
if length(max_DDS_amp) == 1
  % All DDS's have the same maximum
  final_DDS_amp = final_DDS_amp / max(final_DDS_amp(~param.bad_chan_mask)) * max_DDS_amp;
else
  % Each DDS has a separate maximum
  [~,max_DDS_idx] = max(final_DDS_amp(~param.bad_chan_mask) ./ max_DDS_amp(~param.bad_chan_mask));
  good_chan_idxs = find(~param.bad_chan_mask);
  max_DDS_idx = good_chan_idxs(max_DDS_idx);
  final_DDS_amp = final_DDS_amp / final_DDS_amp(max_DDS_idx) * max_DDS_amp(max_DDS_idx);
end
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
    fprintf('\t%.1f', 180/pi*angle(results.DDS_phase_error(file_idx,wf)));
  end
  fprintf('\n');
end
fprintf('DDS_phase (deg)\n');
for file_idx = 1:length(fns)
  fprintf('%s', fns{file_idx});
  for wf = 1:size(results.DDS_phase,2)
    fprintf('\t%.0f', 180/pi*results.DDS_phase(file_idx,wf));
  end
  fprintf('\n');
end
fprintf('Original');
for chan = 1:length(param.DDS_start_phase)
  fprintf('\t%.1f', param.DDS_start_phase(chan));
end
fprintf('\n');
fprintf('Mean');
final_DDS_phase = angle(mean(exp(j*results.DDS_phase),1))*180/pi;
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.1f', final_DDS_phase(:,wf));
end
fprintf('\n');
fprintf('Median');
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.1f', final_DDS_phase(wf) + angle(mean(exp(j*results.DDS_phase(:,wf)) .* exp(-j*final_DDS_phase(wf)/180*pi),1))*180/pi);
end
fprintf('\n');
fprintf('Stdev');
for wf = 1:size(results.DDS_time,2)
  fprintf('\t%.1f', 180/pi*std(angle(exp(j*(results.DDS_phase(:,wf) - mean(results.DDS_phase(:,wf)))) )));
end
fprintf('\n');

if update_delay
  fprintf('Delay Compensated Mean');
  final_DDS_phase_comp = final_DDS_phase + 360*(final_DDS_time/1e9 ...
    - param.DDS_start_time)*(param.pc_param.f0+param.pc_param.f1)/2;
  final_DDS_phase_comp = 180/pi*angle(exp(j*(final_DDS_phase_comp - final_DDS_phase_comp(param.ref_chan))/180*pi));
  final_DDS_phase_comp(logical(param.bad_chan_mask)) = 0;
  for wf = 1:size(results.DDS_time,2)
    fprintf('\t%.1f', final_DDS_phase_comp(:,wf));
  end
  fprintf('\n');
  final_DDS_phase = final_DDS_phase_comp;
end



%% Update XML settings structure
if update_amplitude
  settings_enc.(config_var_enc).(ram_var_enc)(~param.bad_chan_mask) = uint16(final_DDS_amp(~param.bad_chan_mask));
end
for wf = 1:length(settings_enc.(config_var_enc).Waveforms)
  if update_phase
    if all(isreal(param.wf_mapping))
      settings_enc.(config_var_enc).Waveforms(wf).(phase_var_enc)(~param.bad_chan_mask) ...
        = mod(double(final_DDS_phase(~param.bad_chan_mask) + phase_desired(~param.bad_chan_mask)), 360);
    else
      if mod(wf,2)
        settings_enc.(config_var_enc).Waveforms(wf).(phase_var_enc)(~param.bad_chan_mask) ...
          = double(final_DDS_phase(~param.bad_chan_mask) + phase_desired(~param.bad_chan_mask));
      else
        settings_enc.(config_var_enc).Waveforms(wf).(phase_var_enc)(~param.bad_chan_mask) = ...
          angle(exp(j*(double(final_DDS_phase(~param.bad_chan_mask) + phase_desired(~param.bad_chan_mask))/180*pi+pi/2)))*180/pi;
      end
    end
  end
  if update_delay
    settings_enc.(config_var_enc).Waveforms(wf).Delay(~param.bad_chan_mask) ...
      = double(final_DDS_time(~param.bad_chan_mask) + time_delay_desired(~param.bad_chan_mask));
  end
end

%% Write XML file
[xml_fn_dir xml_fn_name xml_fn_ext] = fileparts(settings.fn);
out_xml_fn = fullfile(out_xml_fn_dir, sprintf('txequal_%s%s', xml_fn_name, xml_fn_ext));

settings_enc = rmfield(settings_enc,'fn');
settings_enc = rmfield(settings_enc,'datenum');
settings_enc = rmfield(settings_enc,'FPGAZ20Configuration');

if isfield(settings_enc,'xmlversion') && str2double(settings_enc.xmlversion{1}.values) >= 2.0
  settings_enc.sys.DDCZ20Ctrl = settings_enc.DDCZ20Ctrl;
  settings_enc.sys.DDSZ5FSetup = settings_enc.DDSZ5FSetup;
  settings_enc.sys.XMLZ20FileZ20Path = settings_enc.XMLZ20FileZ20Path;
  settings_enc.sys.xmlversion = settings_enc.xmlversion;
  settings_enc = rmfield(settings_enc,'DDCZ20Ctrl');
  settings_enc = rmfield(settings_enc,'DDSZ5FSetup');
  settings_enc = rmfield(settings_enc,'XMLZ20FileZ20Path');
  settings_enc = rmfield(settings_enc,'xmlversion');
end

fprintf('Writing %s\n', out_xml_fn);
fid = fopen(out_xml_fn,'w');
fprintf(fid,'<?xml version=''1.0'' standalone=''yes'' ?>\n');
fprintf(fid,'<LVData xmlns="http://www.ni.com/LVData">\n');
write_ni_xml_object(settings_enc,fid,true,struct('array_list','Waveforms','enum_list','DDCZ20sel'));
fprintf(fid,'</LVData>');
fclose(fid);

return;








