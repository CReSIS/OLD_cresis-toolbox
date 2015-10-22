% script basic_radiometric_impulse_response
%
% This script is for helping with estimating the radiometric accuracy
% and impulse response of the system for each channel.
%
% Author: John Paden

clear param;
clear pc_param;
physical_constants;

% =======================================================================
% User Settings: Common Settings
% =======================================================================

params = [];
radar_name = 'mcords5';
fs = 1e9/2;
utc_time_correction = 2;
ref_fn_dir = '/mnt/products/mcords4_deconv/ice_shelf/';

if 0
  %% Calgary
  remove_bad_rlines = false;
  
  % Specify pulse compression properties for this waveform
  params(1).pc_param.f0 = 200e6;
  params(1).pc_param.f1 = 450e6;
  params(1).pc_param.Tpd = 3e-6;
  params(1).pc_param.tukey = 0.2;
  params(1).pc_param.window_func = inline('tukeywin_trim(N,0.1)');
  
  % .img = which waveform/adc pairs to load
  params(1).img = cat(2,-j*9*ones(1,1),[1].');
  
  % .rlines = Range lines to process from dataset
  %   These are range lines post presumming
  params(1).rlines = [1 inf];
  
  % .rbins = Range bins to search for surface in
  %   These range bins are post presumming
  params(1).rbins = [3000 3500];
  
  % .noise_rbins,rlines: Use for noise power calculation
  params(1).noise_rbins = [1000 1500];
  params(1).noise_rlines = [1 inf];
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  [params(:).gps_fn] = deal('/cresis/projects/dev/cr1/gps/2013_Antarctica_Basler/gps_20130921.mat');
  
  % base_path = Base path of data (does not include seg directory)
  % set = Which segment directory to load from, leave empty for no segment directory
  % utc_time_correction
  base_path = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/test_flights_20130921/'; % Calgary
  seg = '';
  
  % Optionally restrict search to a particular acquisition number/time
  % (part of the data files' filenames)
  acquisition_num = '20130921_20*03';
  
  % File index in filename
  file_nums = [0 1 2];
  
  % Transmit weights sent to lever_arm_fh
  Z0 = 50;
  params(1).tx_weights = [0 0 0 0 1 0 0 0] * sqrt(250*Z0);
  
  % Map of adcs to receivers for each waveform, used with lever_arm_fh
  params(1).rx_paths{9} = [1 3 5 7 2 4 6 8];
  
  Z0 = 50;
  adc_bits = 12;
  Vpp_scale = 2;
  adc_SNR_dB = 57;
  rx_gain = 10^((50)/20);
  noise_figure = 10^(1.6/10); % Do not include receiver losses
  
elseif 0
  %% 20131216 McMurdo Test flight
  remove_bad_rlines = true;
  params(1).raw_pulse_bins = 5800:6300;
  params(2).raw_pulse_bins = 4800:6300;
  params(3).raw_pulse_bins = 1350:6300;
  
  % Specify pulse compression properties for this waveform
  params(1).pc_param.f0 = 200e6;
  params(1).pc_param.f1 = 450e6;
  params(1).pc_param.Tpd = 1e-6;
  params(1).pc_param.tukey = 0.01;
  params(1).pc_param.window_func = inline('tukeywin_trim(N,0.1)');
  params(2).pc_param.f0 = 200e6;
  params(2).pc_param.f1 = 450e6;
  params(2).pc_param.Tpd = 3e-6;
  params(2).pc_param.tukey = 0.01;
  params(2).pc_param.window_func = inline('tukeywin_trim(N,0.1)');
  params(3).pc_param.f0 = 200e6;
  params(3).pc_param.f1 = 450e6;
  params(3).pc_param.Tpd = 10e-6;
  params(3).pc_param.tukey = 0.01;
  params(3).pc_param.window_func = inline('tukeywin_trim(N,0.1)');
  
  % .img = which waveform/adc pairs to load
  params(1).img = cat(2,-j*5*ones(8,1),[1 2 3 4 5 6 7 8].');
  params(2).img = cat(2,-j*3*ones(8,1),[1 2 3 4 5 6 7 8].');
  params(3).img = cat(2,-j*1*ones(8,1),[1 2 3 4 5 6 7 8].');
  
  % .rlines = Range lines to process from dataset
  %   These are range lines post presumming
  params(1).rlines = [1 inf];
  params(2).rlines = [1 inf];
  params(3).rlines = [1 inf];
  
  % .rbins = Range bins to search for surface in
  %   These range bins are post presumming
  params(1).rbins = [3100 3500];
  params(2).rbins = [3100 3500];
  params(3).rbins = [2500 3500];
  
  % .noise_rbins,rlines: Use for noise power calculation
  params(1).noise_rbins = [1000 1500];
  params(2).noise_rbins = [1000 1500];
  params(3).noise_rbins = [1000 1500];
  params(1).noise_rlines = [1 inf];
  params(2).noise_rlines = [1 inf];
  params(3).noise_rlines = [1 inf];
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  [params(:).gps_fn] = deal('/mnt/products/csarp_support/gps/2013_Antarctica_Basler/gps_20131216.mat');
  
  % base_path = Base path of data (does not include seg directory)
  % set = Which segment directory to load from, leave empty for no segment directory
  % utc_time_correction
  % base_path = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/test_flights_20130921/'; % Calgary
  base_path = '/mnt/backup-iu/array1/20131216/mcords4/';
  seg = '';
  
  % Optionally restrict search to a particular acquisition number/time
  % (part of the data files' filenames)
  acquisition_num = '20131216_04*02';
  
  % File index in filename
  file_num = [53];
  
  % Transmit weights sent to lever_arm_fh
  Z0 = 50;
  params(1).tx_weights = [9826 17844 32780 35264 28835 24243 23180 5515]/60000 * sqrt(250*Z0);
  params(2).tx_weights = [9826 17844 32780 35264 28835 24243 23180 5515]/60000 * sqrt(250*Z0);
  params(3).tx_weights = [9826 17844 32780 35264 28835 24243 23180 5515]/60000 * sqrt(250*Z0);
  
  % Map of adcs to receivers for each waveform, used with lever_arm_fh
  params(1).rx_paths{1} = [1 3 5 7 2 4 6 8];
  params(2).rx_paths{3} = [1 3 5 7 2 4 6 8];
  params(3).rx_paths{5} = [1 3 5 7 2 4 6 8];
  
  Z0 = 50;
  adc_bits = 12;
  Vpp_scale = 2;
  adc_SNR_dB = 57;
  rx_gain = 10^((50-30)/20);
  noise_figure = 10^(1.6/10); % Do not include receiver losses

else
  %% 20131223 Ice Shelf flight
  remove_bad_rlines = false;
  params(1).raw_pulse_bins = 1700:2200;
  params(2).raw_pulse_bins = 4800:6300;
  
  % Specify pulse compression properties for this waveform
  params(1).pc_param.f0 = 200e6;
  params(1).pc_param.f1 = 450e6;
  params(1).pc_param.Tpd = 3e-6;
  params(1).pc_param.tukey = 0.01;
  params(1).pc_param.window_func = inline('tukeywin_trim(N,0.1)');
  params(2).pc_param.f0 = 200e6;
  params(2).pc_param.f1 = 450e6;
  params(2).pc_param.Tpd = 1e-6;
  params(2).pc_param.tukey = 0.01;
  params(2).pc_param.window_func = inline('tukeywin_trim(N,0.1)');
  
  % .img = which waveform/adc pairs to load
  params(1).img = cat(2,-j*1*ones(8,1),[1 2 3 4 5 6 7 8].');
  params(2).img = cat(2,-j*3*ones(8,1),[1 2 3 4 5 6 7 8].');
  
  % .rlines = Range lines to process from dataset
  %   These are range lines post presumming
  params(1).rlines = [3144 3169];
  params(2).rlines = [2900 3300];
  
  % .rbins = Range bins to search for surface in
  %   These range bins are post presumming
  params(1).rbins = [1700 2200];
  params(2).rbins = [3100 3500];
  
  % .noise_rbins,rlines: Use for noise power calculation
  params(1).noise_rbins = [5000 5500];
  params(2).noise_rbins = [5000 5500];
  params(1).noise_rlines = [1 inf];
  params(2).noise_rlines = [1 inf];
  
  % .gps_fn = Optional GPS file name (leave empty to disable)
  [params(:).gps_fn] = deal('/mnt/products/csarp_support/gps/2013_Antarctica_Basler/gps_20131223.mat');
  
  % base_path = Base path of data (does not include seg directory)
  % set = Which segment directory to load from, leave empty for no segment directory
  % utc_time_correction
  % base_path = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/test_flights_20130921/'; % Calgary
  base_path = '/mnt/backup-iu/array1/20131223/mcords4/';
  seg = '';
  
  % Optionally restrict search to a particular acquisition number/time
  % (part of the data files' filenames)
  acquisition_num = '20131223_02*_03';
  
  % File index in filename
  file_num = [78];
  
  % Transmit weights sent to lever_arm_fh
  Z0 = 50;
  params(1).tx_weights = [8e+01 128.122 139.295 126.54 103.408 113.864 165.485 47.1897];
  params(2).tx_weights = [8e+01 128.122 139.295 126.54 103.408 113.864 165.485 47.1897];
  
  % Map of adcs to receivers for each waveform, used with lever_arm_fh
  params(1).rx_paths{1} = [1 3 5 7 2 4 6 8];
  params(2).rx_paths{3} = [1 3 5 7 2 4 6 8];
  
  Z0 = 50;
  adc_bits = 12;
  Vpp_scale = 2;
  adc_SNR_dB = 57;
  rx_gain = 10^((50-30)/20);
  noise_figure = 10^(1.6/10); % Do not include receiver losses
end

% =======================================================================
%% User Settings: Usually not changed
% =======================================================================

% param.type sent to motion_comp.m
[params(:).mocomp_type] = deal(4);

% lever_arm_fh = lever arm
[params(:).lever_arm_fh] = deal(@lever_arm);

% .plot_en = flag to enable plots
[params(:).plot_en] = deal(true);

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR are included in the estimate)
[params(:).snr_threshold] = deal(12);

% presums = Number of presums (coherent averaging) to do
presums = 1;

% .averaging_fh = method to average complex vectors when computing
%   recommended channel compensation (mean is ideal, but median
%   may be necessary to remove outliers)
[params(:).averaging_fh] = deal(@mean);

% Specify if cross correlation should be used (required for finding td)
% or if peak finding should be used when comparing channels
[params(:).cross_correlation_flag] = deal(0);

% Combine channels for surface tracker?
[params(:).combine_channels] = deal(false);

% =======================================================================
%% Automated Section
% =======================================================================

for param_idx = 1:length(params)
  for param_img = 1:size(params(param_idx).img,1)
    param = params(param_idx);
    param.img = param.img(param_img,:);
    ref = [];
    
    clear data;
    clear num_rec;
    epri_intersect = [];
    
    if strcmpi(radar_name,'mcords4')
  % adcs: a list of the adcs that we are loading
  adcs = unique(param.img(:,2));
  
  for adc = reshape(adcs,[1 length(adcs)])
    file_prefix = sprintf('mcords4_%02d_',adc);
    if isempty(acquisition_num)
      file_suffix = sprintf('%04d.bin',file_num);
    else
      file_suffix = sprintf('%s_%04d.bin',acquisition_num,file_num);
    end
    fn_dir = fullfile(base_path, sprintf('chan%d',adc), seg);
    fprintf('  Path: %s\n', fn_dir);
    fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
    fn = get_filename(fn_dir, file_prefix, '', file_suffix);
    fprintf('  Loading file %s\n', fn);
    % Load the data file
    [hdr,data_tmp] = basic_load_mcords4(fn,struct('clk',fs/4));
    % Map each of the read waveforms needed into the correct output
    for wf_adc_idx = 1:size(param.img,1)
      % wf,adc: pair of values for this entry in param.img
      wf = param.img(wf_adc_idx,1);
      if adc == abs(param.img(wf_adc_idx,2));
        % This pair needs to be loaded, insert into output array... handle
        % mismatched EPRIs using intersect function.
        if isempty(epri_intersect)
          epri_intersect = hdr.epri;
          if imag(wf) == 0
            data(:,:,wf_adc_idx) = data_tmp{wf}(:,:);
          else
            data(:,:,wf_adc_idx) = data_tmp{abs(wf)}(:,:) ...
              + sign(wf) * data_tmp{abs(wf)+1}(:,:);
          end
          hdr_utc_time_sod = hdr.utc_time_sod;
        else
          [epri_intersect data_idx data_tmp_idx] = intersect(epri_intersect,hdr.epri);
          data = data(:,data_idx,:);
          hdr_utc_time_sod = hdr.utc_time_sod(data_tmp_idx);
          if imag(wf) == 0
            data(:,:,wf_adc_idx) = data_tmp{wf}(:,data_tmp_idx);
          else
            data(:,:,wf_adc_idx) = data_tmp{abs(wf)}(:,data_tmp_idx) ...
              + sign(wf) * data_tmp{abs(wf)+1}(:,data_tmp_idx);
          end
        end
      end
    end
  end
  hdr.utc_time_sod = hdr_utc_time_sod;
  clear data_tmp;
  
  if any(strcmpi(radar_name,{'mcords2','mcords4'}))
    % =======================================================================
    % Presumming/coherent averaging
    % =======================================================================
    if presums > 1
      fprintf('Coherent averaging (%s)\n', datestr(now,'HH:MM:SS'));
      clear data_out;
      for adc_idx = 1:size(data,3)
        data_out(:,:,adc_idx) = fir_dec(data(:,:,adc_idx),presums);
      end
      data = data_out;
      clear data_out;
    end
    
    % =========================================================================
    % Pulse Compression
    % =========================================================================
    fprintf('Pulse compression (%s)\n', datestr(now,'HH:MM:SS'));
        orig_data_size = size(data);
        clear data_pc;
        pc_param = param.pc_param;
    pc_param.decimate = true;
    pc_param.zero_pad = true;
    for wf_adc_idx = 1:size(data,3)
      % wf,adc: The wf/adc pair of the current wf_adc_idx
      wf = abs(param.img(wf_adc_idx,1));
      adc = param.img(wf_adc_idx,2);
      % pc_param.time: fast time axis for this set of raw data
      pc_param.time = hdr.wfs(wf).t0 + (0:size(data,1)-1)/fs;
      % Pulse compress
      [data_pc(:,:,wf_adc_idx),hdr.wfs(wf).time,hdr.wfs(wf).freq] = pulse_compress(data(:,:,wf_adc_idx),pc_param);
    end
    
    % =========================================================================
        % Remove bad records
        % =========================================================================
        if remove_bad_rlines
          warning('Bad data removal routine running (should only be run on 20131216?)');
          raw_pulse_bins = param.raw_pulse_bins;
          
          % Works for file mcords4_01_20131216_042523_02_0053
          [surf_val surf_bin] = max(data_pc(param.rbins(1):param.rbins(end),:));
          surf_bin = surf_bin + param.rbins(1) - 1;
          imagesc(lp(data_pc));
          hold on;
          plot(surf_bin);
          hold off;
          
          %         pc_param.time(5945) - hdr.wfs(wf).time(3174)
          %         pc_param.time(5922) - hdr.wfs(wf).time(3174)
          
          check_bin = interp1(pc_param.time, 1:length(pc_param.time), ....
            hdr.wfs(wf).time(surf_bin) + 1.9600e-07*pc_param.Tpd/1e-6,'nearest','extrap');
          check_bins = sub2ind(size(data),check_bin.',1:size(data,2));
%           figure(1); clf;
%           plot(lp(real(data(check_bins)),2),'x')
%           hold on;
%           plot(lp(imag(data(check_bins)),2),'r.')
%           hold off;
          good_rlines = lp(real(data(check_bins)),2) > 30 & lp(imag(data(check_bins)),2) > 30;
          %good_rlines = lp(data(check_bins),2) > 40
          
          good_rlines2 = ones(size(good_rlines));
          locs = [];
          for rline = 1:size(data,2)
            [val,lags] = xcorr(real(data(raw_pulse_bins,rline)),imag(data(raw_pulse_bins,rline)));
            %               plot(abs(val))
            [~,locs(1)] = max(abs(val));
            val(locs(1)) = 0;
            [~,locs(2)] = max(abs(val));
            %               locs
            %               lags(locs)
            %               pause
            if any(lags(locs) == 0)
              good_rlines2(rline) = 0;
            end
          end
%           figure(1); clf;
%           plot(good_rlines,'.')
%           hold on;
%           plot(good_rlines2,'ro')
          
          good_rlines_idx = find(good_rlines & good_rlines2);
          
          figure(1); clf;
          imagesc(lp(data(:,good_rlines_idx)));
          ylim(raw_pulse_bins([1 end]))
          drawnow;
          pause(1);
          
          data_pc = data_pc(:,good_rlines_idx);
          data = data(:,good_rlines_idx);
        end
        
        % =======================================================================
        % Convert from quantization to voltage @ ADC
        data_pc = data_pc ...
          * Vpp_scale/2^adc_bits ...
          * 2^hdr.wfs(abs(param.img(1,1))).bit_shifts / hdr.wfs(abs(param.img(1,1))).presums;
        data = data ...
          * Vpp_scale/2^adc_bits ...
          * 2^hdr.wfs(abs(param.img(1,1))).bit_shifts / hdr.wfs(abs(param.img(1,1))).presums;
        
%         %imagesc(lp(data));
%         Pn = mean(abs(data(param.noise_rbins(1):param.noise_rbins(2),:)).^2);
%         figure(2); clf;
%         plot(10*log10(Pn));
%         % Power at radar RF input:
%         Pr_noise = 10*log10(mean(Pn)/Z0) + 30 - 10*log10(rx_gain)
%         Pr_adc = 10*log10((Vpp_scale/2/sqrt(2)).^2/Z0) + 30 - adc_SNR_dB - 10*log10(rx_gain)
%         Pr_noise_expect = 10*log10(BoltzmannConst * 290 * 2 * 250e6 * 10*log10(noise_figure)) + 30
        
        
        %clear data;
        % =========================================================================
    % Select ROI using ginput
    % =========================================================================
    figure(101); clf;
    imagesc(lp(data_pc(:,:,1)));
%         if use_ginput
%           fprintf('\nSelect Signal Region (2 clicks)\n');
%           title('Select Signal Region (2 clicks)');
%           [param.rlines,param.rbins] = ginput(2);
%           fprintf('\nSelect Noise Region (2 clicks)\n');
%           title('Select Noise Region (2 clicks)');
%           [param.noise_rlines,param.noise_rbins] = ginput(2);
%         end
    
  end
  
  % =======================================================================
  % Load GPS
  % =======================================================================
  if ~isempty(param.gps_fn)
    fprintf('Loading GPS (%s)\n', datestr(now,'HH:MM:SS'));
    gps = load(param.gps_fn);
    
    if strcmpi(radar_name,'mcrds')
      % Isolate the section of radar time from gps.radar_time that will be
      % used to interpolate with.
      guard_time = 5;
      good_idxs = find(gps.comp_time >= hdr.comp_time(1)-guard_time ...
        & gps.comp_time <= hdr.comp_time(end)+guard_time);
      good_radar_time = gps.radar_time(good_idxs);
      good_sync_gps_time = gps.sync_gps_time(good_idxs);
      % From these good indexes, remove any repeat radar times (usually caused
      % by there being more than one NMEA string every 1 PPS
      good_idxs = 1+find(diff(good_radar_time) ~= 0);
      good_radar_time = good_radar_time(good_idxs);
      good_sync_gps_time = good_sync_gps_time(good_idxs);
      % Interpolate gps.sync_gps_time to records.gps_time using gps.radar_time
      % and records.radar_time
      hdr.gps_time = interp1(good_radar_time, good_sync_gps_time, ...
        hdr.radar_time,'linear','extrap');
      
      hdr.lat = double(interp1(gps.gps_time,gps.lat,hdr.gps_time));
      hdr.lon = double(mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),hdr.gps_time)*180/pi+180, 360)-180);
      hdr.elev = double(interp1(gps.gps_time,gps.elev,hdr.gps_time));
      hdr.roll = double(interp1(gps.gps_time,gps.roll,hdr.gps_time));
      hdr.pitch = double(interp1(gps.gps_time,gps.pitch,hdr.gps_time));
      hdr.heading = double(mod(interp1(gps.gps_time,unwrap(gps.heading),hdr.gps_time)+pi,2*pi)-pi);
      
    else
      finfo = fname_info_mcords2(fn);
      [year,month,day] = datevec(finfo.datenum);
      hdr.gps_time = datenum_to_epoch(datenum(year,month,day,0,0,hdr.utc_time_sod)); % Still UTC time
      hdr.gps_time = hdr.gps_time + utc_leap_seconds(hdr.gps_time(1)); % Convert from UTC to GPS
      
      if exist('utc_time_correction','var')
        hdr.gps_time_corr = hdr.gps_time + utc_time_correction;
      else
        hdr.gps_time_corr = hdr.gps_time;
      end
      
      hdr.lat = interp1(gps.gps_time, gps.lat, hdr.gps_time_corr);
      hdr.lon = interp1(gps.gps_time, gps.lon, hdr.gps_time_corr);
      hdr.elev = interp1(gps.gps_time, gps.elev, hdr.gps_time_corr);
      hdr.roll = interp1(gps.gps_time, gps.roll, hdr.gps_time_corr);
      hdr.pitch = interp1(gps.gps_time, gps.pitch, hdr.gps_time_corr);
      hdr.heading = interp1(gps.gps_time, gps.heading, hdr.gps_time_corr);
    end
    hdr.gps_time = fir_dec(hdr.gps_time,presums);
    hdr.lat = fir_dec(hdr.lat,presums);
    hdr.lon = fir_dec(hdr.lon,presums);
    hdr.elev = fir_dec(hdr.elev,presums);
    hdr.roll = fir_dec(hdr.roll,presums);
    hdr.pitch = fir_dec(hdr.pitch,presums);
    hdr.heading = fir_dec(hdr.heading,presums);
    hdr.gps_source = gps.gps_source;
    
    figure(100); clf;
    plot(gps.gps_time,gps.roll*180/pi,'r');
    hold on;
    plot(hdr.gps_time,hdr.roll*180/pi);
    hold off;
    xlim([hdr.gps_time(1)-50 hdr.gps_time(end)+50])
    grid on;
    ylabel('Roll (deg)');
  else
    hdr.gps_time = zeros(1,size(data_pc,2));
    hdr.lat = zeros(1,size(data_pc,2));
    hdr.lon = zeros(1,size(data_pc,2));
    hdr.elev = zeros(1,size(data_pc,2));
    hdr.roll = zeros(1,size(data_pc,2));
    hdr.pitch = zeros(1,size(data_pc,2));
    hdr.heading = zeros(1,size(data_pc,2));
  end
    end
    
    
    if remove_bad_rlines
      hdr.elev = hdr.elev(good_rlines_idx);
    else
      data_pc = data_pc(:,param.rlines(1):param.rlines(end));
      hdr.elev = hdr.elev(param.rlines(1):param.rlines(end));
    end
        
    Mt = 20;
    dt = (hdr.wfs(wf).time(2) - hdr.wfs(wf).time(1))/Mt;
    data_pc2 = interpft(data_pc,Mt*size(data_pc,1));
    Nt = size(data_pc2,1);
    time = hdr.wfs(wf).time(1) + dt*(0:Nt-1).';
    
    [max_val max_idx] = max(data_pc2(param.rbins(1)*Mt:param.rbins(end)*Mt,:));
    max_idx = param.rbins(1)*Mt + max_idx - 1;
    plot(angle(max_val),'.')
    
    
    
    % Grx = Receiver gain per channel
    % Nf = Noise figure per channel
    % Pt = Total transmit power
    % Gt = Transmit gain
    % Gr = Receive gain
    % R = range to target (assumed to be specular)
    % gamma = power reflection coefficient
    Pt = sum(abs(param.tx_weights).^2/Z0);
    Gt = 10^(12/10);
    % Pt = 400/8; % Calgary
    % Gt = 10^(12/10)/8; % Calgary
    % Gr = Gt; % Calgary
    Gr = Gt/8;
    R = time(max_idx)*c/2;
    R = mean(R);
    fc = 0.5*(pc_param.f0 + pc_param.f1)
    lambda = c/fc;
    gamma = 10^(-3/10);
    Z0 = 50;
    %Grx = 10^(50/10); % Calgary
    Grx = 10^((50-30)/10);
    
    ref.Pr = 10*log10(max(abs(max_val)).^2/Z0 / Grx)
    
    ref.Pr_expected = 10*log10(Pt*Gt / (4*pi*(2*R)^2) * Gr * lambda.^2 / (4*pi) * gamma)
    
    time(max_idx);
    
    if 0
      dd = unwrap(angle(max_val));
      ee = diff(dd);
      actual_angle = [0 cumsum(ee)];
    end
    
    elev = hdr.elev - hdr.elev(1);
    
    kz = 4*pi*fc/c;

    unwrapped_angle = unwrap(angle(max_val));
%     unwrapped_angle = actual_angle(good_rlines_idx);
    distance = unwrapped_angle / -kz;
    distance = distance - distance(1);
    
    range = max_idx*dt * c/2;
    range = range - range(1);


%     dd = diff((elev-distance) * kz);
%     dd(abs(dd) > 1) = 0;
%     ee = [0 cumsum(dd)];
%     plot(elev - ee/kz)
%     td = (elev - ee/kz) / (c/2);


    td = distance / (c/2);
%     td = elev / (c/2);
    
    plot(distance)
    hold on;
    plot(elev,'r')
    plot(range,'g')
    hold off;
    
    BW = 1/dt;
    df = BW/Nt;
    freq = hdr.wfs(wf).freq(1) + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    
    for rline = 1:length(distance)
      data_pc2(:,rline) = 1/abs(max_val(rline)) * ifft(fft(data_pc2(:,rline)) .* exp(j*2*pi*freq*td(rline)));
    end
    
    [new_max_val new_max_idx] = max(data_pc2(param.rbins(1)*Mt:param.rbins(end)*Mt,:));
    new_max_idx = param.rbins(1)*Mt + new_max_idx - 1;
    
    bad_mask = angle(new_max_val) < -1;
    new_max_val = new_max_val(~bad_mask);
    
    new_max_idx = new_max_idx - new_max_idx(1);
    
    for rline = 1:length(distance)
      data_pc2(:,rline) = circshift(data_pc2(:,rline), [-new_max_idx(rline) 0]);
    end
    
    ref_fn = fullfile(ref_fn_dir, ...
      sprintf('ref_Tpd_%dus_adc_%d.mat',round(param.pc_param.Tpd*1e6),param.img(1,2)));
    fprintf('Writing file %s\n', ref_fn);
    %,'time','freq','signal'
    ref.param = param;
    ref.data = data_pc2(max_idx(1)+(-5000:5000),:);
    ref.freq = freq;
    ref.time = time;
    ref.adc = param.img(1,2);
    ref.Tpd = param.pc_param.Tpd;
    save(ref_fn,'-struct','ref');
    keyboard
    
  end
end

return;


in_fn_dir = '/mnt/products/';

  fns = get_filenames(in_fn_dir,'ref_Tpd_1us','','.mat');

  for fn_idx = 1:length(fns)
    ref = load(fns{fn_idx});
    ref.time = reshape(ref.time,[length(ref.time) 1]);
    ref.ref = mean(ref.data,2);
    % Find the peak
    [peak_val peak_idx] = max(ref.ref);
    ref.ref = ref.ref / peak_val;
    dt = ref.time(2) - ref.time(1);
    Nt = size(ref.data,1);
    ref.time = dt*(0:Nt-1).';
    fc = ref.freq(1);
    df = 1/(Nt*dt);
    ref.freq = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    ref.ref = ifft(fft(ref.ref) .* exp(-j*2*pi*ref.freq*ref.time(peak_idx)));
    ref.ref(1501:end-1500) = 0;
    
    Hwind = tukeywin_trim(3000,0.1);
    
    ref.ref(1:1500) = ref.ref(1:1500) .* Hwind(1500:-1:1);
    ref.ref(end-1499:end) = ref.ref(end-1499:end) .* Hwind(1:1500);
    
    freq_idxs = find(ref.freq >= 205e6 & ref.freq <= 445e6);
    Hwind = hanning(length(freq_idxs));
    Hwind = ifftshift(Hwind);
    ref_wind = zeros(size(ref.ref));
    ref_wind(freq_idxs) = Hwind;
    ref.fref = ref_wind ./ fft(ref.ref);
    ref.fref(isnan(ref.fref)) = 0;
    
    ff = ifft(fft(ref.data) .* repmat(ref.fref,[1 size(ref.data,2)]));
    
    figure(1); clf;
    plot(lp(ref.ref));
    figure(2); clf;
    imagesc(lp(ifft(fft(ref.data) .* repmat(ref.fref,[1 size(ref.data,2)]))))
    title(fns{fn_idx},'Interpreter','none');
    drawnow;
    save(fns{fn_idx}, '-struct', 'ref');
  end











