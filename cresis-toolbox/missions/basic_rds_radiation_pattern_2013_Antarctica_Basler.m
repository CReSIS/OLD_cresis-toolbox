% script mcords_radiation_pattern
%
% This script generates a radiation pattern from raw radar data
% and INS data.
%
% General Operation:
% 1. Set up waveform mapping and other parameters
% 2. Run code
%
% Author: Kyle Byers, John Paden

% README: 
% Three sets of files were created "ALL", "TX", and "BEAMS"
% See these keywords below 

% =======================================================================
% User Settings
% =======================================================================

clear param;
clear beam;

param.radar.fs = 1e9/2;

% .plot_en = flag to enable plots
beam.plot_en = false;

% .rbins = Range bins to search for surface in
beam.rbins = [6000 12000];

% .bin_rng = Range of range bins to grab relative to the threshold
beam.bin_rng = [0:40];

% .imgs = Choose which [wf adc] pairs to include in the radiation pattern generation
% beam.imgs:
%   wf 1-8 is transmit single element 1-8
%   adc 1-8 map to rx antennas [1 3 5 7 2 4 6 8] (e.g. adc 2 is ant 3, adc 3 is rx 5, etc)
%   Antenna mapping is that used in the toolbox (see lever_arm.m)
%   Each of the entries in beam.imgs is in the format of a waveform-adc pair [-j*WAVEFORM ADC]

%beam.imgs = {[-7j 4; -7j 1; -7j 2; -7j 3; -7j 5; -7j 6; -7j 7; -7j 8; -1j 1; -3j 2; -5j 3; -9j 5; -11j 6; -13j 7; -13j 8]}; % ALL
%beam.imgs = {[-1j 4; -3j 4; -5j 4; -7j 4; -9j 4; -11j 4; -13j 4; -15j 4]}; % TX
beam.imgs = {[-1j 4; -3j 4; -5j 4]}; % BEAMS

% Which files to use in radiation pattern generation
param.radar_name = 'mcords4';
%param.vectors.file.start_idx = 1; %0 Individaul antenna mode
%param.vectors.file.stop_idx = 99; %100
param.vectors.file.start_idx = 1; %0 Three beam mode
param.vectors.file.stop_idx = 63; %64
param.vectors.file.base_dir = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_Basler/20131216/';
param.vectors.file.adc_folder_name = 'chan%d/';
param.vectors.file.file_prefix = 'mcords4_';
%param.vectors.file.file_midfix = '_02_'; % ALL and TX
param.vectors.file.file_midfix = '_02[234]*_00_'; % Three BEAM mode
param.records.file.adcs = [1:8];
param.records.file_version = 5;

% .gps_fn = GPS file name
beam.gps_fn = '/cresis/snfs1/dataproducts/csarp_support/gps/2013_Antarctica_Basler/gps_20131216.mat';

% .presums = Number of presums (coherent averaging) to do
beam.presums = 1;

param.radar.adc_bits = 14;
param.radar.Vpp_scale = 2;

% Setup waveforms
param.radar.wfs = [];
for wf = 1:16
  param.radar.wfs(wf).f0 = 200e6;
  param.radar.wfs(wf).f1 = 450e6;
  param.radar.wfs(wf).Tpd = 10e-6;
  param.radar.wfs(wf).tukey = 0.01;
  
  param.radar.wfs(wf).rx_paths = [1 2 3 4 5 6 7 8];
  
  param.radar.wfs(wf).adc_gains = 10.^((52-[25 25 25 25 25 25 25 25])/20);
  param.radar.wfs(wf).chan_equal_dB = [0 0 0 0 0 0 0 0];
  param.radar.wfs(wf).chan_equal_deg = [0 0 0 0 0 0 0 0];
end

param.vectors.gps.time_offset = 0;

%rds_radiation_pattern_fn = '/mnt/products/tmp/rds_radiation_pattern_2013_Antarctica_Basler_all.mat'; % ALL
%rds_radiation_pattern_fn = '/mnt/products/tmp/rds_radiation_pattern_2013_Antarctica_Basler_tx.mat'; %TX
rds_radiation_pattern_fn = '/mnt/products/tmp/rds_radiation_pattern_2013_Antarctica_Basler_beams.mat'; % BEAMS

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

fprintf('Starting radiation pattern generation (%s)\n', datestr(now));

% =======================================================================
% Load GPS data
% =======================================================================
gps = load(beam.gps_fn);
utc_sod = epoch_to_sod(gps.gps_time - utc_leap_seconds(gps.gps_time(1)), gps.gps_time(1));
gps_time = [];
lat = [];
lon = [];
elev = [];
roll = [];
pitch = [];
heading = [];
epri = {};
seconds = {};
fractions = {};
surf_vals = {};
surf_bins = {};

% =====================================================================
% Get the list of files to include in this records file
% =====================================================================
clear wfs hdrs;
boards = adc_to_board(param.radar_name,param.records.file.adcs);
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  fprintf('Getting files for board %d (%d of %d) (%s)\n', ...
    board, board_idx, length(boards), datestr(now));
  
  [base_dir,adc_folder_name,board_fns{board_idx},file_idxs] = get_segment_file_list(param,board);
end

% Parse files
img = 1;
for wf_adc_idx = 1:size(beam.imgs{img},1)
  epri{wf_adc_idx} = [];
  seconds{wf_adc_idx} = [];
  fractions{wf_adc_idx} = [];
  for file_idx = file_idxs
    
    wf = abs(beam.imgs{1}(wf_adc_idx,1));
    adc = beam.imgs{1}(wf_adc_idx,2);
    board = adc_to_board(param.radar_name,adc);
    
    fn = board_fns{find(boards==board)}{file_idx};
    [tmp fn_name] = fileparts(fn);
    
    fprintf('  wf %d adc %d file %s (%s)\n', wf, adc, fn_name, datestr(now))
    if strcmpi(param.radar_name,'mcords4')
      [hdr,data] = basic_load_mcords4(fn);
      if imag(beam.imgs{1}(wf_adc_idx,1)) ~= 0
        data = data{wf} + sign(beam.imgs{1}(wf_adc_idx,1)) * data{wf+1};
      end
      %data = data{wf}(:,:,mod(adc-1,4)+1);
    end
    
    if wf_adc_idx == 1
      % Check for day wrap
      day_wrap_idxs = find(diff(hdr.utc_time_sod) < -85000 & diff(hdr.utc_time_sod) > -88000);
      for day_wrap_idx = day_wrap_idxs
        hdr.utc_time_sod(day_wrap_idx+1) = hdr.utc_time_sod(day_wrap_idx+1) + 86400;
      end
      
      [year,month,day] = datevec(epoch_to_datenum(gps.gps_time));
      day_starts = datenum_to_epoch(unique(datenum(year,month,day)));
      
      day_start_found = false;
      for day_start_idx = 1:length(day_starts)
        hdr.gps_time = day_starts(day_start_idx) + hdr.utc_time_sod + param.vectors.gps.time_offset + utc_leap_seconds(gps.gps_time(1));
        test_gps_time = interp1(gps.gps_time, gps.gps_time, hdr.gps_time);
        if all(~isnan(test_gps_time))
          day_start_found = true;
          break;
        end
      end
      if ~day_start_found
        keyboard
      end
      hdr.gps_time = test_gps_time;
      gps_time = cat(2,gps_time,hdr.gps_time);
      lat = cat(2,lat,interp1(gps.gps_time, gps.lat, hdr.gps_time));
      lon = cat(2,lon,interp1(gps.gps_time, gps.lon, hdr.gps_time));
      elev = cat(2,elev,interp1(gps.gps_time, gps.elev, hdr.gps_time));
      roll = cat(2,roll,interp1(gps.gps_time, gps.roll, hdr.gps_time));
      pitch = cat(2,pitch,interp1(gps.gps_time, gps.pitch, hdr.gps_time));
      heading = cat(2,heading,interp1(gps.gps_time, gps.heading, hdr.gps_time));
      epri{wf_adc_idx} = cat(2,epri{wf_adc_idx},hdr.epri);
      seconds{wf_adc_idx} = cat(2,seconds{wf_adc_idx},hdr.seconds);
      fractions{wf_adc_idx} = cat(2,fractions{wf_adc_idx},hdr.fractions);
    else
      epri{wf_adc_idx} = cat(2,epri{wf_adc_idx},hdr.epri);
      seconds{wf_adc_idx} = cat(2,seconds{wf_adc_idx},hdr.seconds);
      fractions{wf_adc_idx} = cat(2,fractions{wf_adc_idx},hdr.fractions);
    end
    
    % =======================================================================
    % Remove DC-bias (and 1's complement offset) in a robust way
    % =======================================================================
    data = data - median(data(:,1));
    
    % =======================================================================
    % Convert from quantization to voltage @ rx input
    % =======================================================================
    data = data ...
      * param.radar.Vpp_scale/2^param.radar.adc_bits ...
      * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums ...
      / param.radar.wfs(wf).adc_gains(adc);
    
    % =======================================================================
    % Pulse compression
    % =======================================================================
    fprintf('    Pulse compression (%s)\n', datestr(now));
    clear pc_param;
    pc_param.f0 = param.radar.wfs(wf).f0;
    pc_param.f1 = param.radar.wfs(wf).f1;
    pc_param.Tpd = param.radar.wfs(wf).Tpd;
    pc_param.tukey = param.radar.wfs(wf).tukey;
    pc_param.time = hdr.wfs(wf).t0 + (0:size(data,1)-1)/param.radar.fs;
    [data,time] = pulse_compress(data,pc_param);
    
    if beam.rbins(2) > size(data,1)
      rbins = beam.rbins(1):size(data,1);
    else
      rbins = beam.rbins(1):beam.rbins(2);
    end
    
    % =======================================================================
    % Surface tracker
    % =======================================================================
    if file_idx == file_idxs(1)
      surf_vals{wf_adc_idx} = [];
      if wf_adc_idx == 1
        surf_bins{wf_adc_idx} = [];
      end
    end
    
    if wf_adc_idx == 1
      surf_data = filter2(ones(1,10),abs(data(rbins,:).^2));
      
      [max_vals max_bins] = max(surf_data);
      max_bins = rbins(1)-1 + max_bins;
      threshold = max_vals / 10;
      
      % Find first peak after leading edge
      surf_vals_tmp = zeros(length(beam.bin_rng),size(data,2));
      surf_bins_tmp = zeros(1,size(data,2));
      for rline = 1:size(surf_data,2)
        surf_bins_tmp(rline) = find(surf_data(:,rline) > threshold(rline),1);
        surf_bins_tmp(rline) = rbins(1)-1 + surf_bins_tmp(rline);
        surf_vals_tmp(:,rline) = data(surf_bins_tmp(rline)+beam.bin_rng,rline);
      end
      if any(surf_bins_tmp > max_bins) || any(surf_bins_tmp+beam.bin_rng(end) < max_bins)
        warning('Missed peak in %d %s\n', wf_adc_idx, fn);
        if 0
          rline = 35;
          plot(rbins, lp(surf_data(:,rline)))
          hold on
          plot(surf_bins_tmp(rline), lp(surf_data(surf_bins_tmp(rline)+1-rbins(1),rline)),'rx')
          plot(surf_bins_tmp(rline)+beam.bin_rng(end), lp(surf_data(surf_bins_tmp(rline)+1-rbins(1)+beam.bin_rng(end),rline)),'rx')
          hold off;
        end
      end
      
      surf_bins{wf_adc_idx} = cat(2,surf_bins{wf_adc_idx},surf_bins_tmp);
    else
      surf_vals_tmp = zeros(length(beam.bin_rng),size(data,2));
      for rline = 1:size(data,2)
        % Find the corresponding range line in wf_adc_idx==1
        match_rline = find(hdr.epri(rline) == epri{1});
        if isempty(match_rline)
          surf_vals_tmp(:,rline) = NaN;
        else
          surf_vals_tmp(:,rline) = data(surf_bins{1}(match_rline)+beam.bin_rng,rline);
        end
      end
    end
    
    surf_vals{wf_adc_idx} = cat(2,surf_vals{wf_adc_idx},surf_vals_tmp);
    
    
    if beam.plot_en
      figure(wf_adc_idx); clf;
      imagesc(lp(surf_data));
      hold on;
      plot(surf_bins_tmp - rbins(1) + 1,'k');
      hold off;
      figure(wf_adc_idx+100); clf;
      subplot(2,1,1);
      plot(max(lp(surf_vals{wf_adc_idx})));
      subplot(2,1,2);
      plot(roll*180/pi);
      drawnow;
    end
    
    
  end
end

% =======================================================================
% Save Results
% =======================================================================
gps_source = gps.gps_source;
save(rds_radiation_pattern_fn,'lat','lon','elev','roll','pitch','heading','surf_vals','surf_bins','epri','seconds','fractions','gps_time','gps_source','param','beam','hdr');

return
