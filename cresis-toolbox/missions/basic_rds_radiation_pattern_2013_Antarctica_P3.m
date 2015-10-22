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

% =======================================================================
% User Settings
% =======================================================================

clear param;
clear beam;

param.radar.fs = 1e9/9;

% .plot_en = flag to enable plots
beam.plot_en = false;

% .rbins = Range bins to search for surface in
beam.rbins = [3000 6000];

% .bin_rng = Range of range bins to grab relative to the threshold
beam.bin_rng = [0:40];

% .imgs = Choose which [wf adc] pairs to include in the radiation pattern generation
% beam.imgs:
%   wf 1-7 is transmit single element 1-7
%   wf 8 is transmit all elements with tx calibrated weights [1 1 1 1 1 1 1]
%   adc 2-16 map to rx antennas 1-15
%   antenna mapping is that used in the toolbox (see lever_arm.m)
%   Reference pattern is first wf-adc pair (this is the 3rd antenna)
%   Next 6 are with a common receiver (receiver 4 which is the 3rd antenna)
%   Next 6 entries after that are common transmit/receive. e.g. [1 2] is tx element 1, rx element 1
%   Next 8 entries after that are the outer elements
beam.imgs = {[3 4; 3 9; 3 10; 3 11; 3 12; 3 13; 3 14; 3 15; 3 16; 1 4; 2 4; 4 4; 5 4; 6 4; 7 4; 8 4; 1 2; 2 3; 4 5; 5 6; 6 7; 7 8]};

% Which files to use in radiation pattern generation
param.radar_name = 'mcords3';
param.vectors.file.start_idx = 19; %0
param.vectors.file.stop_idx = 107; %107
param.vectors.file.base_dir = '/cresis/snfs1/data/MCoRDS/2013_Antarctica_P3/';
param.vectors.file.adc_folder_name = '20131108/board%b/';
param.vectors.file.file_prefix = 'mcords3_';
param.vectors.file.file_midfix = '_00_';
param.records.file.adcs = [1:16];

% .gps_fn = GPS file name
beam.gps_fn = '/cresis/projects/dev/cr1/gps/2013_Antarctica_P3/gps_20131108.mat';

% .presums = Number of presums (coherent averaging) to do
beam.presums = 1;

param.radar.adc_bits = 14;
param.radar.Vpp_scale = 2;

% Setup waveforms
param.radar.wfs = [];
for wf = 1:8
  param.radar.wfs(wf).f0 = 180e6;
  param.radar.wfs(wf).f1 = 210e6;
  % param.radar.wfs(wf).f0 = 189.15e6;
  % param.radar.wfs(wf).f1 = 198.65e6;
  param.radar.wfs(wf).Tpd = 10e-6;
  param.radar.wfs(wf).tukey = 0.2;
  
  param.radar.wfs(wf).rx_paths = [1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
  
  param.radar.wfs(wf).adc_gains = 10.^((52-[0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0])/20);
  param.radar.wfs(wf).chan_equal_dB = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
  param.radar.wfs(wf).chan_equal_deg = [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0];
end

param.vectors.gps.time_offset = 0;

rds_radiation_pattern_fn = '/cresis/scratch1/paden/rds_radiation_pattern_2013_Antarctica_P3_all.mat';

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
boards = unique(floor((param.records.file.adcs-1)/4));
for board_idx = 1:length(boards)
  board = boards(board_idx);
  
  fprintf('Getting files for board %d (%d of %d) (%s)\n', ...
    board, board_idx, length(boards), datestr(now));
  
  [base_dir,adc_folder_name,board_fns{board_idx},file_idxs] = get_segment_file_list(param,(board+1)*4);
end

% Parse files
img = 1;
for wf_adc_idx = 1:size(beam.imgs{img},1)
  epri{wf_adc_idx} = [];
  seconds{wf_adc_idx} = [];
  fractions{wf_adc_idx} = [];
  for file_idx = file_idxs
    
    wf = beam.imgs{1}(wf_adc_idx,1);
    adc = beam.imgs{1}(wf_adc_idx,2);
    board = adc_to_board(param.radar_name,adc);
    
    fn = board_fns{find(boards==board)}{file_idx};
    [tmp fn_name] = fileparts(fn);
    
    fprintf('  wf %d adc %d file %s (%s)\n', wf, adc, fn_name, datestr(now))
    if strcmpi(param.radar_name,'mcords3')
      if file_idx == file_idxs(1)
        [hdr,data] = basic_load_mcords3(fn,struct('clk',param.radar.fs,'first_byte',0));
      else
        [hdr,data] = basic_load_mcords3(fn,struct('clk',param.radar.fs));
      end
      data = data{wf}(:,:,mod(adc-1,4)+1);
    end
    
    if wf_adc_idx == 1
      gps_time = cat(2,gps_time,interp1(utc_sod, gps.gps_time, hdr.utc_time_sod + param.vectors.gps.time_offset));
      lat = cat(2,lat,interp1(utc_sod, gps.lat, hdr.utc_time_sod + param.vectors.gps.time_offset));
      lon = cat(2,lon,interp1(utc_sod, gps.lon, hdr.utc_time_sod + param.vectors.gps.time_offset));
      elev = cat(2,elev,interp1(utc_sod, gps.elev, hdr.utc_time_sod + param.vectors.gps.time_offset));
      roll = cat(2,roll,interp1(utc_sod, gps.roll, hdr.utc_time_sod + param.vectors.gps.time_offset));
      pitch = cat(2,pitch,interp1(utc_sod, gps.pitch, hdr.utc_time_sod + param.vectors.gps.time_offset));
      heading = cat(2,heading,interp1(utc_sod, gps.heading, hdr.utc_time_sod + param.vectors.gps.time_offset));
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
