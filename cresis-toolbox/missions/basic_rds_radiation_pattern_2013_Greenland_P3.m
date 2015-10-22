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
beam.rbins = [4000 7000];

% .bin_rng = Range of range bins to grab relative to the threshold
beam.bin_rng = [0:40];

% .imgs = Choose which [wf adc] pairs to include in the radiation pattern
% generation
% beam.imgs = {[3 3; 1 3; 2 3; 4 3; 5 3; 6 3; 7 3; 8 3; 1 1; 2 2; 4 4; 5 5; 6 6; 7 7]};
beam.imgs = {[1 1; 1 2; 1 3; 1 4; 1 5; 1 6; 1 7]};

% Which files to use in radiation pattern generation
param.radar_name = 'mcords3';

% Single transmitter mode
param.vectors.file.start_idx = 2;
param.vectors.file.stop_idx = 265;
param.vectors.file.file_midfix = '_06_';

% Ping pong all transmitters mode (left/right beam)
% param.vectors.file.start_idx = 2;
% param.vectors.file.stop_idx = 123;
% param.vectors.file.file_midfix = '_07_';

param.vectors.file.base_dir = '/cresis/snfs1/data/MCoRDS/2013_Greenland_P3/';
param.vectors.file.adc_folder_name = '/20130315/board%b/';
param.vectors.file.file_prefix = 'mcords3_';
param.records.file.adcs = [1:7];

% .gps_fn = GPS file name
beam.gps_fn = '/cresis/projects/dev/cr1/gps/2013_Greenland_P3/gps_20130315.mat';

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
  
  param.radar.wfs(wf).rx_paths = [1 2 3 4 5 6 7];
  
  param.radar.wfs(wf).adc_gains = 10.^((52-[0 0 0 0 0 0 0])/20);
  param.radar.wfs(wf).chan_equal_dB = [0 0 0 0 0 0 0];
  param.radar.wfs(wf).chan_equal_deg = [0 0 0 0 0 0 0];
end

param.vectors.gps.time_offset = 0;

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
epri = [];
seconds = [];
fractions = [];
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
for file_idx = file_idxs
  img = 1;
  
  for wf_adc_idx = 1:size(beam.imgs{img},1)
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
      epri = cat(2,epri,hdr.epri);
      seconds = cat(2,seconds,hdr.seconds);
      fractions = cat(2,fractions,hdr.fractions);
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
      end
      
      surf_bins{wf_adc_idx} = cat(2,surf_bins{wf_adc_idx},surf_bins_tmp);
    else
      for rline = 1:size(surf_data,2)
        surf_vals_tmp(:,rline) = data(surf_bins_tmp(rline)+beam.bin_rng,rline);
      end
    end
    
    surf_vals{wf_adc_idx} = cat(2,surf_vals{wf_adc_idx},surf_vals_tmp);
    
    
    if beam.plot_en
      figure(wf_adc_idx); clf;
      imagesc(lp(surf_data));
      hold on;
      plot(surf_bins_tmp,'k');
      hold off;
      figure(wf_adc_idx+100); clf;
      subplot(2,1,1);
      plot(lp(surf_vals{wf_adc_idx}));
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
rds_radiation_pattern_fn = '/cresis/snfs1/scratch1/yzhu/rds_radiation_pattern_single.mat';
save(rds_radiation_pattern_fn,'lat','lon','elev','roll','pitch','heading','surf_vals','surf_bins','epri','seconds','fractions','gps_time','gps_source','param','beam','hdr');

return


rds_radiation_pattern_fn = '/home/cresis1/rds_radiation_pattern.mat';
rds_radiation_pattern_fn = '/home/cresis1/rds_radiation_pattern_pingpong.mat';
tmp = load(rds_radiation_pattern_fn);
if 0
  gps = load(tmp.beam.gps_fn);
  gps_time_offset = 0;
  roll = interp1(gps.gps_time,gps.roll, tmp.gps_time+gps_time_offset);
end

rlines = 1:length(tmp.surf_vals{1});

figure(1); clf;
legend_txt = {};
ant_vals = [];
for ant_idx = 1:length(tmp.surf_vals)
  
  [sort_rolls sort_idxs] = sort(tmp.roll(rlines));
  
  surf_vals_filt = double(tmp.surf_vals{ant_idx}(rlines)) / 50;
  sort_vals = 10*log10(surf_vals_filt(rlines)) + 30;
  sort_vals = sort_vals(sort_idxs);
  
  ant_vals(:,ant_idx) = medfilt1(sort_vals,301);
  h = plot(sort_rolls*180/pi,ant_vals(:,ant_idx));
  set(h,'Color',color_modes(ant_idx,:));
  hold on;
  
  legend_txt{ant_idx} = sprintf('Antenna %d',ant_idx);
  
  
end
hold off;
xlim([-37 +37])
legend(legend_txt);
grid on;

roll_angle = sort_rolls;



















      figure(1); clf;
      subplot(2,1,1);
      plot(lp( medfilt1(abs(double(surf_vals{8})).^2 / 50,101) ));
      xlim([9e4 10e4]);
      grid on;
      subplot(2,1,2);
      plot(roll*180/pi);
      xlim([9e4 10e4]);
      grid on;
      drawnow;

% Filter data
surf_vals_filt = medfilt1(abs(double(surf_vals{8})).^2 / 50,101);

figure(1001); clf;
plot(roll*180/pi, 10*log10(surf_vals_filt) + 30);
grid on;
xlabel('Roll (deg)');
ylabel('Pr (dBm)');





noise_power = -105;
surf_vals2 = [];
for wf_adc_idx = 1:7
  surf_vals2(wf_adc_idx,:) = surf_vals{wf_adc_idx};
end
tx_snr = tx_powers ./ noise_power;
good_meas = lp(tx_snr) > param.snr_threshold;
good_rlines = zeros(size(param.rlines));
good_rlines(sum(good_meas) == size(tx_snr,1)) = 1;
good_rlines = logical(good_rlines);

num_good_rlines = sum(good_rlines);
fprintf('Number of good range lines: %d out of %d\n', num_good_rlines, length(good_rlines));



















% =======================================================================
% Load and process data
% =======================================================================

% =======================================================================
% Eventually add for loops for ant_idx = 1:length(param.ant) and
% file_idx = 1:length(param.file.data_file_nums)
%     for ant_idx = 1
surf_vals = {};
surf_bins = {};
for ant_idx = 1:length(param.ant)
  gps_time = [];
  lat = [];
  lon = [];
  elev = [];
  roll = [];
  pitch = [];
  heading = [];
  for file_idx = 1:length(param.file.data_file_nums)
    
    % Create short variables for convenience and ease of reading
    wf = param.ant(ant_idx).wf;
    adc = param.ant(ant_idx).adc;
    file_num = param.file.data_file_nums(file_idx);
    
    % Only loading one waveform at a time
    wfs = wf;
    
    
    fprintf('  Done (%.1f sec)\n', toc(basic_mcords_radiation_pattern_tstart));
    
    % =======================================================================
    % Remove DC-bias (and 1's complement offset) in a robust way
    % =======================================================================
    data = data - median(data(:,1));
    
    % =======================================================================
    % Convert from quantization to voltage @ rx input
    % =======================================================================
    data = data ...
      * Vpp_scale/2^adc_bits ...
      * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums ...
      / rx_gain;
    
    if 0
      % Only need to run for Mar 5, 2012 Greenland P3 data (or other
      % datasets where UTC time is not good)
      output_fn = 'D:\data\fix_crappy_utc_time_in_2012rolls.mat';
      hdrs = load(output_fn);
      hdr.utc_time_sod = hdrs.utc_time_sod_corrected(find(hdr.epri(1) == hdrs.epri) + (0:length(hdr.epri)-1));
    end
    
    gps_time = cat(2,gps_time,interp1(utc_sod, gps.gps_time, hdr.utc_time_sod + gps_correction_time));
    lat = cat(2,lat,interp1(utc_sod, gps.lat, hdr.utc_time_sod + gps_correction_time));
    lon = cat(2,lon,interp1(utc_sod, gps.lon, hdr.utc_time_sod + gps_correction_time));
    elev = cat(2,elev,interp1(utc_sod, gps.elev, hdr.utc_time_sod + gps_correction_time));
    roll = cat(2,roll,interp1(utc_sod, gps.roll, hdr.utc_time_sod + gps_correction_time));
    pitch = cat(2,pitch,interp1(utc_sod, gps.pitch, hdr.utc_time_sod + gps_correction_time));
    heading = cat(2,heading,interp1(utc_sod, gps.heading, hdr.utc_time_sod + gps_correction_time));
    
    % =======================================================================
    % Noise removal and presumming
    % =======================================================================
    %     basic_presum_noise_removal;
    
    % =======================================================================
    % Pulse compression
    % =======================================================================
    fprintf('Pulse compression (%.1f sec)\n', toc(basic_mcords_radiation_pattern_tstart));
    clear pc_param;
    pc_param.f0 = param.pc_param.f0;
    pc_param.f1 = param.pc_param.f1;
    pc_param.Tpd = param.pc_param.Tpd;
    pc_param.tukey = param.pc_param.tukey;
    pc_param.time = hdr.wfs(wf).t0 + (0:size(data,1)-1)/fs;
    [data,time] = pulse_compress(data,pc_param);
    
    if param.rbins(2) > size(data,1)
      rbins = param.rbins(1):size(data,1);
    else
      rbins = param.rbins(1):param.rbins(2);
    end
    
    % =======================================================================
    % Surface tracker
    % =======================================================================
    surf_data = filter2(ones(1,10),abs(data.^2));
    [surf_vals_tmp surf_bins_tmp] = max(surf_data(rbins,:));
    surf_bins_tmp = rbins(1)-1 + surf_bins_tmp;
    
    if file_idx  == 1
      surf_vals{ant_idx} = [];
      surf_bins{ant_idx} = [];
    end
    surf_vals{ant_idx} = cat(2,surf_vals{ant_idx},surf_vals_tmp);
    surf_bins{ant_idx} = cat(2,surf_bins{ant_idx},surf_bins_tmp);
    
    if param.plot_en
      figure(ant_idx); clf;
      imagesc(lp(data));
      hold on;
      plot(surf_bins_tmp,'k');
      hold off;
      figure(ant_idx+100); clf;
      subplot(2,1,1);
      plot(lp(surf_vals{ant_idx}));
      subplot(2,1,2);
      plot(roll*180/pi);
      drawnow;
    end
  end
  figure(ant_idx+100); clf;
  plot(medfilt1(double(lp(surf_vals{ant_idx}) - max(lp(surf_vals{ant_idx}))),101));
  hold on;
  plot(roll*180/pi,'r');
  hold off;
  grid on;
  
  % Filter data
  surf_vals_filt = medfilt1(abs(double(surf_vals{ant_idx})).^2 / 50,101);
  
  figure(1001); clf;
  plot(roll*180/pi, 10*log10(surf_vals_filt) + 30);
  grid on;
  xlabel('Roll (deg)');
  ylabel('Pr (dBm)');
  saveas(1001,sprintf('D:/data/Pattern_wf_%d.fig',wf));
end

% =======================================================================
% Save Results
% =======================================================================
gps_source = gps.gps_source;
save('D:/data/roll_measurements.mat','lat','lon','elev','roll','pitch','heading','surf_vals','surf_bins','gps_time','gps_source','param','hdr');

return;


% =======================================================================
% Cull bad data points (low SNR)
% =======================================================================

noise_power = mean(mean(abs(data{param.ref_wf}(param.noise_rbins,param.rlines)).^2));
for wf = 1:5
  for rline_idx = 1:length(param.rlines)
    rline = param.rlines(rline_idx);
    tx_phases(wf,rline_idx) = data(surf_bins(rline_idx),rline);
    tx_powers(wf,rline_idx) = abs(data(surf_bins(rline_idx),rline)).^2;
  end
end
tx_snr = tx_powers ./ noise_power;
good_meas = lp(tx_snr) > param.snr_threshold;
good_rlines = zeros(size(param.rlines));
good_rlines(sum(good_meas) == size(tx_snr,1)) = 1;
good_rlines = logical(good_rlines);

num_good_rlines = sum(good_rlines);
fprintf('Number of good range lines: %d out of %d\n', num_good_rlines, length(good_rlines));

