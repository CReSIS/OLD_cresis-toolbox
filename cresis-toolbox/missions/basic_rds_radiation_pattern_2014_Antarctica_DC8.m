% script basic_rds_radiation_pattern
%
% This script generates a radiation pattern from raw radar data
% and INS data.
%
% General Operation:
% 1. Set up waveform mapping and other parameters
% 2. Run code
% 3. 
%
% Author: Kyle Byers, John Paden

% =======================================================================
% User Settings
% =======================================================================

clear param;
clear beam;

param.radar.fs = 150e6;

% .plot_en = flag to enable plots
beam.plot_en = false;

% .rbins = Range bins to search for surface in
beam.rbins = [1200 2500];

% .bin_rng = Range of range bins to grab relative to the threshold
beam.bin_rng = [0:40];

% .imgs = Choose which [wf adc] pairs to include in the radiation pattern
% generation
% First is reference (center element)
% Then common receive for 6 individual elements (1-6 2)
% Then combined (7 2)
% Then rx/tx common (N N)
beam.imgs = {[2 2; 1 2; 3 2; 4 2; 5 2; 6 2; 7 2; 1 1; 3 3; 4 4; 5 5; 6 6]};

% Which files to use in radiation pattern generation
param.radar_name = 'mcords3';
param.vectors.file.base_dir = 'D:\testflights\';
param.vectors.file.adc_folder_name = 'board%b/';
param.vectors.file.file_prefix = 'mcords3_';
% First set
param.vectors.file.file_midfix = '_07_';
param.vectors.file.start_idx = 1;
param.vectors.file.stop_idx = 127;
output_fn = 'D:\roll_measurements_2014_Antarctica_DC8.mat';
param.vectors.gps.time_offset = 1;
% Second set
% param.vectors.file.file_midfix = '_08_';
% param.vectors.file.start_idx = 1;
% param.vectors.file.stop_idx = 154;
% output_fn = 'D:\roll_measurements_2014_Antarctica_DC8_2.mat';
% param.vectors.gps.time_offset = 1+86400;
% Third set
% param.vectors.file.file_midfix = '_07_';
% param.vectors.file.start_idx = 1;
% param.vectors.file.stop_idx = 127;
% param.vectors.file.file_midfix = '?';
param.records.file.adcs = [1:8];

% Additional fields required for sync_radar_to_gps to work
param.season_name = '2014_Antarctica_DC8';
param.day_seg = '20141007_01';
param.vectors.gps.fn = '';

% .presums = Number of presums (coherent averaging) to do
beam.presums = 1;

param.radar.adc_bits = 14;
param.radar.Vpp_scale = 2;

% Setup waveforms
param.radar.wfs = [];
for wf = 1:7
  param.radar.wfs(wf).f0 = 165e6;
  param.radar.wfs(wf).f1 = 215e6;
  % param.radar.wfs(wf).f0 = 189.15e6;
  % param.radar.wfs(wf).f1 = 198.65e6;
  param.radar.wfs(wf).Tpd = 10e-6;
  param.radar.wfs(wf).tukey = 0.1;
  
  param.radar.wfs(wf).rx_paths = [1 2 3 4 5 6];
  
  param.radar.wfs(wf).adc_gains = 10.^((52-[15 15 15 15 15 15 0 0])/20);
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
%gps = load(beam.gps_fn);
utc_time_sod = [];
epri = cell(1,size(beam.imgs{1},1));
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
    
    epri{wf_adc_idx} = cat(2,epri{wf_adc_idx},hdr.epri);
    if wf_adc_idx == 1
      utc_time_sod = cat(2,utc_time_sod, hdr.utc_time_sod);
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
      % We will capture from the first_peak +/- beam.bin_rng, we want to
      % check to see how often the actual maximum lies out of this region
      num_missed = sum(max_bins < surf_bins_tmp+beam.bin_rng(1) | max_bins > surf_bins_tmp+beam.bin_rng(end));
      if num_missed > 0
        warning('Missed %d of %d peaks in %d %s\n', num_missed, size(surf_data,2), wf_adc_idx, fn);
      end
      
      surf_bins{wf_adc_idx} = cat(2,surf_bins{wf_adc_idx},surf_bins_tmp);
    else
      for rline = 1:size(surf_data,2)
        if rline <= size(data,2)
          surf_vals_tmp(:,rline) = data(surf_bins_tmp(rline)+beam.bin_rng,rline);
        else
          surf_vals_tmp(:,rline) = 0;
        end
      end
    end
    
    surf_vals{wf_adc_idx} = cat(2,surf_vals{wf_adc_idx},surf_vals_tmp);
    
    
    if beam.plot_en
      figure(wf_adc_idx); clf;
      imagesc(lp(surf_data));
      hold on;
      plot(surf_bins_tmp - rbins(1)+1,'k');
      plot(max_bins - rbins(1)+1,'k.');
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
save(output_fn,'surf_vals','surf_bins','utc_time_sod','epri','seconds','fractions','gps_time','gps_source','param','beam','hdr');

records = sync_radar_to_gps(param,[],utc_time_sod);
save(output_fn,'-append','-struct','records');

return
