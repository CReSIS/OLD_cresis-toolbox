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

fs = 1e9/9;

% .plot_en = flag to enable plots
param.plot_en = true;

% .rbins = Range bins to search for surface in
param.rbins = [1 inf];

% .snr_threshold = SNR threshold in dB (range lines exceeding this
%   SNR are included in the estimate)
param.snr_threshold = 10;

% .file = The directory where the data files are located
%    (THIS NEEDS TO BE SET EVERYTIME)
param.radar_name = 'mcords2';
if strcmpi(param.radar_name,'mcords')
  param.file.base_dirs            = '/mnt/';
  param.file.adc_folder_name      = 'mcords_new1/seg13/';
  param.file.file_prefix          = 'mcords.rec008';
  % Vector of data file numbers to load i.e. *.NNNN.dat
  param.file.data_file_nums       = [0:13];
else
  % .board = board number from 0 to 3
  param.board = 0;
  % .adc = the receive channel to use (relative to the board # so that it
  %    is always contained in [1,4] since there are 4 channels per board)
  param.adc = 4;
  param.acquisition_num = 8;
  
  % Parameters to locate specific file of interest
  % (THIS NEEDS TO BE SET EVERYTIME)
  param.seg = 'seg_06';
  param.base_path = 'D:\data\mcords\20120304\';

  % Vector of data file numbers to load i.e. *.NNNN.dat
  param.file.data_file_nums = 9:367;
end

% .gps_fn = GPS file name
param.gps_fn = 'C:/csarp_support/gps/2012_Greenland_P3/gps_20120305.mat';

% .presums = Number of presums (coherent averaging) to do
param.presums = 1;

% Generates an antenna gain pattern for each entry here, but the
% corresponding waveform and adc must be specified for each antenna.

if 0
  % 2011 Antarctica DC8 .ant Parameters
  param.ant           = [];
  param.ant(end+1).wf     = 1;    %Antenna 1
  param.ant(end).adc    = 1;
  param.ant(end+1).wf     = 2;    %Antenna 2
  param.ant(end).adc    = 2;
  param.ant(end+1).wf     = 3;    %Antenna 3
  param.ant(end).adc    = 3;
  param.ant(end+1).wf     = 4;    %Antenna 4
  param.ant(end).adc    = 4;
  param.ant(end+1).wf     = 5;    %Antenna 5
  param.ant(end).adc    = 5;
  param.ant(end+1).wf     = 6;    %Antenna 1-5
  param.ant(end).adc    = 3;
else
  % 2012 Greenland P3 .ant Parameters
  param.ant           = [];
  param.ant(end+1).wf     = 1;    %Antenna 1
  param.ant(end).adc    = 2;
%   param.ant(end+1).wf     = 2;    %Antenna 2
%   param.ant(end).adc    = 3;
%   param.ant(end+1).wf     = 3;    %Antenna 3
%   param.ant(end).adc    = 4;
%   param.ant(end+1).wf     = 4;    %Antenna 4
%   param.ant(end).adc    = 5;
%   param.ant(end+1).wf     = 5;    %Antenna 5
%   param.ant(end).adc    = 6;
%   param.ant(end+1).wf     = 6;    %Antenna 6
%   param.ant(end).adc    = 7;
%   param.ant(end+1).wf     = 7;    %Antenna 7
%   param.ant(end).adc    = 8;
end

% Setup waveforms
param.wfs = [];
param.pc_param.f0 = 180e6;
param.pc_param.f1 = 210e6;
% param.pc_param.f0 = 189.15e6;
% param.pc_param.f1 = 198.65e6;
param.pc_param.Tpd = 10e-6;
param.pc_param.tukey = 0.2;

adc_bits = 14;
Vpp_scale = 2;
adc_SNR_dB = 70;
rx_gain = 10^((72-19)/20);

gps_correction_time = -14;
gps_correction_time = 0;

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

basic_mcords_radiation_pattern_tstart = tic;

% =======================================================================
% Load GPS data
% =======================================================================
gps = load(param.gps_fn);
[year month day hour minute sec] = datevec(epoch_to_datenum(gps.gps_time(1)));
utc_sod = gps.gps_time - datenum_to_epoch(datenum(year,month,day,0,0,0)) - utc_leap_seconds(gps.gps_time(1));

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

    if strcmpi(param.radar_name,'mcords')
      % Create sub-folder name for the particular adc channel
      adc_idx_insert_idxs = strfind(param.file.adc_folder_name,'%d');
      mat_cmd = 'adc_folder_name = sprintf(param.file.adc_folder_name';
      for adc_idx_insert_idx = adc_idx_insert_idxs
        mat_cmd = [mat_cmd sprintf(', %d',adc)];
      end
      mat_cmd = [mat_cmd ');'];
      eval(mat_cmd);

      % Location of data files
      filepath = fullfile(param.file.base_dirs, adc_folder_name);

      % Load in data file
      tic;
      fprintf('Loading data\n');
      fn = get_filename(filepath,param.file.file_prefix,'',sprintf('%04d.dat',file_num));
      if isempty(fn)
        error('Could not find any files which match');
        return;
      end
      fprintf('  Loading file %s\n', fn);
      if file_num == 0
        [hdr,data] = basic_load_mcords(fn, struct('clk',1e9/9,'first_byte',2^26,'wfs',wfs));
      else
        [hdr,data] = basic_load_mcords(fn, struct('clk',1e9/9,'wfs',wfs));
      end
      data = data{1}(1:end-1,:);
    elseif strcmpi(param.radar_name,'mcords2')
      param.board = floor((adc-1)/4);
      board_adc = mod(adc-1,4)+1;
      if isempty(param.seg)
        fn_dir = fullfile(param.base_path, sprintf('board%d',param.board));
      else
        fn_dir = fullfile(param.base_path, sprintf('board%d',param.board), ...
          param.seg);
      end
      file_prefix = sprintf('mcords2_%d_',param.board);
      if isempty(param.acquisition_num)
        file_suffix = sprintf('%04d.bin',file_num);
      else
        file_suffix = sprintf('%02d_%04d.bin',param.acquisition_num,file_num);
      end
      fprintf('  Path: %s\n', fn_dir);
      fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
      fn = get_filename(fn_dir, file_prefix, '', file_suffix);
      if isempty(fn)
        fprintf('  Could not find any files which match\n');
        return;
      end
      fprintf('  Loading file %s\n', fn);
      [hdr,data] = basic_load_mcords2(fn,struct('clk',fs/2));
      data = data{wf}(:,:,board_adc);
    end
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

