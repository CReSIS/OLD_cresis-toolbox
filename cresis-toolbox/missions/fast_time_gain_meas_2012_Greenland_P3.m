% Script fast_time_gain_meas_2012_Greenland_P3.m
%
% Plots the receiver gain versus fast time for each receiver.
% Also stores this into a Matlab file.
%
% Author: John Paden

% =========================================================================
%% Find fast time gain measurements
% This section was used to find out which files were used for the fast time
% gain measurements.
if 0
  % Search through all settings files
  path_dirs = get_filenames('/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/','2012[0-9]','','',struct('type','d'));
  
  for path_dirs_idx = 1:length(path_dirs)
    path_dir = path_dirs{path_dirs_idx};
    fprintf('%s\n', path_dir);
    settings = read_ni_xml_directory(path_dir,'DDS',0);
    % This is how we figured out which ones were fast time gain measurements:
    for settings_idx = 1:length(settings)
      for wf = 1:length(settings(settings_idx).Configuration.Waveforms)
        if settings(settings_idx).Configuration.Waveforms(wf).Start_Freq ...
            - settings(settings_idx).Configuration.Waveforms(wf).Stop_Freq == 0
          fprintf('  %d %d: %s\n', settings_idx, wf, settings(settings_idx).fn);
        end
      end
    end
  end

  % Search through the specific date that seems to be the only one with
  % fast time gain measurements
  path_dir = '/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/';
  settings = read_ni_xml_directory(path_dir,'DDS',0);
  
  % This is how we figured out which ones were fast time gain measurements:
  for settings_idx = 1:length(settings)
    for wf = 1%:length(settings(settings_idx).Configuration.Waveforms)
      if settings(settings_idx).Configuration.Waveforms(wf).Start_Freq == 195000000
        fprintf('%d: %s\n', settings_idx, settings(settings_idx).fn);
      end
    end
  end
  
  % Set the settings in browse_ni_xml_settings script to determine
  % which data files go with which settings
  % browse_ni_xml_settings


  %figure(1); clf;
  for board = 0:3
    % This is how we found the reference measurement (ADC 2 according
    % to the field report):
    %fns = get_filenames(sprintf('/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board%d/seg_04',board),'mcords2','','.bin');
    
    % This is how we found the file_adc_mapping mapping
    fns = get_filenames(sprintf('/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board%d/seg_05',board),'mcords2','','.bin');
    
    for fn_idx = 1:length(fns)
      fn = fns{fn_idx};
      [hdr,data] = basic_load_mcords2(fn,struct('clk',1e9/9,'recs',[0 1]));
      wf = 2;
      for sub_adc = 1:4
        if any(data{wf}(:,1,sub_adc) > 1000)
          % This is how we found the file_adc_mapping comments:
          fprintf('%s (%d): %d %d   %f\n', fn, fn_idx, board, sub_adc, max(data{wf}(:,1,sub_adc)));
          % This is how we found the file_adc_mapping variable:
          % fprintf('[%d %d %d]\n', fn_idx, board, sub_adc);
          if 0
            plot(data{wf}(:,1,sub_adc))
            pause;
          end
        end
      end
    end
    
  end
  
  % TAKEN WITH 9? dB OF PADDING and DIRECT CONNECTION TO ADC
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_04/mcords2_0_20120304_232626_00_0000.bin (1): 0 2   15633.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_04/mcords2_0_20120304_232639_00_0001.bin (2): 0 2   15639.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_04/mcords2_0_20120304_233539_01_0000.bin (3): 0 2   15520.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_04/mcords2_0_20120304_233552_01_0001.bin (4): 0 2   15528.000000
  
  % TAKEN WITH 59? dB OF PADDING and CONNECTION THROUGH EACH TO ADC
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_05/mcords2_0_20120305_131534_01_0010.bin (18): 0 2   4342.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_05/mcords2_0_20120305_131547_01_0011.bin (19): 0 2   4337.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_05/mcords2_0_20120305_131628_01_0012.bin (20): 0 3   4848.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_05/mcords2_0_20120305_131641_01_0013.bin (21): 0 3   4851.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_05/mcords2_0_20120305_131728_01_0014.bin (22): 0 4   6409.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_05/mcords2_0_20120305_131741_01_0015.bin (23): 0 4   6397.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_131850_01_0016.bin (24): 1 1   4218.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_131903_01_0017.bin (25): 1 1   4210.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_132000_01_0018.bin (26): 1 2   3958.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_132013_01_0019.bin (27): 1 2   3945.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_132112_02_0000.bin (28): 1 3   4235.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_132125_02_0001.bin (29): 1 3   4233.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_132243_03_0000.bin (30): 1 4   4238.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board1/seg_05/mcords2_1_20120305_132256_03_0001.bin (31): 1 4   4243.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_130913_01_0002.bin (10): 2 4   5269.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_130926_01_0003.bin (11): 2 4   5268.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_131043_01_0004.bin (12): 2 3   5962.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_131056_01_0005.bin (13): 2 3   5956.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_131207_01_0006.bin (14): 2 2   5451.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_131220_01_0007.bin (15): 2 2   5436.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_131329_01_0008.bin (16): 2 1   6840.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board2/seg_05/mcords2_2_20120305_131342_01_0009.bin (17): 2 1   6847.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130222_00_0000.bin (1): 3 4   5203.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130235_00_0001.bin (2): 3 4   5206.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130248_00_0002.bin (3): 3 4   5198.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130448_00_0003.bin (4): 3 3   6883.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130501_00_0004.bin (5): 3 3   6894.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130552_00_0005.bin (6): 3 2   5534.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130605_00_0006.bin (7): 3 2   5525.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130749_01_0000.bin (8): 3 1   5197.000000
  % /cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board3/seg_05/mcords2_3_20120305_130802_01_0001.bin (9): 3 1   5192.000000
  
end

% =========================================================================
%% Estimate the fast time gain measurements from the data

file_adc_mapping = [[18 0 2]
  [19 0 2]
  [20 0 3]
  [21 0 3]
  [22 0 4]
  [23 0 4]
  [24 1 1]
  [25 1 1]
  [26 1 2]
  [27 1 2]
  [28 1 3]
  [29 1 3]
  [30 1 4]
  [31 1 4]
  [10 2 4]
  [11 2 4]
  [12 2 3]
  [13 2 3]
  [14 2 2]
  [15 2 2]
  [16 2 1]
  [17 2 1]
  [1 3 4]
  [2 3 4]
  [3 3 4]
  [4 3 3]
  [5 3 3]
  [6 3 2]
  [7 3 2]
  [8 3 1]
  [9 3 1]];

adc_bits = 14;
Vpp_scale = 2.23;
test_gain = 10^((-(0))/20);
fc = 195e6;
rx_paths = [1 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15];
% [B,A] = butter(2,0.01); % Original figure parameters
[B,A] = butter(4,0.1);

% Reference power
ref_fn = '/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board0/seg_04/mcords2_0_20120304_233539_01_0000.bin';
fprintf('%s\n', ref_fn);
[hdr,data] = basic_load_mcords2(ref_fn,struct('clk',1e9/9,'recs',[0 1]));
board_adc = 2;
rline = 1;
for wf = 1:2
  data{wf} = data{wf}(:,:,board_adc) - median(data{wf}(:,1,board_adc));
  data{wf} = data{wf}(1:end-1,:);
  data{wf} = data{wf} ...
    * Vpp_scale/2^adc_bits ...
    * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums ...
    / test_gain;
  
  fs = 1e9/9;
  dt = 1/fs;
  Nt = size(data{wf},1);
  rx_time = dt*(0:Nt-1).';
  rx_baseband = filtfilt(B, A, 2*double(data{wf}(:,rline)) .* exp(-j*2*pi*fc*rx_time));
  
  if 0
    figure(101);
    plot(rx_time*1e6,lp(rx_baseband,2), 'k');
    hold on;
    grid on;
    xlabel('fast time (us)');
    ylabel('gain (dB)');
    ylim([-20 60]);
    xlim([0 40]);
    drawnow;
  end
  
  ref_power(wf) = mean(abs(rx_baseband(rx_time > 5e-6 & rx_time < 25e-6)).^2);
end

colors = {'k' 'r' 'y' 'g' 'c' 'b' 'm' 'k:' 'r:' 'y:' 'g:' 'c:' 'b:' 'm:' 'k--' 'r--' 'y--' 'g--' 'c--' 'b--' 'm--'};
label_str = {};
rx_gain = {};
rx_gain_time = {};
rlines = 2001:3000;
test_gain = [10^((-(59.5))/20) 10^((-(54.5))/20)];

for wf = 1:2
  figure(wf); clf;
end


ref_adc_power = 1.0e+05 * ...
  [0.776704757359404   0.760445918296173
   0.974403773234188   0.955191962507858
   1.691165656258723   1.658611430372293
   0.733991620073936   0.720106350710783
   0.653963339435462   0.634647899516298
   0.740284716955099   0.727008776217417
   0.741638644125219   0.727810529338735]
   
rx_gain = {};
rx_gain_time = {};
label_str = {};
adc_power = [];
adcs = 2:8;
for meas_idx = 1:length(adcs)
  adc = adcs(meas_idx);

  board = adc_to_board('mcords2',adc);
  board_adc = mod(adc-1,4)+1;
  
  found = false;
  for fn_idx = 1:length(file_adc_mapping)
    if all(file_adc_mapping(fn_idx,2:3) == [board board_adc])
      found = true;
      break;
    end
  end
  if ~found
    error('Could not find file for adc %d (%d,%d)', adc, board, board_adc);
  end
  fn_idx = file_adc_mapping(fn_idx,1);
  
  fns = get_filenames(sprintf('/cresis/snfs1/data/MCoRDS/2012_Greenland_P3/20120319/board%d/seg_05/',board),'mcords2','','.bin');

  fn = fns{fn_idx};
  fprintf('%d: %s\n', adc, fn);
  [hdr,data] = basic_load_mcords2(fn,struct('clk',1e9/9));
  for wf = 1:2
    data{wf} = data{wf}(:,:,board_adc) - median(data{wf}(:,1,board_adc));
    data{wf} = data{wf}(1:end-1,:);
    data{wf} = data{wf} ...
      * Vpp_scale/2^adc_bits ...
      * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums ...
      / test_gain(wf) / ref_power(wf) / sqrt(ref_adc_power(meas_idx,wf));
    
    fs = 1e9/9;
    dt = 1/fs;
    Nt = size(data{wf},1);
    rx_time = hdr.wfs(wf).t0 - 770e-9 + dt*(0:Nt-1).';
    rx_baseband = filtfilt(B, A, 2*double(mean(data{wf}(:,rlines),2)) .* exp(-j*2*pi*fc*rx_time));
    
    figure(wf);
    h(wf,meas_idx) = plot(rx_time*1e6,lp(rx_baseband,2), colors{meas_idx});
    label_str{wf,meas_idx} = sprintf('chan %d', meas_idx);
    hold on;
    grid on;
    xlabel('fast time (us)');
    ylabel('gain (dB)');
    ylim([-37 60]);
    drawnow;
    
    if wf == 1
      adc_power(meas_idx,wf) = mean(abs(rx_baseband(rx_time > 13e-6 & rx_time < 18e-6)).^2);
      xlim([1 11]);
    else
      adc_power(meas_idx,wf) = mean(abs(rx_baseband(rx_time > 30e-6 & rx_time < 35e-6)).^2);
      xlim([10 20]);
    end
    rx_gain{wf}(:,meas_idx) = reshape(rx_baseband, [length(rx_baseband) 1]);
    rx_gain_time{wf,meas_idx} = rx_time;
    
  end
end

for wf = 1:2
  legend(h(wf,:),label_str(wf,:));
end


Tpd = [1e-6 10e-6];
tukeywin_corr = 1.0;
hanning_corr = 5.9;
Nt_window = round(Tpd/dt);
last_good_time = [17.5e-6 30e-6];
figure(3); clf;
colors = {'r','b'};
for wf = 1:2;
  hold on;
  %plot(rx_gain_time{wf,1}*1e6, 20*log10(abs(rx_gain{wf})));
  %gain_correction{wf} = -(tukeywin_corr+20*log10(conv(mean(abs(rx_gain{wf}),2), tukeywin(Nt_window(wf),0.2) / Nt_window(wf), 'same')));
  gain_correction{wf} = -(hanning_corr+20*log10(conv(mean(abs(rx_gain{wf}),2), hanning(Nt_window(wf)) / Nt_window(wf), 'same')));
  gain_correction_time{wf} = (rx_gain_time{wf,1}-Tpd(wf)/2);
  last_good_idx = find(gain_correction_time{wf} > last_good_time(wf),1);
  gain_correction{wf}(last_good_idx+1:end) = gain_correction{wf}(last_good_idx);
  gain_correction{wf} = gain_correction{wf} - min(gain_correction{wf});
  plot(gain_correction_time{wf}*1e6, gain_correction{wf},colors{wf});
  grid on;
  xlabel('fast time (us)');
  ylabel('power gain correction (dB)');
  ylim([0 60]);
  xlim([0 20])
end
legend('1 us pulse', '10 us pulse');

notes = sprintf('These gain correction curves are for {1} = 1e-6 and {2} = 10e-6 \npulse durations from 2012 Greenland P3 fast time gain measurements on \nMar 4/5. Units are dB power correction and seconds two way travel time. \nGenerated on %s with %s.', datestr(now), mfilename());
save('post_pulse_comp_gain_correction.mat','gain_correction_time','gain_correction','notes');
saveas(3,'post_pulse_comp_gain_correction.fig');


%% Example of fast time gain correction application
% Note that fast time gain before 1-2 us is not accurate and is probably
% just amplifying noise and feed through most of the time.

% Load data
mdata = load('/cresis/snfs1/dataproducts/public/data/rds/2013_Greenland_P3/CSARP_standard/20130420_02/Data_20130420_02_008.mat');
% Load gain correction table
gc_table = load('post_pulse_comp_gain_correction.mat');
% Interpolate the gain correction table onto the radar data fast time axis
wf_gain_corr = {};
for wf = 1:2
  wf_gain_corr{wf} = interp1(gc_table.gain_correction_time{wf}, gc_table.gain_correction{wf}, mdata.Time);
end
first_idx = find(mdata.Time > gc_table.gain_correction_time{1}(1), 1);
last_idx = find(mdata.Time < gc_table.gain_correction_time{2}(end), 1, 'last');
% Plot radar echogram before gain correction
figure(1); clf;
imagesc([],mdata.Time,10*log10(mdata.Data));
% Apply gain correction one range line at a time
mdata.Data_GC = zeros(size(mdata.Data));
for rline = 1:size(mdata.Data,2)
  % Find the point where the two waveforms are combined
  wf_comb_time = find(mdata.Time > mdata.Surface(rline) + mdata.param_combine.combine.img_comb(1),1);
  % Create the combined waveform gain correction curve
  gain_corr = wf_gain_corr{1};
  gain_corr(wf_comb_time:end) = wf_gain_corr{2}(wf_comb_time:end);
  % Extend gain curve to first and last time on the edges so that we have 
  % finite values for all fast time.
  gain_corr(1:first_idx-1) = gain_corr(first_idx);
  gain_corr(last_idx+1:end) = gain_corr(last_idx);
  if 0
    % Plot for debugging purposes
    plot(mdata.Time, gain_corr)
    xlabel('time');
    ylabel('gain correction for combined (dB)');
  end
  % Apply the gain correction
  mdata.Data_GC(:,rline) = mdata.Data(:,rline) .* 10.^(gain_corr/10);
end
% Plot radar echogram after gain correction
figure(2); clf;
imagesc([],mdata.Time,10*log10(mdata.Data_GC));

% Compare a single range line/a-scope before and after gain correction
figure(3); clf;
plot(mdata.Time*1e6,10*log10(mdata.Data(:,1)))
hold on;
plot(mdata.Time*1e6,10*log10(mdata.Data_GC(:,1)),'r')
legend('Before','After')


