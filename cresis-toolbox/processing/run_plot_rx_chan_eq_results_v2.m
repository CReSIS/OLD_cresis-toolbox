% script run_plot_rx_chan_eq_results
%
% Script to plot the results of the receiver equalization process obtained
% using SAR images (f-k images) on a chunk/frame by chunk/frame basis for either a
% frame, segment, a day or an entire season
%
% Author: Peng Seng Tan
%

clear tmp
close all
clc;

% Single segment for 2011 Greenland P3
clear param;
param.radar_name = 'mcords2';
param.season_name = '2011_Greenland_P3';
param.day_seg = '20110506_01';
rx_chan_results_path = ct_filename_out(param, '', 'CSARP_rx_chan_eq');

% Using Waveform 01 and surface layer
% filenames_frames = get_filenames(rx_chan_results_path,'Data_rx_eq_01_2011','','.mat');

% Using Waveform 02 and bottom layer
filenames_frames = get_filenames(rx_chan_results_path,'Data_rx_eq_02_2011','','.mat');

td_out = [];
amp_out = [];
phase_out = [];

td_ave = [];
amp_ave = [];
phase_ave = [];

frame_roll = {};
frame_gps_time = {};
frame_peak_val = {};
frame_peak_offset = {};
frame_id = [];
frame_ref_adc_idx = [];

frame_count = 0;
for array_fn = filenames_frames.'
  tmp = load(array_fn{1});
  frame_count = frame_count + 1;

%   td_out = [td_out transpose(tmp.td_ave)];
%   amp_out = [amp_out transpose(tmp.amp_ave)];
%   phase_out = [phase_out transpose(tmp.phase_ave)];
  size_tmp_out = size(tmp.td_out);
  
  frame_roll{end+1} = tmp.chunk_roll;
  frame_gps_time{end+1} = tmp.chunk_gps_time;
  frame_peak_val{end+1} = tmp.peak_val;
  frame_peak_offset{end+1} = tmp.peak_offset;
  frame_id = [frame_id;array_fn{1}(end-18:end-4)];
  frame_ref_adc_idx = [frame_ref_adc_idx tmp.ref_idx(1)];
  
  for index = 1:size_tmp_out(2)
    td_out = [td_out transpose(tmp.td_out{index}(1:end))];
    amp_out = [amp_out transpose(tmp.amp_out{index}(1:end))];
    phase_out = [phase_out transpose(tmp.phase_out{index}(1:end))];
  end

end

size_results = size(td_out);

finite_index = isfinite(td_out);
num_finite_data = length(td_out(finite_index));

td_ave = mean(reshape(td_out(finite_index),size_results(1),num_finite_data/size_results(1)),2);
amp_ave = mean(reshape(amp_out(finite_index),size_results(1),num_finite_data/size_results(1)),2);
phase_ave = mean(reshape(phase_out(finite_index),size_results(1),num_finite_data/size_results(1)),2);

% Plot the results of the td, amp and phase compensation for each chunk under each frame

% colors: used for plotting
colors = {'-b','-r','-g','-c','-m','-k','--b','--r','--g','--c','--m','--k'};
colors1 = {'.b','.r','.g','.c','.m','.k','ob','or','og','oc','om','ok'};

figure(1)
clf(figure(1))
ha1 = axes;
for chan = 1:size_results(1)
%   plot(1:frame_count,phase_out(chan,1:end),[colors{mod(chan-1,length(colors))+1},'.']);
  plot(1:size_results(2),phase_out(chan,1:end),[colors{mod(chan-1,length(colors))+1},'.']);
  leg_cell{chan} = sprintf('rx %d',chan);
  hold on
end
title(sprintf('Phase offset for each receiver channel under segment: %s',param.day_seg),'Interpreter','none');
% xlabel('Frame ID');
xlabel('Chunk #');
ylabel('Phase offset (degrees)');
legend(ha1,leg_cell);
hold off

figure(2)
clf(2)
ha2 = axes;
for chan = 1:size_results(1)
%   plot(1:frame_count,amp_out(chan,1:end),[colors{mod(chan-1,length(colors))+1},'.']);
    plot(1:size_results(2),amp_out(chan,1:end),[colors{mod(chan-1,length(colors))+1},'.']);
hold on
end
title(sprintf('Amplitude offset for each receiver channel under segment: %s',param.day_seg),'Interpreter','none');
% xlabel('Frame ID');
xlabel('Chunk #');
ylabel('Amp offset (dB)');
legend(ha2,leg_cell);
hold off

figure(3)
clf(3)
ha3 = axes;
for chan = 1:size_results(1)
%   plot(1:frame_count,td_out(chan,1:end)*1e9,[colors{mod(chan-1,length(colors))+1},'.']);
  plot(1:size_results(2),td_out(chan,1:end)*1e9,[colors{mod(chan-1,length(colors))+1},'.']);
  hold on
end
title(sprintf('Time delay for each receiver channel under segment: %s',param.day_seg),'Interpreter','none');
% xlabel('Frame ID');
xlabel('Chunk #');
ylabel('Time delay (ns)');
legend(ha3,leg_cell);
hold off

% Plotting the Phase offset and roll angle of each good range line under
% each frame
size_frame_num = size(frame_roll);
std_phase_angle = zeros(size_results(1),size_frame_num(2));
mean_phase_angle = zeros(size_results(1),size_frame_num(2));
mean_complex = zeros(size_results(1),size_frame_num(2));

for i = 1:size_frame_num(2)
  figure(3+i*2)
  clf(3+i*2)
  ha4 = axes;
  for chan = 1:size_results(1)
    tmp2 = angle(frame_peak_val{i}(chan,:)./frame_peak_val{i}(frame_ref_adc_idx(i),:))*180/pi;
    std_phase_angle(chan,i) = sqrt(var(tmp2));
    mean_phase_angle(chan,i) = mean(tmp2);
    
    tmp3 = (frame_peak_val{i}(chan,:)./frame_peak_val{i}(frame_ref_adc_idx(i),:));
    mean_complex(chan,i) = mean(tmp3);
    tmp4 = angle(tmp3.*exp(-1i*angle(mean_complex(chan,i))))*180/pi;
    std_phase_angle1(chan,i) = sqrt(var(tmp4));
    
    plot(frame_gps_time{i},tmp4+angle(mean_complex(chan,i))*180/pi,[colors1{mod(chan-1,length(colors1))+1}]);
    hold on
  end
  title(sprintf('Phase offset for Good Range Lines under Frame ID: %s',frame_id(i,:)),'Interpreter','none');
  xlabel('GPS Time (seconds) for Good Range Lines');
  ylabel('Phase offset angle (degrees)');
  legend(ha4,leg_cell);
  hold off
  grid;
  
  figure(4+i*2)
  clf(4+i*2)
  plot(frame_gps_time{i},frame_roll{i}*180/pi,'.');
  title(sprintf('Roll angle for Good Range Lines under Frame ID: %s',frame_id(i,:)),'Interpreter','none');
  xlabel('GPS Time (seconds) for Good Range lines');
  ylabel('Roll angle (degrees)');
  grid;

end

% All segments for 2011 Greenland P3
% clear param;
% param.radar_name = 'mcords2';
% param.season_name = '2011_Greenland_P3';
% param.day_seg = '20110331_02';
% rx_chan_results_path = ct_filename_out(param, '', 'CSARP_rx_chan_eq',1');
% dirnames = get_filenames(rx_chan_results_path,'2011','','',struct('type','d'));
