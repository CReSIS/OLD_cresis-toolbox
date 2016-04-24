function basic_rx_chan_equalization(param,defaults)
% basic_rx_chan_equalization(param,defaults)
%
% RUN THIS FUNCTION FROM "run_basic_rx_chan_equalization"
%
% This script is for helping with setting the receiver coefficients
% from raw data. It requires loading one waveform and N receive channels
% and then analyzing these data.
%
% 1. Collect data with so that receiver gain settings are the same as
%    the ordinary data collection, but the air/ice surface response
%    is unsaturated.  Flatter surfaces work best. If INS data is not
%    available, then the platform should be straight and level.
%    If the TR switch is slow, you will need to collect the data from
%    a high enough altitude where the fast-time gain changes of the TR
%    switch are stable (this is 10 us after the switch change for the
%    original MCoRDS system).
% 2. Pulse duration should not matter, but probably affects results
%    slightly
% 3. Time gate should be large enough to include noise-only data which
%    will be used to determine the SNR of the surface return (SNR
%    threshold is used to exclude low SNR data points from the measurement)
%
% Author: John Paden, Logan Smith
%
% See Also: run_basic_rx_chan_equalization

%% basic_rx_chan_equalization preparation

physical_constants;

global g_num_files;
if isempty(g_num_files)
  g_num_files = 1;
end
num_files = input(sprintf('How many files [%d]: ',g_num_files));
if isempty(num_files)
  num_files = g_num_files;
else
  g_num_files = num_files;
end

%% Process the files
file_name_list = {};
for file_idx = 1:num_files
  
  %% Load the files
  if file_idx == 1
    [data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults);
    file_name_list = {fn};
  else
    param.file_search_mode = 'default+1';
    global g_basic_noise_analysis_fn;
    g_basic_noise_analysis_fn = file_name_list{end};
    [data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults);
    file_name_list{end+1} = fn;
    g_basic_noise_analysis_fn = file_name_list{1};
  end
  
  %% Convert from quantization to voltage @ ADC
  wf = abs(param.img(1,1));
  data = data ...
    * default.radar.adc_full_scale/2^default.radar.adc_bits ...
    * 2^hdr.wfs(abs(wf)).bit_shifts / hdr.wfs(wf).presums;
  
  %% Additional software presums
  for wf_adc = 1:size(data,3)
    data(:,1:floor(size(data,2)/param.presums),wf_adc) = fir_dec(data(:,:,wf_adc),param.presums);
  end
  data = data(:,1:floor(size(data,2)/param.presums),:);
  hdr.radar_time = fir_dec(hdr.radar_time,param.presums);
  hdr.lat = fir_dec(hdr.lat,param.presums);
  hdr.lon = fir_dec(hdr.lon,param.presums);
  hdr.elev = fir_dec(hdr.elev,param.presums);
  hdr.roll = fir_dec(hdr.roll,param.presums);
  hdr.pitch = fir_dec(hdr.pitch,param.presums);
  hdr.heading = fir_dec(hdr.heading,param.presums);
  
  %% Pulse compression
  [pc_signal,pc_time] = pulse_compress(data,pc_param);

  %% Track surface
  surf_bin = NaN*zeros(1,size(pc_signal,2));
  ml_data = lp(fir_dec(abs(pc_signal(:,:,param.ref_wf_adc)).^2,ones(1,5)/5,1));
  for rline = 1:size(ml_data,2)
    % Threshold is hard coded to max of 7 dB of noise or 13 below peak
    cur_threshold = max([ml_data(1,rline)+7; ml_data(:,rline)-13]);
    tmp = find(ml_data(1:end-2,rline) > cur_threshold,1);
    if ~isempty(tmp)
      [~,max_offset] = max(ml_data(tmp+(0:2),rline));
      tmp = tmp-1 + max_offset;
      surf_bin(rline) = tmp;
    end
  end
  
  param.noise_rlines = 1:size(ml_data,2);
  param.noise_rbins = min(surf_bin)-100 : min(surf_bin)-40;
  
  param.rlines = 1:size(ml_data,2);
  param.rbins= min(surf_bin)-10 : max(surf_bin)+10;

  %% Perform receiver channel equalization
  param.averaging_fh = @mean;
  param.time = pc_time;
  dt = pc_time(2) - pc_time(1);
  Nt = length(pc_time);
  df = 1/(Nt*dt);
  param.freq = pc_param.DDC_freq + (-floor(Nt/2)*df : df : floor((Nt-1)/2)*df).';
  
  param.mocomp_type = 4;
  param.tx_weights = double(settings.DDS_Setup.Ram_Amplitude);
  param.rx_paths = {}; param.rx_paths{wf} = default.radar.rx_paths;
  param.lever_arm_fh = @lever_arm;
  
  param.combine_channels = false;
  param.snr_threshold = 10;
  param.phase = default.radar.wfs(1).chan_equal_dB;
  param.amp = default.radar.wfs(1).chan_equal_deg;
  param.td = default.radar.wfs(1).chan_equal_Tsys;

  [td_out(:,file_idx),amp_out(:,file_idx),phase_out(:,file_idx), full_out] = rx_chan_equal(pc_signal,param,hdr);
  
  %% Collate results
  if file_idx == 1
    peak_ref = full_out.peak_ref;
    roll = gps.roll;
    gps_time = gps.gps_time;
    peak_offset = full_out.peak_offset;
    roll_est_theta = full_out.roll_est_theta;
    roll_est_val = full_out.roll_est_val;
    roll_est_gps_time = full_out.gps_time;
  else
    peak_ref = cat(2,peak_ref,full_out.peak_ref);
    roll = cat(2,roll,gps.roll);
    gps_time = cat(2,gps_time,gps.gps_time);
    peak_offset = cat(2,peak_offset,full_out.peak_offset);
    roll_est_theta = cat(2,roll_est_theta,full_out.roll_est_theta);
    roll_est_val = cat(2,roll_est_val,full_out.roll_est_val);
    roll_est_gps_time = cat(2,roll_est_gps_time,full_out.gps_time);
  end
  
end

%% Print Results
fprintf('========================================================\n');
fprintf('Recommended equalization coefficients (averaged results)\n');

fprintf('  Date of processing: %s, mocomp %d, wf/adc %d/%d bins %d-%d\n', ...
  datestr(now), param.mocomp_type, param.img(param.ref_wf_adc,1), ...
  param.img(param.ref_wf_adc,2), param.rbins(1), param.rbins(end));
fprintf('td settings\n');
for file_idx = 1:num_files
  [~,fn] = fileparts(file_name_list{file_idx});
  fprintf('%s', fn);
  fprintf('\t%.2f', td_out(:,file_idx)*1e9);
  fprintf('\n');
end
fprintf('amp settings\n');
for file_idx = 1:num_files
  [~,fn] = fileparts(file_name_list{file_idx});
  fprintf('%s', fn);
  fprintf('\t%.1f', amp_out(:,file_idx));
  fprintf('\n');
end
fprintf('phase settings\n');
for file_idx = 1:num_files
  [~,fn] = fileparts(file_name_list{file_idx});
  fprintf('%s', fn);
  fprintf('\t%.1f', phase_out(:,file_idx));
  fprintf('\n');
end

td_ave = mean(td_out,2);
amp_ave = mean(amp_out,2);
phase_ave = angle(mean(exp(j*phase_out/180*pi),2))*180/pi;

fprintf('Rx Path\n');
for wf_adc = 1:size(param.img,1)
  wf = abs(param.img(wf_adc,1));
  adc = param.img(wf_adc,2);
  if wf_adc < size(param.img,1)
    fprintf('%d\t', param.rx_paths{wf}(adc));
  else
    fprintf('%d', param.rx_paths{wf}(adc));
  end
end
fprintf('\n');

fprintf('Original/Recommended/Difference td settings (ns):\n');
fprintf('%.2f\t', param.td(1:end-1)*1e9);
fprintf('%.2f', param.td(end)*1e9);
fprintf('\n');
fprintf('%.2f\t', td_ave(1:end-1)*1e9);
fprintf('%.2f', td_ave(end)*1e9);
fprintf('\n');
fprintf('%.2f\t', (td_ave(:).' - param.td(:).')*1e9);
fprintf('\n');

fprintf('Original/Recommended/Difference amp settings (dB):\n');
fprintf('%.1f\t', param.amp(1:end-1));
fprintf('%.1f', param.amp(end));
fprintf('\n');
fprintf('%.1f\t', amp_ave(1:end-1));
fprintf('%.1f', amp_ave(end));
fprintf('\n');
fprintf('%.1f\t', amp_ave(:).' - param.amp(:).');
fprintf('\n');

% Rewrap each phase so that the output does not print +355 deg and -5 deg
fprintf('Original/Recommended/Difference phase settings (deg):\n');
phase_rewrapped = angle(exp(j*param.phase/180*pi)) * 180/pi;
fprintf('%.1f\t', phase_rewrapped(1:end-1));
fprintf('%.1f', phase_rewrapped(end));
fprintf('\n');
phase_rewrapped = angle(exp(j*phase_ave/180*pi)) * 180/pi;
fprintf('%.1f\t', phase_rewrapped(1:end-1));
fprintf('%.1f', phase_rewrapped(end));
fprintf('\n');
phase_rewrapped = angle(exp(j*(phase_ave(:).' - param.phase)/180*pi)) * 180/pi;
fprintf('%.1f\t', phase_rewrapped(1:end-1));
fprintf('%.1f', phase_rewrapped(end));
fprintf('\n');

return;



