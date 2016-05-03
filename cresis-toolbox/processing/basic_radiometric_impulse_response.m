function basic_radiometric_impulse_response(param,defaults)
% basic_radiometric_impulse_response(param,defaults)
%
% This script evaluates the signal level from a perfect reflector. It plots
% the raw impulse response, pulse compressed response, and the raw
% spectrum.
%
% 1. Collect data with deconv mode over flat water.
%
% Author: John Paden

physical_constants;

%% Load the files

% Load a single file for evaluation
original_img = param.img;
param.img = param.img(1,:);
[data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults);
param.img = original_img;

%% basic_radiometric_impulse_response preparation
[pc_signal,pc_time] = pulse_compress(data,pc_param);

%% Plot the echogram and track surface
good_time_bins = find(pc_time > pc_param.Tpd*1.1);
[max_value,surf_bin] = max(pc_signal(good_time_bins,:));
surf_bin = surf_bin + good_time_bins(1)-1;

done = false;
while ~done
  figure(1000); clf;
  imagesc(lp(pc_signal));
  hold on
  plot(good_time_bins([1 1 end end 1]), [1 size(pc_signal,2) size(pc_signal,2) 1 1]);
  hold off
  try
    min_range_bin = input(sprintf('\nEnter min range bin for surface [%d]:', good_time_bins(1)));
    if isempty(min_range_bin)
      done = true;
    else
      good_time_bins = min_range_bin(1):size(pc_signal,1);
      done = true;
    end
  end
end

%% Elevation compensation
wf_adc = 1;
elev_filt_len = round(length(hdr.elev(wf_adc,:))/20)*2+1;
hdr.elev(wf_adc,:) = sgolayfilt(hdr.elev(wf_adc,:), 2, elev_filt_len, hanning(elev_filt_len));

wf = abs(param.img(wf_adc,1));
adc = abs(param.img(wf_adc,2));

% Create frequency axis
dt = pc_time(2) - pc_time(1);
Nt = length(pc_time);
T = Nt*dt;
df = 1/T;
freq = pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

% Correct all the data to a constant elevation (no zero padding is
% applied so wrap around could be an issue)
for rline = 1:size(pc_signal,2)
  elev_dt = (hdr.elev(rline) - hdr.elev(1)) / (c/2);
  pc_signal(:,rline,wf_adc) = ifft(fft(pc_signal(:,rline,wf_adc)) .* exp(1i*2*pi*freq*elev_dt));
end

% Perform STFT (short time Fourier transform) (i.e. overlapping short FFTs in slow-time)
param.analysis.specular.ave = 128;
param.analysis.specular.signal_doppler_bins = [1:4,125:128];
param.analysis.specular.noise_doppler_bins = [8:121];
H = spectrogram(double(max_value),hanning(param.analysis.specular.ave),param.analysis.specular.ave/2,param.analysis.specular.ave);

% Since there may be a little slope in the ice, we sum the powers from
% the lower frequency doppler bins rather than just taking DC. It seems to help
% a lot to normalize by the sum of the middle/high-frequency Doppler bins.   A coherent/specular
% surface will have high power in the low bins and low power in the high bins
% so this ratio makes sense.
peakiness = lp(max(abs(H(param.analysis.specular.signal_doppler_bins,:)).^2) ./ mean(abs(H(param.analysis.specular.noise_doppler_bins,:)).^2));

[best_val,best_idx] = max(peakiness);
best_idx = round((best_idx-1) * param.analysis.specular.ave/2);
fprintf('\nBest rangeline: %d with peakiness: %.1f dB\n', best_idx, best_val);
try
  user_best_idx = input(sprintf('Enter range line [%d]: ',best_idx));
  if ~isempty(user_best_idx)
    best_idx = user_best_idx;
  end
end

param.recs = [best_idx param.analysis.specular.ave];

param.file_search_mode = 'default';
[data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults);

%% Convert from quantization to voltage @ ADC
data = data ...
  * default.radar.adc_full_scale/2^default.radar.adc_bits ...
  * 2^hdr.wfs(wf).bit_shifts / hdr.wfs(wf).presums;

%% Additional software presums
for wf_adc = 1:size(data,3)
  data(:,:,wf_adc) = fir_dec(data(:,:,wf_adc),param.presums);
end

%% Elevation Compensation
% =====================================================================

try
  for wf_adc = 1:size(param.img,1)
    % Create frequency axis
    dt = 1/fs;
    Nt = size(data,1);
    T = Nt*dt;
    df = 1/T;
    freq = pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    
    % Correct all the data to a constant elevation (no zero padding is
    % applied so wrap around could be an issue)
    for rline = 1:size(data,2)
      elev_dt = (hdr.elev(rline) - hdr.elev(1)) / (c/2);
      data(:,rline,wf_adc) = ifft(fft(data(:,rline,wf_adc)) .* exp(1i*2*pi*freq*elev_dt));
    end
  end
end

%% Waveform power
% =====================================================================

signal_power_dBm = [];
expected_power_received_dBm = [];
for wf_adc = 1:size(param.img,1)
  wf = abs(param.img(wf_adc,1));
  adc = abs(param.img(wf_adc,2));
  dt = 1/hdr.fs;
  Nt = size(data,1);
  h_fig = figure(wf_adc); clf;
  set(h_fig,'WindowStyle','Docked');
  set(h_fig,'NumberTitle','off');
  set(h_fig,'Name',sprintf('WF %d ADC %d',wf,adc));

  if 1
    % Determine the phase and delay offset of each range line
    [pc_signal,pc_time] = pulse_compress(data(:,:,wf_adc),pc_param);
    
    max_value = zeros(1,size(pc_signal,2));
    max_idx_unfilt = zeros(1,size(pc_signal,2));
    Mt = 100;
    for idx = 1:size(pc_signal,2)
      oversampled_rline = interpft(pc_signal(:,idx),size(pc_signal,1)*Mt);
      start_bin = surf_bin(best_idx)*Mt-1000;
      [max_value(idx),max_idx_unfilt(idx)] ...
        = max(oversampled_rline(start_bin:end));
      max_idx_unfilt(idx) = max_idx_unfilt(idx) + start_bin - 1;
    end
    
    % Filter the delay (max_idx) and phase (max_value) vectors
    max_idx = sgolayfilt(max_idx_unfilt/Mt,3,17);
    max_idx = max_idx - mean(max_idx);
    phase_corr = sgolayfilt(double(unwrap(angle(max_value))),3,17);
    
    % Compensate range lines for phase and delay variance
    % in the peak value
    T = Nt*dt;
    df = 1/T;
    fc = (pc_param.f0+pc_param.f1)/2;
    freq = pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    
    % Apply true time delay shift to flatten surface
    comp_data = ifft(fft(data(:,:,wf_adc)) .* exp(1i*2*pi*freq*max_idx/Mt*dt) );
    % Apply phase correction (compensating for phase from time delay shift)
    comp_data = comp_data .* repmat(exp(-1i*(phase_corr + 2*pi*fc*max_idx/Mt*dt)), [Nt 1]);
  
  else
    % Do not do phase and delay corrections
    comp_data = data(:,:,wf_adc);
  end
  
  % Plot over-interpolated results
  Mt = 10;
  time_Mt = pc_param.time(1) + dt/Mt * (0:Nt*Mt-1).';
  data_Mt = interpft(comp_data, Nt*Mt);

  good_idxs_Mt = find(time_Mt<pc_time(surf_bin(best_idx)),1,'last') + (0:ceil(pc_param.Tpd/dt*Mt));
  rel_time = time_Mt(good_idxs_Mt);
  rel_time = rel_time - rel_time(1);

  plot(rel_time*1e6, real(mean(data_Mt(good_idxs_Mt,:), 2)) )
  
  if 1
    hold on
    plot(rel_time*1e6, real(data_Mt(good_idxs_Mt,1)))
  elseif 0
    hold on
    plot(rel_time*1e6, imag(data_Mt(good_idxs,1)))
  end
    
  good_idxs = find(pc_param.time<pc_time(surf_bin(best_idx)),1,'last') + (0:ceil(pc_param.Tpd/dt));
  title(sprintf('WF %d ADC %d',wf,adc));
  xlabel('Time (us)');
  ylabel('Voltage (V)');
  grid on;
  xlim(rel_time([1 end])*1e6)
  ylim([ 1.1*min(min(real(data(good_idxs,1,:))))  1.1*max(max(real(data(good_idxs,1,:)))) ]);
  if 0
    saveas(h_fig,sprintf('radiometric_wf_%d_adc_%d.fig',wf,adc));
  end
  
  signal_power_dBm(wf_adc) = 10*log10((max(abs(data_Mt(good_idxs_Mt,1)))/sqrt(2))^2/50)+30;
  range = pc_time(surf_bin(best_idx))*c/2;
  power_reflectance = 1;
  expected_power_received_dBm(wf_adc) = 10*log10(default.Pt * default.Gt * default.Ae / (4*pi*(2*range).^2) * default.system_loss_dB * power_reflectance) + 30;
end
fprintf('Expected and measured receive power in dBm\n');
fprintf('ADC     \t'); fprintf('%d\t', param.img(:,2)); fprintf('\n');
fprintf('Expected\t'); fprintf('%.0f\t', expected_power_received_dBm); fprintf('\n');
fprintf('Measured\t'); fprintf('%.0f\t', signal_power_dBm); fprintf('\n');

%% Waveform PSD
% =====================================================================

signal_psd = fftshift(lp(mean(abs(fft(data(good_idxs,:,:))).^2,2),1), 1);
ylims = [min(signal_psd(:))-1 max(signal_psd(:))+1];
ylims = max(signal_psd(:)) + [-25 0];

signal_power_dBm = [];
expected_power_received_dBm = [];
for wf_adc = 1:size(param.img,1)
  wf = abs(param.img(wf_adc,1));
  adc = abs(param.img(wf_adc,2));
  dt = 1/hdr.fs;
  good_idxs = find(pc_param.time<pc_time(surf_bin(best_idx)),1,'last') + (0:ceil(pc_param.Tpd/dt));
  h_fig = figure(100+wf_adc); clf;
  set(h_fig,'WindowStyle','Docked');
  set(h_fig,'NumberTitle','off');
  set(h_fig,'Name',sprintf('WF %d ADC %d',wf,adc));
  
  % Create frequency axis
  dt = 1/hdr.fs;
  Nt = length(good_idxs);
  T = Nt*dt;
  df = 1/T;
  freq = fftshift(pc_param.DDC_freq + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).');
  
  plot(freq/1e6, signal_psd(:,1,wf_adc));
  
  title(sprintf('WF %d ADC %d',wf,adc));
  xlabel('Frequency (MHz)');
  ylabel('Relative power (dB)');
  grid on;
  xlim((pc_param.DDC_freq + 1.2*[-hdr.BW/2 hdr.BW/2])/1e6);
  ylim(ylims);
  if 0
    saveas(h_fig,sprintf('signal_psd_wf_%d_adc_%d.fig',wf,adc));
  end
end

%% Done
return;

