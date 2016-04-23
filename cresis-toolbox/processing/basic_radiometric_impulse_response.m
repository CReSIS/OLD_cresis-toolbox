% script basic_radiometric_impulse_response
%
% This script evaluates the signal level from a perfect reflector. It plots
% the raw impulse response, pulse compressed response, and the raw
% spectrum.
%
% 1. Collect data with deconv mode over flat water.
%
% Author: John Paden

%% Load the files

% Load a single file for evaluation
original_img = param.img;
param.img = param.img(1,:);
basic_file_loader;
param.img = original_img;

%% basic_radiometric_impulse_response preparation
dt = 1/fs;
Nt = size(data,1);
time = dt*(0:Nt-1);
clear pc_param;
pc_param.DDC_mode = DDC_mode;
pc_param.DDC_freq = DDC_freq;
pc_param.f0 = f0;
pc_param.f1 = f1;
pc_param.Tpd = Tpd;
pc_param.zero_pad = 1;
pc_param.decimate = true;
pc_param.window_func = @hanning;
pc_param.time = t0 + (0:dt:(Nt-1)*dt).';
pc_param.tukey = tukey;
[pc_signal,pc_time] = pulse_compress(data,pc_param);

%% Plot the echogram
imagesc(lp(pc_signal));

%% Elevation compensation
wf_adc = 1;
try
  elev_filt_len = round(length(elev(wf_adc,:))/20)*2+1;
  elev(wf_adc,:) = sgolayfilt(elev(wf_adc,:), 2, elev_filt_len, hanning(elev_filt_len));
end

wf = abs(param.img(wf_adc,1));
adc = abs(param.img(wf_adc,2));

% Create frequency axis
dt = pc_time(2) - pc_time(1);
Nt = length(pc_time);
T = Nt*dt;
df = 1/T;
freq = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';

% Correct all the data to a constant elevation (no zero padding is
% applied so wrap around could be an issue)
for rline = 1:size(pc_signal,2)
  elev_dt = (elev(rline) - elev(1)) / (c/2);
  pc_signal(:,rline,wf_adc) = ifft(fft(pc_signal(:,rline,wf_adc)) .* exp(1i*2*pi*freq*elev_dt));
end

[max_value,surf_bin] = max(pc_signal);

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
fprintf('Best rangeline: %d with peakiness: %.1f dB\n', best_idx, best_val);
try
  user_best_idx = input(sprintf('Enter range line [%d]: ',best_idx));
  if ~isempty(user_best_idx)
    best_idx = user_best_idx;
  end
end

param.recs = [best_idx param.analysis.specular.ave];

param.file_search_mode = 'default';
basic_file_loader;

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
    freq = fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).';
    
    % Correct all the data to a constant elevation (no zero padding is
    % applied so wrap around could be an issue)
    for rline = 1:size(data,2)
      elev_dt = (elev(rline) - elev(1)) / (c/2);
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
  dt = 1/fs;
  good_idxs = find(pc_param.time<pc_time(surf_bin(best_idx)),1,'last') + (0:ceil(Tpd/dt));
  h_fig = figure(wf_adc); clf;
  set(h_fig,'WindowStyle','Docked');
  set(h_fig,'NumberTitle','off');
  set(h_fig,'Name',sprintf('WF %d ADC %d',wf,adc));
  rel_time = pc_param.time(good_idxs);
  rel_time = rel_time - rel_time(1);
  
  plot(rel_time*1e6, real(mean(data(good_idxs,:,wf_adc), 2)) )
  
  if 0
    hold on
    plot(rel_time*1e6, imag(data(good_idxs,1,wf_adc)))
  elseif 0
    hold on
    plot(rel_time*1e6, real(data(good_idxs,end,wf_adc)))
  end
  
  title(sprintf('WF %d ADC %d',wf,adc));
  xlabel('Time (us)');
  ylabel('Voltage (V)');
  grid on;
  xlim(rel_time([1 end])*1e6)
  ylim([ min(min(real(data(good_idxs,1,:))))  max(max(real(data(good_idxs,1,:)))) ]);
  if 1
    saveas(h_fig,sprintf('radiometric_wf_%d_adc_%d.fig',wf,adc));
  end
  
  signal_power_dBm(wf_adc) = 10*log10((max(abs(data(good_idxs,1,wf_adc)))/sqrt(2))^2/50)+30;
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
  dt = 1/fs;
  good_idxs = find(pc_param.time<pc_time(surf_bin(best_idx)),1,'last') + (0:ceil(Tpd/dt));
  h_fig = figure(100+wf_adc); clf;
  set(h_fig,'WindowStyle','Docked');
  set(h_fig,'NumberTitle','off');
  set(h_fig,'Name',sprintf('WF %d ADC %d',wf,adc));
  
  % Create frequency axis
  dt = 1/fs;
  Nt = length(good_idxs);
  T = Nt*dt;
  df = 1/T;
  freq = fftshift(fc + ifftshift( -floor(Nt/2)*df : df : floor((Nt-1)/2)*df ).');
  
  plot(freq/1e6, signal_psd(:,1,wf_adc));
  
  title(sprintf('WF %d ADC %d',wf,adc));
  xlabel('Frequency (MHz)');
  ylabel('Relative power (dB)');
  grid on;
  xlim((fc + 1.2*[-BW/2 BW/2])/1e6);
  ylim(ylims);
  if 1
    saveas(h_fig,sprintf('signal_psd_wf_%d_adc_%d.fig',wf,adc));
  end
end

%% Done
return;

