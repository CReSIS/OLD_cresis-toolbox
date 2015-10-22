% script radiometric_cal_2012_Antarctica_DC8
%
% Take tx equalization data and plot the Fresnel reflection coefficient.
% Do this for one file, single waveform, all receivers?
%
% Author: John Paden

physical_constants;
close all;
tstart = tic;

% =======================================================================
% User Settings
% =======================================================================

param.fs = 1e9/9;
param.rlines = [];
param.noise_rbins = [4000 5000];
param.wf = 2;
param.radar_name = 'mcords2';

% .gps_fn = Optional GPS file name (leave empty to disable)
param.gps_fn = '/home/polargrid/csarp_support/gps/2010_Antarctica_DC8_GPS/20101026_ALL_pos.mat';
param.gps_fn = '';

% .presums = Number of presums (coherent averaging) to do
param.presums = 10;

% .noise_removal_en = mcords noise removal (i.e. should be false for
%    mcrds, mcords2, etc)
param.noise_removal_en = false;

param.pc_param.f0 = 189.15e6;
param.pc_param.f1 = 198.65e6;
% param.pc_param.f0 = 180e6;
% param.pc_param.f1 = 210e6;
param.pc_param.Tpd = 10e-6;
param.pc_param.tukey = 0.2;

% =======================================================================
% User Settings
% =======================================================================

% ======================================================================
%    (THIS NEEDS TO BE SET EVERYTIME)
% ======================================================================
% .adcs = the receive channels to use
param.adcs = [1 2 3 4];

load_param = [];
idx = 1;
load_param(idx).base_path = 'B:\';
load_param(idx).acquisition_num = 0; 
load_param(idx).seg = ''; % 0,0 settings w/ 0.7 dB loss post-filter
load_param(idx).file_num = 12;
% load_param(idx).seg = 'seg_00';
% load_param(idx).file_num = 3;
% load_param(idx).seg = 'seg_14';
% load_param(idx).file_num = 2;
load_param(idx).legend = 'Ocean, level, 20000 ft';
load_param(idx).marker = 'k';

% 21:35:50 UTC Alluvial Deposits (first site on testflight) seg_07
% 21:53:42 UTC Dry lake site (second site) seg_08
% 22:10:59 UTC Chelsa (sp?) Dunes site (third site): seg_09
% 22:23:22 UTC third site: seg_10
% 22:41:12 UTC second site: seg_11
% 22:57:07 UTC first site: seg_12 (windmills turning

% .presums = Number of presums (coherent averaging) to do
param.presums = 1;

adc_bits = 14;
Vpp_scale = 2.23;
adc_SNR_dB = 72;
rx_gain = 10^((50.25-0.7)/20);
noise_figure = 10^(1.96/10); % Do not include receiver losses
% bandwidth_range = [188e6 202e6]; % Leave empty to use all the bandwidth (~30 MHz)
% bandwidth = 10e6;
bandwidth_range = []; % Leave empty to use all the bandwidth (~30 MHz)
bandwidth = 30e6;

rline = 1;

param.pc_param.f0 = 189.15e6;
param.pc_param.f1 = 198.65e6;
% param.pc_param.f0 = 180e6;
% param.pc_param.f1 = 210e6;
param.pc_param.Tpd = 1e-6;
param.pc_param.tukey = 0.2;

% =======================================================================
% =======================================================================
% Automated Section
% =======================================================================
% =======================================================================

% =======================================================================
% Load data
% =======================================================================
fprintf('Vpp Scale %f Vpp\n', Vpp_scale);
fprintf('ADC SNR %f dB\n', adc_SNR_dB);
fprintf('Noise figure %f dB\n', 10*log10(noise_figure));
fprintf('Rx gain %f dB\n', 10*log10(rx_gain));
fprintf('Bandwidth %f MHz\n', bandwidth/1e6);
fprintf('Expected ADC noise floor @ ADC %.1f dBm\n', lp((Vpp_scale/2/sqrt(2))^2/50)+30 - adc_SNR_dB );
fprintf('Expected Rx noise floor @ ADC %.1f dBm\n', lp(BoltzmannConst*290*bandwidth*noise_figure*rx_gain^2)+30);
fprintf('Expected levels only valid for param.presums == 1, param.presums is %i\n', param.presums);
fprintf('Noise power is in dBm:\n')

legend_label = {};
for load_idx = 1:length(load_param)
  param.base_path = load_param(load_idx).base_path;
  param.acquisition_num = load_param(load_idx).acquisition_num;
  param.seg = load_param(load_idx).seg;
  param.file_num = load_param(load_idx).file_num;
  [hdr,data,fn_name] = noise_meas_2012_Antarctica_DC8_func(param);
  
  if ~isfield(param,'rlines') || isempty(param.rlines)
    rlines = 1:size(data,2);
  elseif param.rlines(end) > size(data,2)
    rlines = param.rlines(1):size(data,2);
  else
    rlines = param.rlines(1):param.rlines(end);
  end
  if ~isfield(param,'noise_rbins') || isempty(param.noise_rbins)
    noise_rbins = 1:size(data,1);
  elseif param.noise_rbins(end) > size(data,1)
    noise_rbins = param.noise_rbins(1):size(data,1);
  else
    noise_rbins = param.noise_rbins(1):param.noise_rbins(end);
  end
  
  % =====================================================================
  % Select only the desired data
  % Additional software presums
%   fir_data = [];
%   for adc_idx = 1:size(data,3)
%     if isempty(bandwidth_range)
%       % Use all the bandwidth
%       fir_data(:,:,adc_idx) = fir_dec(data(noise_rbins,rlines,adc_idx),param.presums);
%     else
%       [B,A] = butter(2,sort(abs(bandwidth_range-2*param.fs))/(param.fs/2));
%       fir_data(:,:,adc_idx) = fir_dec(filtfilt(B,A,double(data(noise_rbins,rlines,adc_idx))),param.presums);
%     end
%   end
  fir_data = data;
  
  % =======================================================================
  % Convert from quantization to voltage @ ADC
  fir_data = fir_data ...
    * Vpp_scale/2^adc_bits ...
    * 2^hdr.wfs(param.wf).bit_shifts/hdr.wfs(param.wf).presums;

  fprintf('Pulse compression (%.1f sec)\n', toc(tstart));
clear pc_param pc_time pc_data;
for wf = 1:size(fir_data,3)
  wf
  pc_param(wf).f0 = param.pc_param.f0;
  pc_param(wf).f1 = param.pc_param.f1;
  pc_param(wf).Tpd = param.pc_param.Tpd;
  pc_param(wf).tukey = param.pc_param.tukey;
  pc_param(wf).time = hdr.wfs(param.wf).t0 + (0:size(fir_data,1)-1)/fs;
  [pc_data(:,:,wf),pc_time] = pulse_compress(fir_data(:,:,wf),pc_param(wf));
end
imagesc(lp(pc_data(:,:,1)))
return;

  % =====================================================================
  % Noise power
  % =====================================================================
  
  dt = 1/param.fs;
  Nt = size(data,1);
  time = dt*(0:(Nt-1)).';
  % Calculate noise power as Vrms (assuming param.presums = 1)
  fprintf('%s\n  ', load_param(load_idx).legend);
  for adc_idx = 1:size(fir_data,3)
    if adc_idx > 1
      fprintf('\t');
    end
    fprintf('%.1f ', ...
      lp(mean(mean(abs(fir_data(:,:,adc_idx)).^2/50)) * hdr.wfs(param.wf).presums,1)+30 );
  end
  fprintf('\n');
  
  
  % =====================================================================
  % Power Spectrum
  % =====================================================================
  for adc_idx = 1:size(fir_data,3)
    adc = param.adcs(adc_idx);
    
    clear pc_param;
    pc_param.time = hdr.wfs(param.wf).t0 + (0:size(fir_data,1)-1)/param.fs;
    dt = pc_param.time(2) - pc_param.time(1);
    Nt = length(pc_param.time);
    df = 1/(Nt*dt);
    freq = param.fs + (0:df:(Nt-1)*df).';
    
    figure(120+adc);
    if load_idx == 1
      clf;
    else
      hold on;
    end
    %     set(120+adc,'WindowStyle','docked','NumberTitle','off','Name',sprintf('M%d',adc_idx));
    plot(freq/1e6, lp(mean(abs(fft(fir_data(:,:,adc_idx))).^2*2^2 / 50,2)/size(fir_data,1)) + 30, load_param(load_idx).marker)
    title(sprintf('MeanFFT adc%d ave%d %s/%s', adc, param.presums, param.seg, fn_name),'Interpreter','none');
    ylabel('Relative power (dB)');
    xlabel('Frequency (MHz)');
    xlim(param.fs/1e6*[1.5 2]);
    grid on;
  end
  legend_label{load_idx} = load_param(load_idx).legend;
end

for adc_idx = 1:size(fir_data,3)
  adc = param.adcs(adc_idx);
  figure(120+adc);
  legend(legend_label,'Location','best');
  ylim([-59 -52])
end

return;
