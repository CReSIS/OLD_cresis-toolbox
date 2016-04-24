function basic_noise_analysis(param,defaults)
% basic_noise_analysis(param,defaults)
%
% Use run_basic_noise_analysis to run.
%
% This script evaluates 50 ohm term and receive only data.
% It loads all receive channels and then analyzes these data.
% Normally called from run_basic_noise_analysis.
%
% 1. Collect data with ? waveforms in whatever noise configuration
%    you want to measure.
% 2. Characterization should be done for 50 ohm and receive only at
%    least.
% 3. If transmitting, time gate should be large enough to include
%    noise-only data at the end of the record.
%
% Author: John Paden

physical_constants;

%% Load the files
[data,fn,settings,default,gps,hdr,pc_param] = basic_file_loader(param,defaults);

%% basic_noise_analysis preparation
[fn_dir fn_name] = fileparts(fn);
if ~isfield(param,'seg') || isempty(param.seg)
  param.seg = -1;
end
clear data_tmp;

default_noise_rbins = 1:size(data,1);
if ~isfield(param,'noise_rbins') || isempty(param.noise_rbins)
  noise_rbins = default_noise_rbins;
elseif param.noise_rbins(end) > size(data,1)
  noise_rbins = param.noise_rbins(1):size(data,1);
elseif param.noise_rbins(1) <= 0
  noise_rbins = size(data,1)+param.noise_rbins;
else
  noise_rbins = param.noise_rbins(1):param.noise_rbins(end);
end
noise_rbins = intersect(noise_rbins,default_noise_rbins);

wf = abs(param.img(1,1));

%% Convert from quantization to voltage @ ADC
data = data ...
  * default.radar.adc_full_scale/2^default.radar.adc_bits ...
  * 2^hdr.wfs(abs(wf)).bit_shifts / hdr.wfs(wf).presums;

%% Additional software presums
for wf_adc = 1:size(data,3)
  data(:,:,wf_adc) = fir_dec(data(:,:,wf_adc),param.presums);
end

%% Noise power
% =====================================================================

% Calculate noise power as Vrms (assuming param.presums = 1)
if strcmp(param.radar_name,'mcords5') && isfield(hdr,'DDC') && hdr.DDC(1) >= 2
  % Add 3 dB for IQ combination
  fprintf('Expected ADC noise floor @ ADC %.1f dBm\n', lp((default.radar.adc_full_scale/2/sqrt(2))^2/50)+30 - default.radar.adc_SNR_dB + 3 );
  fprintf('Expected Rx noise floor @ ADC %.1f dBm\n', lp(BoltzmannConst*290*hdr.BW_noise*default.radar.noise_figure*hdr.rx_gain^2) + 3 +30);
else
  fprintf('Expected ADC noise floor @ ADC %.1f dBm\n', lp((default.radar.adc_full_scale/2/sqrt(2))^2/50)+30 - default.radar.adc_SNR_dB );
  fprintf('Expected Rx noise floor @ ADC %.1f dBm\n', lp(BoltzmannConst*290*hdr.BW_noise*default.radar.noise_figure*hdr.rx_gain^2)+30);
end
fprintf('Expected levels only valid for param.presums = 1\n');
fprintf('Noise power (dBm) at each ADC rx input and relative to 50 ohm (dB):\n')
noise_power_dBm = permute(lp(mean(mean(abs(data(noise_rbins,:,:)).^2/50, 1), 2) ...
  / 2^hdr.wfs(abs(wf)).bit_shifts * hdr.wfs(wf).presums * param.presums, 1) + 30, [1 3 2]);
fprintf('ADC\t'); fprintf('%.0f\t', 1:24); fprintf('\n');
fprintf('Noise\t'); fprintf('%+.1f\t', noise_power_dBm); fprintf('\n');
fprintf('Rel\t'); fprintf('%+.1f\t', noise_power_dBm - default.noise_50ohm); fprintf('\n');

%% Quantization analysis
% =====================================================================
if param.pdf_en
  for adc_idx = 1:size(data,3)
    figure(adc_idx); clf;
    plot(real(data(:,1,adc_idx)),'.');
    grid on;
    xlabel('Range bin');
    ylabel('Quantization level');
    
    figure(100+adc_idx);
    imagesc(lp(data(:,:,adc_idx),2));
    colorbar;
    
    % Plot estimated pdf and approximate a gaussian to it
    figure(200+adc_idx);
    ROI = data(noise_rbins,:,adc_idx);
    ROI = real(ROI(:));
    [n,x] = hist(ROI,64);    
    bar(x,n);
    mean_x = mean(ROI);
    var_x = var(ROI);
    hold on;
    plot(x, numel(ROI)*(x(2)-x(1)) * 1/sqrt(2*pi*var_x) * exp(-(x - mean_x).^2 / (2*var_x)),'r');
    hold off;
  end
  for adc_idx = 1:size(data,3)
    set(adc_idx,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Q%d',adc_idx));
  end  
  for adc_idx = 1:size(data,3)
    set(100+adc_idx,'WindowStyle','docked','NumberTitle','off','Name',sprintf('E%d',adc_idx));
  end  
  for adc_idx = 1:size(data,3)
    set(200+adc_idx,'WindowStyle','docked','NumberTitle','off','Name',sprintf('P%d',adc_idx));
  end  
end

%% Power Spectrum
% =====================================================================
if param.psd_en
  plot_combined_psd = true;
  combined_psd_cmap = hsv(size(data,3));
  combined_psd_legend = {};
  if plot_combined_psd
    h_psd_fig = figure(500); clf; h_psd_axes = axes('parent',h_psd_fig); hold(h_psd_axes,'on'); grid(h_psd_axes,'on'); xlabel('Frequency (MHz)','parent',h_psd_axes); ylabel('Relative noise power (dB)','parent',h_psd_axes);
  end
  
  for adc_idx = 1:size(data,3)
    adc = param.img(adc_idx,2);
    
    fir_data = fir_dec(data(noise_rbins(1):noise_rbins(end),:,adc_idx),param.presums);
    
    if strcmp(param.radar_name,'mcords5') && isfield(hdr,'DDC') && hdr.DDC(1) >= 2
      dt = pc_param.time(2) - pc_param.time(1);
      Nt = size(fir_data,1);
      df = 1/(Nt*dt);
      freq = pc_param.DDC_freq + (-floor(Nt/2)*df : df : floor((Nt-1)/2)*df);
      
      figure(400+adc); clf;
      set(400+adc,'WindowStyle','docked','NumberTitle','off','Name',sprintf('M%d',adc_idx));
      plot(freq/1e6, lp(mean(abs(fftshift(fft(fir_data))).^2*2^2 / 50,2)/size(fir_data,1)) + 30)
      title(sprintf('MeanFFT adc%d ave%d %s/%s', adc, param.presums, param.seg, fn_name),'Interpreter','none');
      ylabel('Relative power (dB)');
      xlabel('Frequency (MHz)');
      grid on;
      ylims = ylim;
      
      figure(300+adc); clf;
      set(300+adc,'WindowStyle','docked','NumberTitle','off','Name',sprintf('FFT%d',adc_idx));
      imagesc([], freq/1e6, lp(fftshift(fft(fir_data))) + 30 + 10*log10(2^2/50/size(fir_data,1)) )
      title(sprintf('Freq-space adc%d ave%d %s/%s', adc, param.presums, param.seg, fn_name),'Interpreter','none');
      xlabel('Range line');
      ylabel('Frequency (MHz)');
      h = colorbar;
      set(get(h,'YLabel'),'String','Relative power (dB)');
      caxis_lims = caxis;
      caxis([ylims(1) caxis_lims(2)]);
      
      if plot_combined_psd
        plot(freq/1e6, lp(mean(abs(fftshift(fft(fir_data))).^2*2^2 / 50,2)/size(fir_data,1)) + 30, 'parent',h_psd_axes,'color',combined_psd_cmap(adc_idx,:))
        combined_psd_legend{adc_idx} = sprintf('chan %d', adc);
      end
      
    else
      pc_param.time = hdr.wfs(wf).t0 + (0:size(fir_data,1)-1)/fs;
      dt = pc_param.time(2) - pc_param.time(1);
      Nt = length(pc_param.time);
      df = 1/(Nt*dt);
      freq = (0:df:(Nt-1)*df).';
      
      figure(300+adc); clf;
      set(300+adc,'WindowStyle','docked','NumberTitle','off','Name',sprintf('FFT%d',adc_idx));
      imagesc([], freq/1e6, lp(fft(fir_data)) + 30 + 10*log10(2^2/50/size(fir_data,1)) )
      title(sprintf('Freq-space adc%d ave%d %s/%s', adc, param.presums, param.seg, fn_name),'Interpreter','none');
      xlabel('Range line');
      ylabel('Frequency (MHz)');
      ylim(fs/1e6*[0 1]);
      h = colorbar;
      set(get(h,'YLabel'),'String','Relative power (dB)');
      
      figure(400+adc); clf;
      set(400+adc,'WindowStyle','docked','NumberTitle','off','Name',sprintf('M%d',adc_idx));
      plot(freq/1e6, lp(mean(abs(fft(fir_data)).^2*2^2 / 50,2)/size(fir_data,1)) + 30)
      title(sprintf('MeanFFT adc%d ave%d %s/%s', adc, param.presums, param.seg, fn_name),'Interpreter','none');
      ylabel('Relative power (dB)');
      xlabel('Frequency (MHz)');
      xlim(fs/1e6*[0 1]);
      grid on;
    end
  end
end
figure(500);
legend(combined_psd_legend,'location','best')

%% Done
return;

