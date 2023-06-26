% Script run_basic_load
%
% Several examples of using basic_load_mcords and basic_load_mcords2.
% Set "run_example" before runng the script.
%
% Author: John Paden
%
% See also basic_load_mcords, basic_load_mcords2, pulse_compress

physical_constants;
close all;
tstart = tic;

% =======================================================================
% User Settings: Select Example
% =======================================================================
run_example = 2;

% =======================================================================
% User Settings: To Load Data
% =======================================================================
fs = 1e9/9;
param.radar_name = 'mcords';
if strcmpi(param.radar_name,'mcords')
  param.radar_num = 4;
  % .adc = the adc channel to use
  param.adc = 1;
  
  % Parameters to locate specific file of interest
  %    (THIS NEEDS TO BE SET EVERYTIME)
  param.data_file_num = 1;
  param.base_path = '/home/polargrid/mcords/2010_Antarctica_DC8/20101106/seg3/';
elseif strcmpi(param.radar_name,'mcords2')
  % .adc = the receive channel to use (relative to the board # so that it
  %    is always contained in [1,4] since there are 4 channels per board)
  param.adc = 1;
  param.acquisition_num = [];
  param.board = 3;
  
  % Parameters to locate specific file of interest
  % (THIS NEEDS TO BE SET EVERYTIME)
  param.seg = 25;
  param.file_num = 19;
  param.base_path = 'D:\mcords2\2011_Greenland_P3\20110411\';
  param.base_path = '/cresis/data4/MCoRDS/2011_Greenland_P3/20110314/';
end

% =======================================================================
% Load data
% =======================================================================
fprintf('============================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(tstart));

clear data;
if strcmpi(param.radar_name,'mcords')
  file_midfix = sprintf('r%d-%d.',param.radar_num,param.adc);
  file_suffix = sprintf('.%04d.dat',param.data_file_num);
  fprintf('  Path: %s\n', param.base_path);
  fprintf('  Match: mcords*%s*%s\n', file_midfix, file_suffix);
  fn = get_filename(param.base_path,'mcords',file_midfix,file_suffix);
  if isempty(fn)
    fprintf('  Could not find any files which match\n');
    return;
  end
  [hdr,data] = basic_load_mcords(fn, struct('clk',fs,'first_byte',2^26));

  fn = get_filename(fn_dir, file_prefix, '', file_suffix);
  if isempty(fn)
    fprintf('  Could not find any files which match\n');
    return;
  end
  fprintf('  Loading file %s\n', fn);
  [hdr,data] = basic_load_mcords2(fn,struct('clk',fs));
  basic_remove_mcords_digital_errors;

elseif strcmpi(param.radar_name,'mcords2')
  fn_dir = fullfile(param.base_path, sprintf('board%d',param.board), ...
    sprintf('seg_%02d',param.seg));
  file_prefix = sprintf('mcords2_%d_',param.board);
  if isempty(param.acquisition_num)
    file_suffix = sprintf('%04d.bin',param.file_num);
  else
    file_suffix = sprintf('%02d_%04d.bin',param.acquisition_num,param.file_num);
  end
  fprintf('  Path: %s\n', fn_dir);
  fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
  fn = get_filename(fn_dir, file_prefix, '', file_suffix);
  if isempty(fn)
    fprintf('  Could not find any files which match\n');
    return;
  end
  fprintf('  Loading file %s\n', fn);
  [hdr,data] = basic_load_mcords2(fn,struct('clk',fs));
end

[fn_dir fn_name] = fileparts(fn);

if run_example == 1
  % =======================================================================
  % =======================================================================
  % Setup loading parameters for example 1
  %  - Plot data in time-space, time-wavenumber, freq-space, and
  %    freq-wavenumber domains
  %  - Print out some basic statistics
  % =======================================================================
  % =======================================================================

  % =======================================================================
  % User Settings: Radar parameters and What to Plot
  % =======================================================================
  adc_bits = 14;
  Vpp_scale = 2;
  adc_SNR_dB = 70;
  rx_gain = 10^((72-20)/20);
  noise_figure = 10^(1.6/10); % Do not include receiver losses
  % Choose which ADC and waveform to plot and other settings
  plot_adc = 2;
  plot_wf = 2;
  noise_bins = 3000:4000;

  % =======================================================================
  % Convert from quantization to voltage @ ADC
  for wf = 1:length(data)
    data{wf} = (data{wf} - median(data{wf}(:,1))) ...
      * Vpp_scale/2^adc_bits ...
      * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums;
  end
  clear wf;
    
  % =======================================================================
  % Plot data
  clear pc_param;
  
  % Create time and freq axes
  pc_param.time = hdr.wfs(plot_wf).t0 + (0:size(data{plot_wf},1)-1)/fs;
  dt = pc_param.time(2) - pc_param.time(1);
  Nt = length(pc_param.time);
  df = 1/(Nt*dt);
  freq = (0:df:(Nt-1)*df).';
  
  figure(1); clf; set(1,'WindowStyle','docked','NumberTitle','off','Name','tx raw');
  imagesc([],[],lp( abs(data{plot_wf}(:,:,plot_adc)).^2/50 ) + 30);
  title(sprintf('Time-space raw wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  grid on;
  h = colorbar;
  set(get(h,'YLabel'),'String','Power (dBm)');

  % Calculate signal power assuming a constant modulus sinusoid
  Vpp = max(data{plot_wf}(:,1,plot_adc)) - min(data{plot_wf}(:,1,plot_adc));
  fprintf('Signal peak in record 1 @ ADC is %.1f mVpp %.1f dBm sinusoid\n', Vpp*1000, lp(abs(Vpp/2/sqrt(2)).^2/50)+30);
  
  % Calculate noise power as Vrms
  fprintf('Expected ADC noise floor @ ADC %.1f dBm\n', lp((Vpp_scale/2/sqrt(2))^2/50)+30 - adc_SNR_dB );
  fprintf('Expected Rx noise floor @ ADC %.1f dBm\n', lp(BoltzmannConst*290*30e6*noise_figure*rx_gain^2)+30);
  fprintf('Noise power is %.1f dBm\n', lp(mean(mean(abs(data{plot_wf}(noise_bins,:,plot_adc)).^2*hdr.wfs(plot_wf).presums/50)))+30 );
  
  figure(2); clf; set(2,'WindowStyle','docked','NumberTitle','off','Name','fx noise');
  imagesc([],freq/1e6,lp(fft(data{plot_wf}(noise_bins,:,plot_adc))));
  title(sprintf('Freq-space noise wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  ylabel('Frequency (MHz)');
  grid on;
  colorbar

  figure(3); clf; set(3,'WindowStyle','docked','NumberTitle','off','Name','tk');
  imagesc([],pc_param.time*1e6, ...
    lp(fft(data{plot_wf}(:,:,plot_adc),[],2)));
  title(sprintf('Time-k wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  grid on;
  colorbar

  figure(4); clf; set(4,'WindowStyle','docked','NumberTitle','off','Name','fk noise');
  imagesc([],freq/1e6,lp(fft2(data{plot_wf}(noise_bins,:,plot_adc))));
  title(sprintf('Freq-k noise wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  ylabel('Frequency (MHz)');
  grid on;
  colorbar

elseif run_example == 2
  % =======================================================================
  % =======================================================================
  % Setup loading parameters for example 2
  %  - Plot raw data in time-space, pulse compress externally with
  %    pulse_compress and plot this
  %  - Prints out expected power from specular surface
  % =======================================================================
  % =======================================================================
  
  
  % =====================================================================
  % User Settings: What to plot
  % =====================================================================
  adc_bits = 14;
  Vpp_scale = 2;
  adc_SNR_dB = 70;
  rx_gain = 10^((72-20)/20);
  noise_figure = 10^(1.6/10); % Do not include receiver losses
  % Choose which ADC and waveform to plot and other settings
  plot_adc = 1;
  plot_wf = 2;
  rline = 60;

  % =======================================================================
  % Convert from quantization to voltage @ ADC
  for wf = 1:length(data)
    data{wf} = (data{wf} - median(data{wf}(:,1))) ...
      * Vpp_scale/2^adc_bits ...
      * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums;
  end
  
  % =========================================================================
  % User Settings: Pulse Compression
  % =========================================================================  
  clear pc_param; wf = 0;
  if 0
    % Normal 2 waveform
    wf = wf + 1;
    pc_param(wf).f0 = 180e6;
    pc_param(wf).f1 = 210e6;
    pc_param(wf).Tpd = 1e-6;
    pc_param(wf).time = hdr.wfs(wf).t0 + (0:size(data{wf},1)-1).'/fs;
    pc_param(wf).tukey = 0.2;
    wf = wf + 1;
    pc_param(wf).f0 = 180e6;
    pc_param(wf).f1 = 210e6;
    pc_param(wf).Tpd = 10e-6;
    pc_param(wf).time = hdr.wfs(wf).t0 + (0:size(data{wf},1)-1).'/fs;
    pc_param(wf).tukey = 0.2;
  else
    % Tx Calibration
    for wf = 1:8
      pc_param(wf).f0 = 180e6;
      pc_param(wf).f1 = 210e6;
      pc_param(wf).Tpd = 10e-6;
      pc_param(wf).time = hdr.wfs(wf).t0 + (0:size(data{wf},1)-1).'/fs;
      pc_param(wf).tukey = 0.2;
    end
  end

  % =========================================================================
  % Pulse compress data
  % =========================================================================  
  clear data_pc time;
  for wf = 1:length(data)
    for adc = 1:size(data{wf},3)
      [data_pc{wf}(:,:,adc),time{wf}] = pulse_compress(data{wf}(:,:,adc),pc_param(wf));
    end
  end
  
  % =====================================================================
  % Plot echogram raw
  figure(1); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name','echogram (V)');
  imagesc([], pc_param(plot_wf).time*1e6, data{plot_wf}(:,:,plot_adc));
  title(sprintf('Time-space raw wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  ylabel('Time (us)');
  grid on;
  h = colorbar;
  set(get(h,'YLabel'),'String','Volts (V)');
  
  % =====================================================================
  % Plot echogram pulse compressed
  figure(2); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name','echogram (dBm)');
  imagesc([],time{plot_wf}*1e6, lp(abs(data_pc{plot_wf}(:,:,plot_adc)).^2/50) + 30);
  title(sprintf('Time-space PC wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  ylabel('Time (us)');
  grid on;
  h = colorbar;
  set(get(h,'YLabel'),'String','Power at ADC (dBm)');
  
  % =====================================================================
  % Radar equation for specular target
  num_tx = 7; power_per_tx = 25; ground_plane_gain = 4; tx_gain = 1;
  lambda_fc = c/195e6; refl_dB = -5; system_dB = -4.5; range = 21160*12*2.54/100;
  fprintf('Expected power from surface (all tx): %.1f dBm\n', ...
    lp(num_tx*power_per_tx*num_tx*tx_gain*ground_plane_gain ...
    *ground_plane_gain*lambda_fc^2/(4*pi)*10^(refl_dB/10)*10^(system_dB/10) / (4*pi*(2*range)^2)) + 30);
  num_tx = 1;
  fprintf('Expected power from surface (tx-cal): %.1f dBm\n', ...
    lp(num_tx*power_per_tx*num_tx*tx_gain*ground_plane_gain ...
    *ground_plane_gain*lambda_fc^2/(4*pi)*10^(refl_dB/10)*10^(system_dB/10) / (4*pi*(2*range)^2)) + 30);
    
  % =====================================================================
  % Plot pulse compressed A-scope
  figure(3); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name','a-scope pc');
  plot(pc_param(plot_wf).time*1e6, lp(abs(data{plot_wf}(:,rline,plot_adc)/rx_gain).^2/50) + 30);
  hold on;
  plot(time{plot_wf}*1e6, lp(abs(data_pc{plot_wf}(:,rline,plot_adc)/rx_gain).^2/50) + 30,'r');
  hold off;
  title(sprintf('A-Scope PC wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  grid on;
  xlabel('Time (us)');
  ylabel('Power at receiver input (dBm)');
  
  % =====================================================================
  % Plot raw A-scope (compare wf 1 and 2)
  figure(4); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name','a-scope raw');
  color_list = {'k' 'r' 'y' 'g' 'c' 'b' 'm'};
  for wf = 1:min(2,length(data))
    plot(pc_param(wf).time*1e6, data{wf}(:,rline,plot_adc),color_list{wf});
    hold on;
  end
  hold off;
  title(sprintf('A-Scope Raw wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  grid on;
  xlabel('Time (us)');
  ylabel('Volts at receiver input');
  
  % =====================================================================
  % Digital down conversion and envelope plot
  fc = (pc_param(wf).f1+pc_param(wf).f0)/2;
  BW = pc_param(wf).f1-pc_param(wf).f0;
  data_BB = double(data{wf}(:,rline,plot_adc) .* exp(1i*2*pi*fc*pc_param(wf).time));
  [Blpf,Alpf] = butter(2,BW/fs);
  figure(5); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name','envelope');
  plot(abs(filtfilt(Blpf,Alpf,data_BB)));
  title(sprintf('Envelope wf%d adc%d %s', plot_wf, plot_adc, fn_name),'Interpreter','none');
  %xlim([1800 2850]);

end





