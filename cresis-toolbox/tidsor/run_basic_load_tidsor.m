% Script run_basic_load_tidsor
%
% Several examples of using basic_load_tidsor.
%
% Author: John Paden
%
% See also basic_load_tidsor

physical_constants;
close all;
run_basic_load_tidsor_tstart = tic;

% =======================================================================
% User Settings: Select Example
% =======================================================================
run_example = 2;

% =======================================================================
% User Settings: To Load Data
% =======================================================================
fs = 1e9/8;
base_path = '/cresis/data2/TIDSoR/2011_Chile/Fligth_test_Sat_Jul23_2011/';
fn_name = 'data.07232011.0001.dat';
coh_ave = 100;

% =======================================================================
% Load data
% =======================================================================
fprintf('============================================================\n');
fprintf('Loading data (%.1f sec)\n', toc(run_basic_load_tidsor_tstart));

fn = fullfile(base_path,fn_name);
[hdr,data] = basic_load_tidsor(fn, struct('clk',fs,'first_byte',2^26,'coh_ave',coh_ave));

[fn_dir fn_name] = fileparts(fn);

fprintf('Running example %d (%.1f sec)\n', run_example, toc(run_basic_load_tidsor_tstart));

if run_example == 1
  % =======================================================================
  % =======================================================================
  % Example 1
  %  - Plot data in time and frequency
  % =======================================================================
  % =======================================================================
  
  % =======================================================================
  % User Settings: Radar parameters and What to Plot
  % =======================================================================
  adc_bits = 14;
  Vpp_scale = 2;
  adc_SNR_dB = 70;
  rx_gain = 10^((35)/20);
  % Choose which range line to plot and other settings
  plot_rline = 2;
  
  % =======================================================================
  % Convert from quantization to voltage @ ADC
  for wf = 1:size(data,3)
    data(:,:,wf) = data(:,:,wf) ...
      * Vpp_scale/2^adc_bits ...
      * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums;
  end
  
  % Remove DC-bias and ADC unsigned sample offset (e.g. unsigned 16 bit
  % numbers will have 2^15 offset). Use the first record to determine
  % the correct offset.
  for wf = 1:size(data,3)
    data(:,:,wf) = data(:,:,wf) - median(data(:,1,wf));
  end
  
  for wf = 1:size(data,3)
    
    % Create time and freq axes
    time = hdr.wfs(wf).t0 + (0:size(data,1)-1).'/fs;
    dt = time(2) - time(1);
    Nt = length(time);
    df = 1/(Nt*dt);
    freq = (0:df:(Nt-1)*df).';
    
    % Plot data
    figure(wf); clf;
    plot(time*1e6, lp( abs(data(:,plot_rline,wf)).^2/50 ) + 30);
    title(sprintf('A-scope raw rline%d wf%d %s', plot_rline, wf, fn_name),'Interpreter','none');
    xlabel('Fast time (us)');
    ylabel('Power at rx input (dBm)');
    grid on;
    xlim([time(1) time(end)]*1e6);
    %saveas(wf,sprintf('/cresis/scratch1/paden/TIDSoR/file_1_wf_%d_raw.fig',wf));
    
    figure(10+wf); clf;
    plot(freq/1e6,lp( fft(data(:,plot_rline,wf)) ));
    title(sprintf('A-scope raw-fft rline%d wf%d %s', plot_rline, wf, fn_name),'Interpreter','none');
    xlabel('Frequency (MHz)');
    ylabel('Relative power (dB)');
    grid on;
    xlim([freq(1) fs/2]/1e6);
    %saveas(10+wf,sprintf('/cresis/scratch1/paden/TIDSoR/file_1_wf_%d_rawfft.fig',wf));
    
  end
  
elseif run_example == 2
  % =======================================================================
  % =======================================================================
  % Example 2
  %  - Basic processing and plot of echogram
  % =======================================================================
  % =======================================================================
  
  % =======================================================================
  % User Settings: Radar parameters and What to Plot
  % =======================================================================
  adc_bits = 14;
  Vpp_scale = 2;
  adc_SNR_dB = 70;
  rx_gain = 10^((45)/20);
  incoh_ave = 10;
  plot_rline = 2;
  
  % =======================================================================
  % Convert from quantization to voltage @ ADC
  for wf = 1:size(data,3)
    data(:,:,wf) = data(:,:,wf) ...
      * Vpp_scale/2^adc_bits ...
      * 2^hdr.wfs(wf).bit_shifts/hdr.wfs(wf).presums;
  end
  
  % Remove DC-bias and ADC unsigned sample offset (e.g. unsigned 16 bit
  % numbers will have 2^15 offset). Use the first record to determine
  % the correct offset.
  for wf = 1:size(data,3)
    data(:,:,wf) = data(:,:,wf) - median(data(:,1,wf));
  end
  
  fc(1) = 8.4e6;
  fc(2) = 12.5e6;
  BW = 1e6;
  
  fir_order = round(2*fs/BW);
  fir_order = 128;
  B = fir1(fir_order, BW/(fs/2));
  A = 1;
  
  clear out_data;
  for wf = 1:size(data,3)
    % Create time and freq axes
    time = hdr.wfs(wf).t0 + (0:size(data,1)-1).'/fs;
    dt = time(2) - time(1);
    Nt = length(time);
    df = 1/(Nt*dt);
    freq = (0:df:(Nt-1)*df).';
    for rline = 1:size(data,2)
      % Complex baseband the data
      data(:,rline,wf) = data(:,rline,wf) .* exp(-j*2*pi*fc(wf)*time);
      % LPF the baseband data
      data(:,rline,wf) = filtfilt(B,A,data(:,rline,wf));
    end
    out_data(:,:,wf) = fir_dec(abs(data(:,:,wf)).^2,incoh_ave);
    
    % Plot echogram
    figure(wf); clf;
    imagesc([],time*1e6,lp(out_data(:,:,wf),1));
    xlabel('Range line');
    ylabel('Fast time (us)');
    title(sprintf('Echogram wf%d %s', wf, fn_name),'Interpreter','none');
    grid on;
    %saveas(wf,sprintf('/cresis/scratch1/paden/TIDSoR/file_1_wf_%d_echo.fig',wf));
    
    % Plot A-scope
    figure(10+wf); clf;
    plot(time*1e6,lp(out_data(:,plot_rline,wf),1));
    xlabel('Fast time (us)');
    ylabel('Relative power (dB)');
    title(sprintf('A-scope rline%d wf%d %s', plot_rline, wf, fn_name),'Interpreter','none');
    grid on;
    %saveas(10+wf,sprintf('/cresis/scratch1/paden/TIDSoR/file_1_wf_%d_ascope.fig',wf));
  end
    
end

fprintf('Example done (%.1f sec)\n', toc(run_basic_load_tidsor_tstart));

return;
