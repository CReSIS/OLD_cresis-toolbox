% script basic_noise_analysis
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

%% Determine which file to load
good_mask = logical(zeros(size(param.base_dir_search)));
for base_dir_idx = 1:length(param.base_dir_search)
  base_dir = param.base_dir_search{base_dir_idx};
  if exist(base_dir,'dir')
    good_mask(base_dir_idx) = true;
  end
end
param.base_dir_search = param.base_dir_search(good_mask);

global g_basic_noise_analysis_base_dir;
base_dir = [];
if length(param.base_dir_search) >= 1
  default_base_dir_idx = [];
  if isempty(g_basic_noise_analysis_base_dir)
    g_basic_noise_analysis_base_dir = param.base_dir_search{1};
  end
  for base_dir_idx = 1:length(param.base_dir_search)
    fprintf('(%d): %s', base_dir_idx, param.base_dir_search{base_dir_idx});
    if strcmp(param.base_dir_search{base_dir_idx},g_basic_noise_analysis_base_dir)
      fprintf(' *');
      default_base_dir_idx = base_dir_idx;
    end
    fprintf('\n');
  end
  fprintf('(%d): Custom', base_dir_idx+1);
  if isempty(default_base_dir_idx)
      fprintf(' *');
  end
  fprintf('\n');
  base_dir_idx = input('More than one base directory exists, choose one: ');
  if isempty(base_dir_idx) && ~isempty(default_base_dir_idx)
    base_dir_idx = default_base_dir_idx;
  end
  if base_dir_idx <= length(param.base_dir_search)
    base_dir = param.base_dir_search{base_dir_idx};
  end
end

if isempty(base_dir)
  while ~exist(base_dir,'dir')
    base_dir = input('Enter custom directory path: ','s');
    if ~exist(base_dir,'dir')
      warning('Does not exist: %s', base_dir);
    end
  end
end
g_basic_noise_analysis_base_dir = base_dir;

global g_file_select;
global g_basic_noise_analysis_fn;
if strcmpi(param.file_search_mode,'last_file')
  fns = get_filenames(base_dir,'','','.bin',struct('recursive',true));
  if isempty(fns)
    error('No data files: %s\n', base_dir);
  end
  fns_idxs = max(1,length(fns)-9) : length(fns);
  if isempty(g_file_select)
    g_file_select = max(1,length(fns_idxs)-1);
  end
  for fn_idx = 1:length(fns_idxs)
    fprintf('(%d): %s', fn_idx, fns{fns_idxs(fn_idx)});
    if g_file_select == fn_idx
      fprintf(' *');
    end
    fprintf('\n');
  end
  done = false;
  while ~done
    try
      user_fn_idx = input('Choose one: ');
      if isempty(user_fn_idx)
        user_fn_idx = g_file_select;
      end
      fn = fns{fns_idxs(user_fn_idx)};
      done = true;
    catch
      user_fn_idx = g_file_select;
      fn = fns{fns_idxs(user_fn_idx)};
      done = true;
    end
  end
else
  fn = '';
  while ~exist(fn,'file')
    fn = input(sprintf('Filename [%s]: ', g_basic_noise_analysis_fn),'s');
    if isempty(fn)
      fn = g_basic_noise_analysis_fn;
    end
    fn = get_filename(base_dir,'',fn,'.bin',struct('recursive',true));
  end
end
g_basic_noise_analysis_fn = fn;

%% Load the chosen file(s)
tstart = tic;
fprintf('Loading data (%.1f sec)\n', toc(tstart));
% Load the data (disable if you have already loaded)
clear data;
clear num_rec;
if strcmpi(param.radar_name,'mcords')
  for adc_idx = 1:length(param.adcs)
    adc = param.adcs(adc_idx);

    % May need to adjust base_path for non-standard directory structures
    base_path = fullfile(param.base_path, sprintf('chan%d',adc), ...
      param.seg);
    file_midfix = sprintf('r%d-%d.',param.radar_num,adc);
    file_suffix = sprintf('.%04d.dat',param.data_file_num);
    fprintf('  Path: %s\n', base_path);
    fprintf('  Match: mcords*%s*%s\n', file_midfix, file_suffix);
    fn = get_filename(base_path,'mcords',file_midfix,file_suffix);
    if isempty(fn)
      fprintf('  Could not find any files which match\n');
      return;
    end
    fprintf('  Loading file %s\n', fn);
    [hdr,data_tmp] = basic_load_mcords(fn, struct('clk',1e9/9,'first_byte',2^26));
    data(:,:,adc_idx) = data_tmp{param.wf}(1:end-1,1:min(size(data_tmp{param.wf},2),param.rlines(2)));
  end
  data = data - median(data(:,1));
%   basic_remove_mcords_digital_errors;
elseif any(strcmpi(param.radar_name,{'mcords2','mcords3'}))
  % test1_1.dat0
  %   testA_N.datB
  %   A = acquisition number
  %   N = file number
  %   B = board number
  % Map ADCs to board numbers
  for board = 0:3
    if any(board == floor((param.img(:,2)-1)/4))
      get_adcs = board*4 + (1:4);
      file_prefix = sprintf('%s_%d_',param.radar_name,board);
      if isempty(acquisition_num)
        file_suffix = sprintf('%04d.bin',file_num);
      else
        file_suffix = sprintf('%s_%04d.bin',acquisition_num,file_num);
      end
      base_path = fullfile(param.base_path, sprintf('board%d',board), ...
        seg);
      fprintf('  Path: %s\n', base_path);
      fprintf('  Match: %s*%s\n', file_prefix, file_suffix);
      fn = get_filenames(base_path, file_prefix, '', file_suffix);
      if isempty(fn)
        fprintf('  Could not find any files which match\n');
        return;
      end
      fn = fn{end};
      fprintf('  Loading file %s\n', fn);
      % Fix get_filenames     'The filename, directory name, or volume label syntax is incorrect.'
      if strcmpi(param.radar_name,'mcords2')
        [hdr,data_tmp] = basic_load_mcords2(fn,struct('clk',fs));
      else
        [hdr,data_tmp] = basic_load_mcords3(fn,struct('clk',fs));
      end
      for get_adc_idx = 1:length(get_adcs)
        adc = get_adcs(get_adc_idx);
        for wf_adc_idx = find(param.img(:,2) == adc)
          wf = param.img(wf_adc_idx,1);
          fprintf('  Loading wf %d, adc %d\n', wf, adc);
          if ~exist('num_rec','var')
            % Since each file may have slightly different numbers of
            % records we do this
            num_rec = size(data_tmp{wf},2) - 1;
          end
          data(:,:,wf_adc_idx) = data_tmp{wf}(:,1:num_rec,get_adc_idx);
        end
      end
    end
  end
  
elseif any(strcmpi(param.radar_name,{'mcords4','mcords5'}))
  file_idx = 1;
  epri_intersect = [];
  % adcs: a list of the adcs that we are loading
  adcs = unique(param.img(:,2));
  
  for adc = reshape(adcs,[1 length(adcs)])
    
    [fn_dir,fn_name] = fileparts(g_basic_noise_analysis_fn);
    fn_dir = fileparts(fn_dir);
    fn_name(9:10) = sprintf('%02d',adc);
    fn = fullfile(fn_dir,sprintf('chan%d',adc),[fn_name,'.bin']);
    
    file_name_list{file_idx} = fn;
    fprintf('  Loading file %s\n', fn);
    % Load the data file
    if strcmp(param.radar_name,'mcords4')
      [hdr,data_tmp] = basic_load_mcords4(fn,struct('clk',default.radar.fs/4,'recs',param.recs));
    else
      [hdr,data_tmp] = basic_load_mcords5(fn,struct('clk',default.radar.fs,'recs',param.recs));
    end
    % Remove extra records to help reduce total memory usage
%     if isfield(param,'rlines') && ~isempty(param.rlines)
%       for wf = 1:length(data_tmp)
%         data_tmp{wf} = data_tmp{wf}(:,param.rlines,:);
%       end
%       hdr.utc_time_sod = hdr.utc_time_sod(param.rlines);
%       hdr.epri = hdr.epri(param.rlines);
%     end
    % Map each of the read waveforms needed into the correct output
    for wf_adc_idx = 1:size(param.img,1)
      % wf,adc: pair of values for this entry in param.img
      wf = param.img(wf_adc_idx,1);
      if adc == abs(param.img(wf_adc_idx,2));
        % This pair needs to be loaded, insert into output array... handle
        % mismatched EPRIs using intersect function.
        if isempty(epri_intersect)
          epri_intersect = hdr.epri;
          if imag(wf) == 0
            data(:,:,wf_adc_idx) = data_tmp{wf}(:,:);
          else
            data(:,:,wf_adc_idx) = data_tmp{abs(wf)}(:,:) ...
              + sign(wf) * data_tmp{abs(wf)+1}(:,:);
          end
          hdr_utc_time_sod = hdr.utc_time_sod;
        else
          [epri_intersect data_idx data_tmp_idx] = intersect(epri_intersect,hdr.epri);
          data = data(:,data_idx,:);
          hdr_utc_time_sod = hdr.utc_time_sod(data_tmp_idx);
          if imag(wf) == 0
            data(:,:,wf_adc_idx) = data_tmp{wf}(:,data_tmp_idx);
          else
            data(:,:,wf_adc_idx) = data_tmp{abs(wf)}(:,data_tmp_idx) ...
              + sign(wf) * data_tmp{abs(wf)+1}(:,data_tmp_idx);
          end
        end
      end
    end
  end
  hdr.utc_time_sod = hdr_utc_time_sod;
  
  xml_version = 2.0;
  cresis_xml_mapping;
  
    %% Read XML files in this directory
  settings = read_ni_xml_directory(fn_dir,'',false);
  finfo = fname_info_mcords2(fn);
  
  settings_idx = find(cell2mat({settings.datenum}) < finfo.datenum,1,'last');
  settings = settings(settings_idx);

  BW = settings.DDS_Setup.Waveforms(wf).Stop_Freq(1) - settings.DDS_Setup.Waveforms(wf).Start_Freq(1);
  atten = double(settings.DDS_Setup.Waveforms(wf).Attenuator_1(1) + settings.DDS_Setup.Waveforms(wf).Attenuator_2(1));
  rx_gain = default.radar.rx_gain .* 10.^(-atten/20);
  
end


[fn_dir fn_name] = fileparts(fn);
if ~isfield(param,'seg') || isempty(param.seg)
  param.seg = -1;
end
clear data_tmp;

if ~isfield(param,'noise_rbins') || isempty(param.noise_rbins)
  noise_rbins = 1:size(data,1);
elseif param.noise_rbins(end) > size(data,1)
  noise_rbins = param.noise_rbins(1):size(data,1);
elseif param.noise_rbins(1) <= 0
  noise_rbins = size(data,1)+param.noise_rbins;
else
  noise_rbins = param.noise_rbins(1):param.noise_rbins(end);
end

%% Convert from quantization to voltage @ ADC
data = data ...
  * default.radar.adc_full_scale/2^default.radar.adc_bits ...
  * 2^hdr.wfs(abs(param.img(1,1))).bit_shifts / hdr.wfs(abs(param.img(1,1))).presums;

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
  fprintf('Expected Rx noise floor @ ADC %.1f dBm\n', lp(BoltzmannConst*290*BW*default.radar.noise_figure*rx_gain^2) + 3 +30);
else
  fprintf('Expected ADC noise floor @ ADC %.1f dBm\n', lp((default.radar.adc_full_scale/2/sqrt(2))^2/50)+30 - default.radar.adc_SNR_dB );
  fprintf('Expected Rx noise floor @ ADC %.1f dBm\n', lp(BoltzmannConst*290*BW*default.radar.noise_figure*rx_gain^2)+30);
end
fprintf('Expected levels only valid for param.presums = 1\n');
fprintf('Noise power is in dBm:\n')
for adc_idx = 1:size(data,3)
  if adc_idx > 1
    fprintf('\t');
  end
  fprintf('%.1f', ...
    lp(mean(mean(abs(data(noise_rbins,:,adc_idx)).^2/50)) * hdr.wfs(abs(param.img(1,1))).presums, 1) +30 );
end
fprintf('\n');

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
    
    clear pc_param;
    if strcmp(param.radar_name,'mcords5') && isfield(hdr,'DDC') && hdr.DDC(1) >= 2
      pc_param.time = hdr.wfs(abs(param.img(1,1))).t0 + (0:size(fir_data,1)-1)/default.radar.fs*2^hdr.DDC(1);
      dt = pc_param.time(2) - pc_param.time(1);
      Nt = length(pc_param.time);
      df = 1/(Nt*dt);
      freq = double(settings.DDC_Ctrl.NCO_freq)*1e6 + (-floor(Nt/2)*df : df : floor((Nt-1)/2)*df);
      
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
      pc_param.time = hdr.wfs(abs(param.img(1,1))).t0 + (0:size(fir_data,1)-1)/fs;
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

