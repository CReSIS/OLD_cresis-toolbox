% Script run_load_data
%
% Several examples of using load_data
%
% Author: John Paden
%
% See also load_data, pulse_compress

run_example = 1;
clear data;

if run_example == 1
  % =======================================================================
  % Setup loading parameters for example 1
  %  - Plot data in time-space, time-wavenumber, freq-space, and
  %    freq-wavenumber domains
  % =======================================================================
  
  param = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140401_03');
  
  % Determine which records you want to load:
  frames_fn = '';
  frames_fn = ct_filename_support(param,frames_fn,'frames');
  load(frames_fn);
  frm = 37;
  param.load_data.recs = frames.frame_idxs(frm) + [0 99];
  
  if ~isfield(param.records,'records_fn')
    param.records.records_fn = '';
  end
  param.load_data.records_fn = param.records.records_fn;
  %   param.load_data.imgs = {[-1j 5]};
  %   param.load_data.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 2 9; 2 10; 2 11; 2 12; 2 13; 2 14; 2 15; 2 16]};
  param.load_data.imgs = {[1 2]};
  param.load_data.pulse_comp         = true;
  param.load_data.ft_dec             = true; % Fast-time decimation
  param.load_data.ft_wind            = @hanning;
  param.load_data.ft_wind_time       = false; % Apply window on time domain chirp
  param.load_data.presums            = 1; % Coherent averaging
  param.load_data.combine_rx         = false;
  param.load_data.pulse_rfi.en       = false;
  param.load_data.pulse_rfi.inc_ave  = 101;
  param.load_data.pulse_rfi.thresh_scale = 10^(13/10);
  param.load_data.trim_vals          = [0 0];
  
  % Load data
  [data,hdr] = load_data(param);
  
  % Plot data
  img = 1;
  wf_adc_idx = 1;
  wf = abs(param.load_data.imgs{img}(wf_adc_idx,1));
  
  figure(1); clf;
  imagesc([],hdr.wfs(wf).time, ...
    lp(data{img}(:,:,wf_adc_idx)));
  title('Time-space domain');
  grid on;
  colorbar
  
  noise_bins = 190:230;
  
  figure(2); clf;
  imagesc([],fftshift(hdr.wfs(wf).freq), ...
    lp(fftshift(fft(data{img}(noise_bins,:,wf_adc_idx)),1)));
  title('Frequency-space domain (noise)');
  grid on;
  colorbar
  
  figure(3); clf;
  imagesc([],hdr.wfs(wf).time, ...
    lp(fft(data{img}(:,:,wf_adc_idx),[],2)));
  title('Time-wavenumber domain');
  grid on;
  colorbar
  
  figure(4); clf;
  imagesc([],fftshift(hdr.wfs(wf).freq), ...
    lp(fftshift(fft2(data{img}(noise_bins,:,wf_adc_idx)),1)));
  title('Frequency-wavenumber domain (noise)');
  grid on;
  colorbar
  
elseif run_example == 2
  %% Setup loading parameters for example 2
  %  - Plot raw data in time-space, pulse compress externally with
  %    pulse_compress and plot this
  % =======================================================================
  
  param = read_param_xls(ct_filename_param('rds_param_2016_Greenland_Polar6.xls'),'20160413_04');
  
  % Determine which records you want to load:
  frames_fn = '';
  frames_fn = ct_filename_support(param,frames_fn,'frames');
  load(frames_fn);
  frm = 1;
  param.load_data.recs = frames.frame_idxs(frm) - 1 + [10000 10250];
  
  param.load_data.records_fn = param.records.records_fn;
  param.load_data.imgs = {[ones(24,1) (1:24).'],[2*ones(24,1) (1:24).'],[3*ones(24,1) (1:24).']};
  param.load_data.pulse_comp         = false;
  param.load_data.ft_dec             = false; % Fast-time decimation
  param.load_data.ft_wind            = @hanning;
  param.load_data.ft_wind_time       = false; % Apply window on time domain chirp
  param.load_data.presums            = 1; % Coherent averaging
  param.load_data.combine_rx         = false;
  param.load_data.pulse_rfi.en       = false;
  param.load_data.pulse_rfi.inc_ave  = 101;
  param.load_data.pulse_rfi.thresh_scale = 10^(13/10);
  param.load_data.trim_vals          = [1 1];
  param.load_data.raw_data           = true;
  
  %% Load data
  [data,hdr] = load_data(param);

  %% Print out DC values and create DC adjust files
  if 0
    adc_mean = [];
    for img = 1:length(param.load_data.imgs)
      wf = param.load_data.imgs{img}(1,1);
      adc_mean = [];
      for wf_adc = 1:length(param.load_data.imgs{img})
        adc = param.load_data.imgs{img}(wf_adc,2);
        dd=data{img}(end-499:end,:,adc); % Assumes bottom 500 range bins contain only noise (this should always be verified by looking at the echogram)
        adc_mean{img}(wf_adc) = mean(dd(:));
      end
      fprintf('%s wf %d\n', param.day_seg, wf);
      fprintf('%.3f\t', real(adc_mean{img}(1:end-1)));
      fprintf('\n');
      fprintf('%.3f\t', imag(adc_mean{img}(1:end-1)));
      fprintf('\n');
      fprintf('[');
      fprintf('%.3f%+.3fj ', real(adc_mean{img}(1:end-1)), imag(adc_mean{img}(1:end-1)));
      fprintf('%.3f%+.3fj', real(adc_mean{img}(end)), imag(adc_mean{img}(end)));
      fprintf(']\n');
      DC_adjust = adc_mean{img};
      DC_adjust_file_en = false;
      if DC_adjust_file_en
        %% Create the DC adjust files
        fn = fullfile(ct_filename_out(param,'noise','',1),sprintf('DC_%s_wf%d.mat',param.day_seg,wf));
        save(fn, '-v6', 'DC_adjust');
      end
    end
  end

  %% Plot data
  img = 1;
  wf_adc_idx = 1;
  wf = abs(param.load_data.imgs{img}(wf_adc_idx,1));
  
  figure(1); clf;
  imagesc([],hdr.wfs(wf).time, ...
    lp(data{img}(:,:,wf_adc_idx)));
  title('Time-space domain');
  grid on;
  colorbar
  
  % Pulse compress data
  clear pc_param;
  pc_param.f0 = hdr.wfs(wf).f0;
  pc_param.f1 = hdr.wfs(wf).f1;
  pc_param.Tpd = hdr.wfs(wf).Tpd;
  pc_param.time = hdr.wfs(wf).time;
  pc_param.tukey = param.radar.wfs(wf).tukey;
  pc_param.presums = hdr.wfs(wf).presums;
  pc_param.bit_shifts = max(0,(14+log(hdr.wfs(wf).presums)/log(2))-16);
  [data_pc,time] = pulse_compress(data{img}(:,:,wf_adc_idx),pc_param);
  
  figure(2); clf;
  imagesc([],time, lp(data_pc));
  title('Time-space domain');
  grid on;
  colorbar
  
elseif run_example == 3
  % =======================================================================
  % Setup loading parameters for example 3
  %  - Examines waveform imbalance
  % =======================================================================
  
  param_fn = '/cresis/projects/dev/csarp_support/documents/mcords_param_2010_Greenland_DC8.xls';
  param = read_param_xls(param_fn,'20100324_01');
  
  % Determine which records you want to load:
  frames_fn = '';
  frames_fn = ct_filename_support(param,frames_fn,'frames');
  load(frames_fn);
  frm = 10;
  param.load_data.recs = frames.frame_idxs(frm) - 1 + [4001 8000];
  
  param.load_data.records_fn = param.records.records_fn;
  param.load_data.imgs = {[1 1; 1 2; 1 3; 1 4; 1 5], [2 1; 2 2; 2 3; 2 4; 2 5]};
  param.load_data.imgs = {[1 1], [2 1]};
  param.load_data.pulse_comp         = true;
  param.load_data.ft_dec             = true; % Fast-time decimation
  param.load_data.ft_wind            = @hanning;
  param.load_data.ft_wind_time       = false; % Apply window on time domain chirp
  param.load_data.presums            = 50; % Coherent averaging
  param.load_data.combine_rx         = true;
  param.load_data.pulse_rfi.en       = true;
  param.load_data.pulse_rfi.inc_ave  = 101;
  param.load_data.pulse_rfi.thresh_scale = 10^(13/10);
  param.load_data.trim_vals          = [1 1];
  
  % Load data
  [data,hdr] = load_data(param);
  
  % Plot data
  sig_time = [17e-6 20e-6];
  noise_time = [60e-6 70e-6];
  for img = 1:length(data)
    for wf_adc_idx = 1:size(data{img},3)
      wf = param.load_data.imgs{img}(wf_adc_idx,1);
      
      figure(10*(img-1) + wf_adc_idx); clf;
      imagesc([],hdr.wfs(wf).time, ...
        lp(data{img}(:,:,wf_adc_idx)));
      title('Time-along track domain');
      grid on;
      colorbar
      
      figure(100 + 10*(img-1) + wf_adc_idx); clf;
      sig_bins = find(hdr.wfs(wf).time > sig_time(1) & hdr.wfs(wf).time < sig_time(2));
      noise_bins = find(hdr.wfs(wf).time > noise_time(1) & hdr.wfs(wf).time < noise_time(2));
      sig_val = max(abs(data{img}(sig_bins,:,wf_adc_idx)).^2);
      noise_val = mean(abs(data{img}(noise_bins,:,wf_adc_idx)).^2);
      plot(lp(sig_val));
      hold on;
      plot(lp(noise_val),'r');
      hold off;
      title('SNR-along track domain');
      grid on;
      colorbar
    end
  end
  
elseif run_example == 4
  % =======================================================================
  % Setup loading parameters for example 4
  %  - Receiver channel equalization
  % =======================================================================
  
  %param_fn = 'C:\csarp_support\documents\mcords_param_2011_Greenland_P3.xls';
  %param_fn = '/mnt/scratch1/csarp_support/documents/mcords_param_2011_Greenland_P3.xls';
  %param_fn = '/mnt/scratch1/csarp_support/documents/seaice_param_2011_Greenland_P3.xls';
  param_fn = 'C:\csarp_support\documents\mcords_param_2011_Greenland_P3.xls';
  param = read_param_xls(param_fn,'20110317_01');
  
  % Determine which records you want to load:
  frames_fn = '';
  frames_fn = ct_filename_support(param,frames_fn,'frames');
  load(frames_fn);
  frm = 2;
  param.load_data.recs = frames.frame_idxs(frm) - 1 + [5000 17000];
  
  param.load_data.records_fn = param.records.records_fn;
  %   param.load_data.imgs = {[-1j 5]};
  param.load_data.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 2 9; 2 10; 2 11; 2 12; 2 13; 2 14; 2 15; 2 16]};
  %   param.load_data.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8]};
  %   param.load_data.imgs = {[1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8]};
  %   param.load_data.imgs = {[2 9; 2 10; 2 11; 2 12]};
  %   param.load_data.imgs = {[2 13; 2 14; 2 15; 2 16]};
  param.load_data.pulse_comp         = true;
  param.load_data.ft_dec             = true; % Fast-time decimation
  param.load_data.ft_wind            = @hanning;
  param.load_data.ft_wind_time       = false; % Apply window on time domain chirp
  param.load_data.presums            = 10; % Coherent averaging
  param.load_data.combine_rx         = false;
  param.load_data.pulse_rfi.en       = false;
  param.load_data.pulse_rfi.inc_ave  = 101;
  param.load_data.pulse_rfi.thresh_scale = 10^(13/10);
  param.load_data.trim_vals          = [0 0];
  
  % Load data
  [data,hdr] = load_data(param);
  data = data{1};
  
  old_param = param;
  clear param;
  param.plot_en = true;
  param.adcs = old_param.load_data.imgs{1}(:,2);
  fn_name = sprintf('%s recs %d to %d', old_param.day_seg, ...
    old_param.load_data.recs);
  
  % =====================================================================
  % Setup parameters for
  
  % .caxis = Color axis limits (leave empty first time since this causes it to use the
  % defaults).
  param.caxis = [];
  %param.caxis = [50 120];
  
  % .ylim = Leave empty the first time (it just uses the default limits then)
  param.ylim = [];
  %param.ylim = [1 500];
  
  % .xlim = Leave empty the first time (it just uses the default limits then)
  %   These limits are in range lines post presumming
  param.xlim = [];
  
  % .ref_adc = Reference receive channel (surface location determined from this
  %   channel and all phase measurements made relative to it)
  param.ref_adc = 3;
  % .rlines = Range lines to process
  %   These are range lines post presumming
  param.rlines = [1 inf];
  % .rbins = Range bins to search for surface in (NEEDS TO BE SET IF NOT BLANKED)
  param.rbins = [500 inf];
  % .noise_rbins = Range bins to use for noise power calculation (THIS OFTEN NEEDS TO BE SET)
  param.noise_rbins = 700:1000;
  
  % .snr_threshold = SNR threshold in dB (range lines exceeding this
  %   SNR are included in the estimate)
  param.snr_threshold = 25;
  
  if param.rbins(2) > size(data,1)
    param.rbins(2) = size(data,1);
  end
  if param.rlines(2) > size(data,2)
    param.rlines(2) = size(data,2);
  end
  param.rbins = param.rbins(1):param.rbins(2);
  param.rlines = param.rlines(1):param.rlines(2);
  
  % =======================================================================
  % Echogram plots
  % =======================================================================
  if param.plot_en
    for adc_idx = 1:size(data,3)
      figure(adc_idx); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('E %d',adc_idx));
      imagesc(lp(data(:,:,adc_idx)));
      title(sprintf('ADC %d File %s\nTime-Space Relative Power', param.adcs(adc_idx), fn_name),'Interpreter','none');
      grid on;
      colorbar
      if ~isempty(param.caxis)
        caxis(param.caxis);
      end
      if ~isempty(param.ylim)
        ylim(param.ylim);
      end
      if ~isempty(param.xlim)
        xlim(param.xlim);
      end
    end
  end
  
  % =======================================================================
  % Surface tracker
  % =======================================================================
  surf_data = filter2(ones(1,10),abs(data(:,:,param.ref_adc).^2));
  [surf_vals surf_bins] = max(surf_data(param.rbins,param.rlines));
  surf_bins = param.rbins(1)-1 + surf_bins;
  
  if param.plot_en
    for adc_idx = 1:size(data,3)
      figure(adc_idx);
      hold on;
      plot(param.rlines, surf_bins,'k');
      hold off;
    end
  end
  
  % =======================================================================
  % Noise power estimate and SNR threshold
  % =======================================================================
  noise_power = mean(mean(abs(data(param.noise_rbins,param.rlines)).^2));
  wf = 1;
  clear tx_phases tx_powers;
  for adc_idx = 1:size(data,3)
    for rline_idx = 1:length(param.rlines)
      rline = param.rlines(rline_idx);
      tx_phases(adc_idx,rline_idx) = data(surf_bins(rline_idx),rline,adc_idx);
      tx_powers(adc_idx,rline_idx) = abs(data(surf_bins(rline_idx),rline,adc_idx)).^2;
    end
  end
  tx_snr = tx_powers ./ noise_power;
  good_meas = lp(tx_snr) > param.snr_threshold;
  good_rlines = zeros(size(param.rlines));
  good_rlines(sum(good_meas) == size(tx_snr,1)) = 1;
  good_rlines = logical(good_rlines);
  
  num_good_rlines = sum(good_rlines);
  fprintf('Number of good range lines: %d out of %d\n', num_good_rlines, length(good_rlines));
  fprintf('========================================================\n');
  
  % =======================================================================
  % Amplitude Settings
  % =======================================================================
  clear delta_power;
  fprintf('Relative power for each waveform (dB)\n');
  for adc_idx = 1:size(data,3)
    ref_power = tx_powers(adc_idx,:)./tx_powers(param.ref_adc,:);
    if param.plot_en
      figure(10+adc_idx); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Pow %d',adc_idx));
      plot(lp(ref_power));
    end
    delta_power(adc_idx) = lp(median(ref_power(good_rlines)));
    fprintf('%10.2f\n', delta_power(adc_idx));
    %   fprintf('WF %d: relative power: %10.2f dB\n', wf, lp(median(ref_power(good_rlines))));
    %   fprintf('    Std. dev. power: %.2f dB\n', lp(std(ref_power(good_rlines))));
    ref_power(~good_rlines) = NaN;
    if param.plot_en
      hold on;
      plot(lp(ref_power),'ro');
      hold off;
      title(sprintf('Relative Power (%d to ref %d)', adc_idx, param.ref_adc));
    end
  end
  fprintf('Recommended amplitude settings:\n');
  fprintf('%.1f\t', delta_power(1:end-1));
  fprintf('%.1f', delta_power(end));
  fprintf('\n');
  fprintf('========================================================\n');
  
  % =======================================================================
  % Phase Settings
  % =======================================================================
  clear ref_phase_median;
  for adc_idx = 1:size(data,3)
    ref_phase = angle(tx_phases(adc_idx,:)./tx_phases(param.ref_adc,:));
    if param.plot_en
      figure(20+adc_idx); clf; set(gcf,'WindowStyle','docked','NumberTitle','off','Name',sprintf('Ang %d',adc_idx));
      plot(ref_phase);
      ylim([-pi pi]);
    end
    ref_phase_median(adc_idx) = median(ref_phase(good_rlines));
    fprintf('ADC %d: relative phase: %10.4f rad, %10.1f deg\n', adc_idx, ...
      median(ref_phase(good_rlines)), median(ref_phase(good_rlines))*180/pi);
    fprintf('    Std. dev. phase: %.4f rad, %.1f deg\n', ...
      std(ref_phase(good_rlines)), std(ref_phase(good_rlines))*180/pi);
    ref_phase(~good_rlines) = NaN;
    if param.plot_en
      hold on;
      plot(ref_phase,'ro');
      hold off;
      title(sprintf('Relative Phase (%d to ref %d)', adc_idx, param.ref_adc));
    end
  end
  
  fprintf('Recommended phase settings:\n');
  fprintf('%.1f\t', ref_phase_median(1:end-1)*180/pi);
  fprintf('%.1f', ref_phase_median(end)*180/pi);
  fprintf('\n');
  fprintf('========================================================\n');
  
elseif run_example == 5
  % =======================================================================
  % Setup loading parameters for example 5
  %  - Layer data for Joe MagGregor (pulsed compressed, no presums)
  % =======================================================================
  
  base_dir = '/cresis/snfs1/scratch/paden/winnie_chu/';
  % /cresis/scratch2/mdce/OIB_internal_layer_example/';
  
  param_fn = ct_filename_param('rds_param_2011_Greenland_P3.xls');
  %   param = read_param_xls(param_fn,'20110502_01');
  params = read_param_xls(param_fn,'20110416_02'); % USE GENERIC FIELD IN params(1).cmd worksheet
  params(1).cmd.frms = 19;
  params(1).cmd.generic = 1;
  %specify the total numer of segment
  
  for idx = 1:length(params)
    param = params(idx);
    if param.cmd.generic==1
      
      fprintf('Processing segment %s (%i of %i)\n', param.day_seg, idx, length(params));
      
      % Determine which records you want to load:
      frames_fn = '';
      frames_fn = ct_filename_support(param,frames_fn,'frames');
      load(frames_fn); % GET NUMBER OF FRAMES HERE (look at variable named "frames")
      
      if isempty(param.cmd.frms)
        frms = 1:length(frames.frame_idxs);
      else
        frms = param.cmd.frms;
      end
      for frm = frms
        fprintf('  Frame %d of %d\n', frm, length(frms));
        
        if frm == length(frames.frame_idxs)
          records_fn = '';
          records_fn = ct_filename_support(param,records_fn,'records');
          load(records_fn);
          all_recs = [frames.frame_idxs(frm) length(records.lat)];
          clear records;
        else
          all_recs = [frames.frame_idxs(frm) frames.frame_idxs(frm+1)-1];
        end
        
        block_size = 8000;
        blocks = [0:block_size:diff(all_recs)+1];
        for block_idx = 1:length(blocks)
          for wf = 1:length(param.radar.wfs)
            for adc = 1:length(param.radar.wfs(wf).rx_paths)
              param.load_data.recs = all_recs(1) + blocks(block_idx) + [0 block_size-1];
              param.load_data.recs(2) = min(param.load_data.recs(2),all_recs(2));
              param.load_data.records_fn = param.records.records_fn;
              param.load_data.imgs = {[wf adc]};
              param.load_data.pulse_comp         = false;
              param.load_data.ft_dec             = false; % Fast-time decimation
              param.load_data.ft_wind            = @boxcar;
              param.load_data.ft_wind_time       = false; % Apply window on time domain chirp
              param.load_data.presums            = 1; % Coherent averaging
              param.load_data.combine_rx         = true;
              param.load_data.pulse_rfi.en       = false;
              param.load_data.pulse_rfi.inc_ave  = 101;
              param.load_data.pulse_rfi.thresh_scale = 10^(13/10);
              param.load_data.trim_vals          = [0 0];
              param.load_data.raw_data          = true;
              
              % Load data
              [data,hdr] = load_data(param);
              data = data{1};
              
              out_dir = fullfile(base_dir,param.season_name,param.day_seg);
              if ~exist(out_dir,'dir')
                mkdir(out_dir);
              end
              fn = fullfile(out_dir,sprintf('Data_%s_%03.0f_wf_%d_adc_%02d_block_%02d', param.day_seg, frm, wf, adc, block_idx));
              fprintf('  Saving data %s\n', fn);
              save(fn,'-v6','data','hdr');
              if 0
                % Plot data
                imagesc(lp(data));
                keyboard
              end
            end
          end
        end
      end
    end
  end
  
end