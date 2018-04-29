% Script collate_deconv.m
%
% This scripts takes the results from coh_noise_tracker spectral analysis
% and creates deconvolution waveforms from spectral targets. The spectral
% analysis has actually already created the deconv_H filter and this
% function groups these filters together using correlation statistics,
% finds the best one for each twtt, and then stores these best waveforms
% into a file. This is done per segment, but if a segment does not have
% any spectral targets for a particular twtt bin, waveforms from other
% days are used.
%
% Sections:
% 1. Debug section for visualizing results
% 2. Stage 1: Remove bad waveforms, group good ones based on statistical
%    similarity (but require that groups are contiguous in time). Also
%    put each waveform into an altitude bin.
% 3. Stage 2: Often segments will not have a good deconvolution waveform
%    for all twtt bins. This stage grabs good deconv waveforms from
%    neighboring segments to fill in for the missing twtts.
%
% Problems/Short-comings:
% 1. The same deconvolution waveform is used for the entire segment.
%    This means that falling edge sidelobes from internal reflections
%    tend to not be suppressed well.
% 2. Current application of deconvolution waveform assumes the surface
%    is at 0 m elevation (i.e. so that the surface does not have to be
%    tracked). This works okay for sea ice, but will not work for land ice.
%    This weakness is in load_fmcw_data.m and not this script
%    (coh_noise_tracker_task does track the surface for the twtt estimate).
% 3. ONLY WORKS ON SEA ICE WHERE SPECULAR TARGETS ARE COMMON
%
% Example:
%  See run_collate_deconv.m to run.
%
% Author: Jilu Li, John Paden

% For debugging, leave empty otherwise
day_seg_debug = ''; % Set to 'YYYYMMDD_SS' to debug one segment

%%
% =========================================================================

if stage_one_en
  %% STAGE ONE
  % =========================================================================
  % Loads all the specular_* files and creates tmp_* files that have the
  % poorly performing deconvolution waveforms removed and the good ones
  % grouped and averaged.
  % Outputs:
  %   final.deconv_H = Nt by Nx matrix of fast-time matched filters for deconvolution
  %     Nx is the number of detected leads that pass the metric thresholds.
  %     Assumption is that frequency domain window has been applied already,
  %     so this is a relative change from that initial window.
  %   final.deconv_gps_time = 1 by Nx vector containing GPS times of each lead
  %   final.deconv_twtt_min = 1 by Nx vector containing the minimum TWTT
  %      bin that this waveform should be used for;
  %   final.deconv_twtt_max = 1 by Nx vector with max TWTT bin.
  %   final.deconv_impulse_response = Nt by Nx (impulse response after
  %     deconvolution filter is applied to the deconv_sample range line)
  % =========================================================================
  deconv_DDC_Mt = {};
  deconv_H = {};
  deconv_elev = [];
  deconv_f0 = [];
  deconv_f1 = [];
  deconv_metric = [];
  deconv_gps_time = [];
  
  for param_idx = 1:length(params)
    %% Is this segment selected in the param spreadsheet
    param = params(param_idx);
    if isempty(day_seg_debug) && (~isfield(param.cmd,'generic') || iscell(param.cmd.generic) ...
        || ischar(param.cmd.generic) || ~param.cmd.generic) ...
        || (~isempty(day_seg_debug) && ~strcmp(day_seg_debug,param.day_seg))
      continue;
    end
    
    if iscell(param.analysis.specular.rbins)
      param.analysis.specular.rbins = param.analysis.specular.rbins{img};
    end
    
    if iscell(param.analysis.specular.Nt_shorten)
      param.analysis.specular.Nt_shorten = param.analysis.specular.Nt_shorten{img};
    end

    %% Load the specular surface file
    fprintf('Loading %s img %d wf_adc %d\n', param.day_seg, img, wf_adc);
    fn_dir = fileparts(ct_filename_out(param,spec_file_input_type, ''));
    fn = fullfile(fn_dir,sprintf('specular_%s_img_%02d.mat', param.day_seg, img));
    spec = load(fn);
    
    %% Create the frequency spectrum axis
    wf = spec.param_analysis.analysis.imgs{img}(1);
    if isempty(spec.deconv_mean)
      Nt = length(spec.wfs(wf).time);
    else
      Nt = mode(cellfun(@length,spec.deconv_mean)); % HACK: Force Nt to be constant... need to handle differently for multiple NZ in same processing block and DDC
    end
    [output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
    if strcmpi(radar_type,'deramp')
      spec.freq = spec.wf_freq;
      spec.Tpd = spec.wfs(wf).Tpd;
    end
    
    %% Handle the case where no specular targets were found
    if size(spec.deconv_mean,2) == 0
      warning('This segment has no deconvolution waveforms');
      final = [];
      final.param_collate = param;
      final.param_analysis = spec.param_analysis;
      final.metric = [];
      final.num_response = [];
      if any(strcmpi(ct_output_dir(radar_name),{'kuband','snow'}))
        final.match_freq = (spec.wfs(wf).fc + 2*spec.wfs(wf).chirp_rate / spec.wfs(wf).fs_raw * ((0:Nt-1) - floor(Nt/2))).';
        final.Tpd = spec.Tpd;
      end
      final.freq = {};
      final.deconv_DDC_Mt = [];
      final.deconv_H = {};
      final.deconv_gps_time = [];
      final.deconv_twtt_min = [];
      final.deconv_twtt_max = [];
      final.deconv_impulse_response = {};
      fn_out = fullfile(fn_dir,sprintf('deconv_tmp_%s.mat',spec.param_analysis.day_seg));
      save(fn_out,'-struct','final');
      continue;
    end
    
    % =====================================================================
    % =====================================================================
    % First Step:
    % Go through each set of specular range lines and extract a deconvolution
    % filter response and metrics.
    % =====================================================================
    % =====================================================================
    %% Preallocate Outputs from this first step
    num_rlines = size(spec.deconv_mean,2);
    spec.deconv_H = {};
    spec.deconv_raw = {};
    spec.rising_edge_SL = zeros(1,num_rlines);
    spec.falling_edge_SL = zeros(1,num_rlines);
    spec.rising_edge_ISL = zeros(1,num_rlines);
    spec.falling_edge_ISL = zeros(1,num_rlines);
    spec.width_ML = zeros(1,num_rlines);
    spec.peak = zeros(1,num_rlines);
    spec.metric = inf(6,num_rlines);
    mask = zeros(1,num_rlines);
    
    for rline = 1:num_rlines
      sig = spec.deconv_mean{rline};
      sig_std = spec.deconv_std{rline};
      sig_sample = spec.deconv_sample{rline};

      %% Move time zero to the center to make indexing easier
      sig = fftshift(sig,1);
      
      %% Estimate SNR as a function of range bin
      SNR = abs(sig).^2 ./ sig_std.^2;
      
      %% Apply windowing or time-gating to the specular target
      max_bin_search_bins = round(length(sig)/2) + (-50:50);
      [~,max_bin] = max(sig(max_bin_search_bins));
      max_bin = max_bin + max_bin_search_bins(1)-1;
      good_bins = max_bin + param.analysis.specular.rbins;
      sig_tg = zeros(size(sig));
      sig_tg(good_bins) = sig(good_bins) ...
        .* repmat(tukeywin_trim(length(good_bins),0.1), [1 size(sig,2)]);
      
      if any(strcmpi(radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3','snow8'}))
        %% Perform required FFT shifts to bring fc into center and time zero back to start
        sig_tg = ifft(fftshift(fft(ifftshift(sig_tg,1)),1));
        spec.freq{rline} = fftshift(spec.freq{rline});
      else
        %% Perform FFT shift to bring time zero back to start
        sig_tg = ifft(fftshift(fft(ifftshift(sig_tg,1)),1));
      end
      
      %% Extract the deconvolution information
      sig_tg_fft = fft(sig_tg);
      spec.deconv_H{rline} = zeros(size(sig_tg_fft));
      Nt_new = Nt-sum(param.analysis.specular.Nt_shorten);
      
      spec.deconv_H{rline}(param.analysis.specular.Nt_shorten(1)+(0:Nt_new-1)) ...
        = param.get_heights.ft_wind(Nt_new) ./ sig_tg_fft(param.analysis.specular.Nt_shorten(1)+(0:Nt_new-1));
      
      if debug_level == 3 || debug_level >= 0 && rline == 1
        %% DEBUG CODE: For setting Nt_shorten
        figure(1); clf;
        plot(lp(sig_tg_fft));
        title(sprintf('%s: Without Nt_shorten %.0f meters',param.day_seg,interp1(spec.gps_time,spec.elev,spec.deconv_gps_time(rline))),'interpreter','none');
        grid on;
        figure(2); clf;
        plot(lp(param.get_heights.ft_wind(Nt_new)));
        hold on
        plot(lp(sig_tg_fft(param.analysis.specular.Nt_shorten(1)+(0:Nt_new-1))), 'r');
        plot(lp(spec.deconv_H{rline}(param.analysis.specular.Nt_shorten(1)+(0:Nt_new-1))), 'g');
        hold off;
        title(sprintf('%s: With Nt_shorten',param.day_seg),'interpreter','none');
        legend('window','raw','correction','location','best');
        wf = spec.param_analysis.analysis.imgs{img}(wf_adc,1);
        adc = spec.param_analysis.analysis.imgs{img}(wf_adc,2);
        fig1_fn = [ct_filename_tmp(param,'','deconv','Nt_shorten') sprintf('_wf%d_adc%d_without.fig',wf,adc)];
        fig1_fn_dir = fileparts(fig1_fn);
        if ~exist(fig1_fn_dir)
          mkdir(fig1_fn_dir);
        end
        saveas(1,fig1_fn);
        fig2_fn = [ct_filename_tmp(param,'','deconv','Nt_shorten') sprintf('_wf%d_adc%d_with.fig',wf,adc)];
        saveas(2,fig2_fn);
        
        grid on;
        axis tight;
        debug_Nt_shorten_threshold = 35;
        max_val = max(lp(sig_tg_fft));
        Nt_shorten = find(lp(sig_tg_fft)>max_val-debug_Nt_shorten_threshold,1) - 1;
        Nt_shorten(2) = length(lp(sig_tg_fft)) - find(lp(sig_tg_fft)>max_val-debug_Nt_shorten_threshold,1,'last');
        
        fprintf('%s Nt_shorten\n\t%.0f\t%.0f\n', param.day_seg, Nt_shorten);
        
        if debug_level == 3
          fprintf('Usually set so deconv_H (green) does not have large values relative to zero where the signal is weak. Regions of the FFT waveform that are not stable from one deconv waveform to the next should be clipped if possible too. specular.interp_rbins can be used to interpolate across bad FFT bins in the middle of the waveform.\n');
          keyboard;
        end
      end
      
      %% Interpolating deconvolution spectrum where specified
      deconv_H = spec.deconv_H{rline};
      good_mask = logical(ones(size(spec.deconv_H{rline})));
      good_mask(param.analysis.specular.interp_rbins) = 0;
      interp_mag = abs(deconv_H);
      interp_mag(~good_mask) = interp1(find(good_mask),interp_mag(good_mask),find(~good_mask),'spline');
      interp_angle = angle(deconv_H);
      interp_angle = unwrap(angle(deconv_H));
      interp_angle(~good_mask) = interp1(find(good_mask),interp_angle(good_mask),find(~good_mask),'spline');
      
      spec.deconv_H{rline} = interp_mag .* exp(j*interp_angle);
      
      %% Apply deconvolution to sample range line
      sig_deconv = ifft(fftshift(fft(sig_sample),1) .* spec.deconv_H{rline});
      
      %% Normalize
      sample_Mt = lp(interpft(sig_sample,Mt*length(sig_deconv)));
      sample_peak = max(sample_Mt);
      sig_deconv_Mt = lp(interpft(ifft(ifftshift(fft(sig_deconv))),Mt*length(sig_deconv)));
      [sig_deconv_peak,peak_idx] = max(sig_deconv_Mt);
      rising_edge_bins = param.analysis.specular.rbins(1)*Mt : -param.analysis.specular.SL_guard_bins*Mt;
      falling_edge_bins = param.analysis.specular.SL_guard_bins*Mt : param.analysis.specular.rbins(end)*Mt;
      
      % Check to make sure peak is not too close to start/stop of range line
      if peak_idx + rising_edge_bins(1) < 1 || peak_idx+falling_edge_bins(end) > length(sig_deconv_Mt)
        spec.metric(:,rline) = NaN(6,1);
        continue;
      end
      
      % Determine the width of the main lobe
      rising_idx = peak_idx-1;
      while rising_idx >= 1 ...
          && sig_deconv_Mt(rising_idx) > sig_deconv_peak-param.analysis.specular.ML_threshold
        rising_idx = rising_idx - 1;
      end
      falling_idx = peak_idx+1;
      while falling_idx < length(sig_deconv_Mt) ...
          && sig_deconv_Mt(falling_idx) > sig_deconv_peak-param.analysis.specular.ML_threshold
        falling_idx = falling_idx + 1;
      end
      % We don't care so much about the falling edge... so we'll ignore it
      %spec.width_ML(rline) = (falling_idx - rising_idx) / Mt;
      spec.width_ML(rline) = 2*(peak_idx - rising_idx) / Mt;
      
      % Determine the rising and falling edge side lobes
      spec.rising_edge_SL(rline) = max(sig_deconv_Mt(peak_idx+rising_edge_bins)) - sig_deconv_peak;
      spec.falling_edge_SL(rline) = max(sig_deconv_Mt(peak_idx+falling_edge_bins)) - sig_deconv_peak;
      if 0
        figure(1); clf;
        plot(sig_deconv_Mt - sig_deconv_peak)
        hold on;
        plot(peak_idx+rising_edge_bins,sig_deconv_Mt(peak_idx+rising_edge_bins)-sig_deconv_peak,'r');
        xlims = peak_idx+rising_edge_bins;
        xlim(xlims([1 end]) + [-20 40]);
        grid on
        ylim([-60 0]);
        keyboard
      end
      
      % Skip bad waveforms if the rising/falling edge are at the edge of
      % the A-scope
      peak_idx = round(peak_idx/Mt);
      rising_idx = floor(rising_idx/Mt);
      falling_idx = ceil(falling_idx/Mt);
      if rising_idx + param.analysis.specular.rbins(1) < 1 ...
          || falling_idx + param.analysis.specular.rbins(end) >= length(sig_deconv)
        warning('waveform %d may not be a good one, skipped',rline);
        continue
      end

      % Calculate the integrated sidelobe level
      spec.rising_edge_ISL(rline) = sum(abs(sig_deconv(rising_idx + param.analysis.specular.rbins(1):rising_idx-1)).^2)/Mt;
      spec.rising_edge_ISL(rline) = lp(spec.rising_edge_ISL(rline))-sig_deconv_peak;
      spec.falling_edge_ISL(rline) = sum(abs(sig_deconv(falling_idx+1:falling_idx + param.analysis.specular.rbins(end))).^2)/Mt;
      spec.falling_edge_ISL(rline) = lp(spec.falling_edge_ISL(rline))-sig_deconv_peak;
      normal_factor = sample_peak - sig_deconv_peak;
      deconv_H = deconv_H * 10^(normal_factor/20);
      spec.deconv_H{rline} = spec.deconv_H{rline} * 10^(normal_factor/20);
      sig_deconv = sig_deconv * 10^(normal_factor/20);
      
      spec.peak(rline) = sample_peak;
      spec.deconv_raw{rline} = sig_deconv;
      
      % Define spec.metric so that lower is better
      spec.metric(:,rline) = [-spec.peak(rline) spec.width_ML(rline) spec.falling_edge_SL(rline) ...
        spec.rising_edge_SL(rline),spec.falling_edge_ISL(rline), spec.rising_edge_ISL(rline)];
      
      % Set the absolute metric levels that will be used to threshold
      % deconvolution waveforms as good or bad
      abs_metric = param.analysis.specular.abs_metric;
      %abs_metric = [-4 5.375 -27 -32 inf inf]; % DEBUG LINE FOR OVERRIDING PARAM SPREADSHEET
      
      %% DEBUG CODE
      if all(spec.metric(:,rline) <= abs_metric.')
        mask(rline) = 1;
        if debug_level > 2
          
          rline
          
          figure(1); clf;
          subplot(2,1,1);
          plot(lp(deconv_H),'r')
          hold on;
          time_gate_window = zeros(size(sig_deconv));
          [~,time_gate_bins] = max(sig_deconv);
          time_gate_bins = time_gate_bins + param.analysis.specular.rbins;
          time_gate_sig_deconv = zeros(size(sig_deconv));
          time_gate_sig_deconv(time_gate_bins) = sig_deconv(time_gate_bins) ...
            .* tukeywin_trim(length(param.analysis.specular.rbins),0.2);
          sig_deconv_fft = fft(time_gate_sig_deconv);
          plot(lp(sig_deconv_fft),'g')
          hold off;
          grid on;
          ylabel('amplitude')
          legend('inverse deconvolution filter','deconvolved sample signal','location','best')
          h_axes = gca;
          subplot(2,1,2);
          plot(angle(spec.deconv_H{rline}),'r')
          grid on;
          xlabel('freq index')
          ylabel('phase(rad)')
          ylim([-pi pi]);
          h_axes(2) = gca;
          linkaxes(h_axes,'x');
          
          %% Extract metrics from deconvolved sample range line
          
          [max_val,max_bin] = max(sig_sample);
          good_bins = max_bin + param.analysis.specular.rbins;
          figure(2); clf;
          plot((1:length(sig_deconv)) - max_bin, lp(sig_deconv)-lp(max_val));
          hold on
          plot((1:length(sig_deconv)) - max_bin, lp(sig_sample)-lp(max_val),'r');
          plot((1:length(sig_deconv)) - max_bin, lp(+max_val)+circshift(lp(sig_tg),[max_bin-1 0])-lp(max_val),'g')
          hold off;
          grid on;
          xlim(good_bins([1 end]) - max_bin)
          xlabel('relative range bin')
          ylabel('power(dB)')
          legend('deconvolved ice lead signal','ice lead signal','averaged ice lead signal','location','best')
          ylim([-80 0]);
          keyboard
        end
      end
      
    end
    
    if debug_level > 0
      %% DEBUG CODE
      mask = logical(mask);
      fprintf('Table 1. Minimum metric of all waveforms.\n');
      fprintf('%12s\t%12s\t%12s\t%12s\t%12s\t%12s\n', '-peak', 'ML_width', ...
        'Falling SL','Rising SL','Falling ISL','Rising ISL');
      spec.metric(1,:) = -spec.metric(1,:);
      fprintf('%12.1f\t%12.3f\t%12.1f\t%12.1f\t%12.1f\t%12.1f\n', ...
        min(spec.metric,[],2));
      spec.metric(1,:) = -spec.metric(1,:);
      fprintf('Table 2. Median metric of all waveforms.\n');
      fprintf('%12s\t%12s\t%12s\t%12s\t%12s\t%12s\n', '-peak', 'ML_width', ...
        'Falling SL','Rising SL','Falling ISL','Rising ISL');
      fprintf('%12.1f\t%12.3f\t%12.1f\t%12.1f\t%12.1f\t%12.1f\n', ...
        nanmedian(spec.metric,2));
      
      figure(1); clf;
      plot(spec.metric(1,mask),'k');
      hold on;
      plot(spec.metric(2,mask),'r');
      plot(spec.metric(3,mask),'g');
      plot(spec.metric(4,mask),'c');
      plot(spec.metric(5,mask),'b');
      plot(spec.metric(6,mask),'m');
      hold off;
      grid on;
      legend('P','ML','FSL','RSL','IFSL','IRSL')
      aa = gca;
      title(param.day_seg,'Interpreter','none')
      xlabel('Deconv waveform index');
      ylabel('Metric (lower is better)');
      
      % Aligns data to frames for helping to adjust parameters
      % -------------------------------------------------------
      % to find good deconvolution waveforms (e.g. specular.gps_times,
      % specular.metric).
      
      % Debug code which determines which frames each deconvolution waveform
      % comes from.
      records = load(ct_filename_support(spec.param_analysis,'','records'));
      load(ct_filename_support(spec.param_analysis,'','frames'));
      for idx = 1:length(spec.deconv_gps_time)
        spec.deconv_frame(idx) = find(spec.deconv_gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      end
      for idx = 1:length(spec.gps_time)
        spec.frame(idx) = find(spec.gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      end
      
      figure(2); clf;
      plot(spec.deconv_frame(mask),'.-');
      grid on;
      aa(2) = gca;
      linkaxes(aa,'x');
      title(param.day_seg,'Interpreter','none')
      grid on;
      xlabel('Deconv waveform index');
      ylabel('Frame');
      
      frame_time = gps_time_to_frame(records.gps_time,frames.frame_idxs);
      fig_h = 4; figure(fig_h); clf;
      plot(frame_time, records.elev)
      hold on;
      deconv_frame_time = interp1(records.gps_time,frame_time,spec.deconv_gps_time);
      plot(deconv_frame_time, interp1(frame_time,records.elev,deconv_frame_time),'rx','LineWidth',2);
      grid on;
      xlabel('Frame');
      ylabel('Elevation (m)');
      legend('Elev','Waveform','Location','best');
      title(param.day_seg,'Interpreter','none')
      
      mask_idxs = find(mask);
      if ~any(mask) || debug_level == 2
        warning('No waveforms passed. Relax specular.abs_metric requirements.\n');
        figure(1); clf;
        plot(spec.metric(1,:),'k');
        hold on;
        plot(spec.metric(2,:),'r');
        plot(spec.metric(3,:),'g');
        plot(spec.metric(4,:),'c');
        plot(spec.metric(5,:),'b');
        plot(spec.metric(6,:),'m');
        hold off;
        grid on;
        legend('P','ML','FSL','RSL','IFSL','IRSL')
        aa = gca;
        title(param.day_seg,'Interpreter','none')
        xlabel('Deconv waveform index');
        ylabel('Metric (lower is better)');

        
        figure(2); clf;
        plot(spec.deconv_frame(:),'.-');
        grid on;
        aa(2) = gca;
        linkaxes(aa,'x');
        title(param.day_seg,'Interpreter','none')
        xlabel('Deconv waveform index');
        ylabel('Frame');
        
        figure(3); clf;
        peakiness_frame_time = interp1(records.gps_time,frame_time,spec.gps_time);
        plot(peakiness_frame_time, spec.peakiness, '.-')
        xlabel('Frame');
        ylabel('Peakiness');
        title(param.day_seg,'Interpreter','none')
        
        [day_seg,frm_id,recs] = get_frame_id(param,spec.deconv_gps_time);

        % Nyquist zone twtt barriers
        % - For each group of twtt execute the relative metric
        % - For all groups force the absolute metric
        wfs = spec.param_analysis.radar.wfs;
        BW = (wfs.f1-wfs.f0)*wfs.fmult;
        chirp_rate = BW/wfs.Tpd;
        nz_twtt = abs(spec.param_analysis.radar.fs / chirp_rate / 2 / TWTT_GROUPS_PER_NZ);
        twtt_bin_spacing = nz_twtt;
        
        twtt_zone = 1 + floor(spec.deconv_twtt / nz_twtt);
        
        % Useful debug fprintf
        fprintf('Table 3. Metrics for each waveform where lower is better.\n');
        fprintf('%5s\t%5s\t%12s\t%12s\t%12s\t%12s\t%12s\t%12s\t%20s\t%12s\t%6s\n', 'Index','Frame', ...
          '-peak','ML_width','Falling SL','Rising SL','Falling ISL','Rising ISL','GPS time','Record','twtt');
        fprintf('%5.0f\t%5.0f\t%12.1f\t%12.3f\t%12.1f\t%12.1f\t%12.1f\t%12.1f%20.3f\t%12.0f\t%6.0f\n', ...
          [1:length(spec.deconv_frame); spec.deconv_frame; spec.metric; spec.deconv_gps_time; recs; twtt_zone]);
      end
      
      keyboard
    end
    
    metric = spec.metric;
    
    %% Final RDS Deconvolution Waveform Generation
    if any(strcmpi(output_dir,{'accum','rds'}))
      % Choose a reference function
      best_score = sum(spec.metric(2:4,:));
      [~,best_idx] = min(best_score);
      spec.metric(:,best_idx)
      
      wf = spec.param_analysis.analysis.imgs{img}(wf_adc,1);
      adc = spec.param_analysis.analysis.imgs{img}(wf_adc,2);
      
      % Prepare frequency axis
      Nt = spec.wfs(wf).Nt_pc;
      fc = spec.wfs(wf).fc;
      if spec.wfs(wf).DDC_mode == 0
        fs = param.radar.fs;
        dt = 1/fs;
        df = 1/(Nt*dt);
        freq = fs*floor(fc/fs) + (0:df:(Nt-1)*df).';
      else
        fs = param.radar.fs / 2^(1+spec.wfs(wf).DDC_mode);
        dt = 1/fs;
        df = 1/(Nt*dt);
        freq = spec.wfs(wf).DDC_freq + df*ifftshift( -floor(Nt/2) : floor((Nt-1)/2) ).';
      end

      % Start with the original reference function that was used to
      % compress the specular lead data
      ref = spec.wfs(wf).ref{adc};
      % Then add in the correction
      ref = ref .* interp1(fftshift(spec.wfs(wf).freq), spec.deconv_H{best_idx}, freq, 'linear', 0);
      
      % Estimate delay and phase shift caused by deconvolution process
      % relative to reference waveform
      
      deconv_test = ifft(spec.wfs(wf).ref{adc} .* conj(ref));
      ideal_test = ifft(spec.wfs(wf).ref{adc} .* conj(spec.wfs(wf).ref{adc}));
      Mt = 20;
      deconv_test = interpft(deconv_test,Mt*length(deconv_test));
      ideal_test = interpft(ideal_test,Mt*length(ideal_test));
      figure(1); clf;
      plot(lp(deconv_test));
      hold on;
      plot(lp(ideal_test));
      
      [~,idx_error] = max(deconv_test);
      idx_error = (idx_error-1)/Mt
      
      phase_error = angle(deconv_test(1))*180/pi
      
      freq_norm = ifftshift( -floor(Nt/2) : floor((Nt-1)/2) ).';
      ref = ref .* exp(-1i*2*pi*idx_error*freq_norm/Nt);
      
      deconv_test = ifft(spec.wfs(wf).ref{adc} .* conj(ref));
      ideal_test = ifft(spec.wfs(wf).ref{adc} .* conj(spec.wfs(wf).ref{adc}));
      Mt = 20;
      deconv_test = interpft(deconv_test,Mt*length(deconv_test));
      ideal_test = interpft(ideal_test,Mt*length(ideal_test));
      figure(1); clf;
      plot(lp(deconv_test));
      hold on;
      plot(lp(ideal_test));
      
      [~,idx_error] = max(deconv_test);
      idx_error = (idx_error-1)/Mt
      
      phase_error = angle(deconv_test(1))*180/pi
      
      ref = ref .* exp(1i*phase_error/180*pi);
      
      deconv_test = ifft(spec.wfs(wf).ref{adc} .* conj(ref));
      ideal_test = ifft(spec.wfs(wf).ref{adc} .* conj(spec.wfs(wf).ref{adc}));
      Mt = 20;
      deconv_test = interpft(deconv_test,Mt*length(deconv_test));
      ideal_test = interpft(ideal_test,Mt*length(ideal_test));
      figure(1); clf;
      plot(lp(deconv_test));
      hold on;
      plot(lp(ideal_test));
      
      [~,idx_error] = max(deconv_test);
      idx_error = (idx_error-1)/Mt % Should be zero
      
      phase_error = angle(deconv_test(1))*180/pi % Should be ~zero
      
      % Remove time delay corrections that were applied to the original
      % reference function (these will be reapplied at processing time).
      ref = ref .* exp(-1i*2*pi*freq*spec.param_analysis.radar.wfs(wf).Tsys(spec.param_analysis.radar.wfs(wf).rx_paths(adc)));
      ref = ref .* exp(-1i*2*pi*freq*spec.wfs(wf).time_correction);
      ref = ifft(conj(ref));
      
      [max_val,max_bin] = max(spec.deconv_sample{best_idx});
      good_bins = max_bin + param.analysis.specular.rbins;
      figure(2); clf;
      tmp = lp( ifft( fft(spec.deconv_sample{best_idx}) .* ifftshift(spec.deconv_H{best_idx}) ) );
      plot(tmp-max(tmp));
      hold on
      tmp = lp(spec.deconv_sample{best_idx});
      plot(tmp-max(tmp),'r');
      tmp = lp(+max_val)+circshift(lp(spec.deconv_mean{best_idx}),[max_bin-1 0]);
      plot(tmp-max(tmp),'g')
      hold off;
      grid on;
      xlim(good_bins([1 end]))
      xlabel('range bin')
      ylabel('power(dB)')
      legend('deconvolved ice lead signal','ice lead signal','averaged ice lead signal','location','best')
      title(sprintf('Impulse Response %d to %d range bins',param.analysis.specular.rbins([1 end])));
      
      figure(1); clf;
      plot(lp(ref));
      grid on;
      xlabel('range bin')
      ylabel('power(dB)')
      title(sprintf('Inverse Filter Impulse Response %d to %d range bins', ...
        param.analysis.specular.ref_negative{img}(1), ...
        param.analysis.specular.ref_nonnegative{img}(end)));
      
      fig1_fn = [ct_filename_tmp(param,'','deconv','inverse_filter') sprintf('_wf%d_adc%d.fig',wf,adc)];
      fig1_fn_dir = fileparts(fig1_fn);
      if ~exist(fig1_fn_dir)
        mkdir(fig1_fn_dir);
      end
      saveas(1,fig1_fn);
      fig2_fn = [ct_filename_tmp(param,'','deconv','impulse_response') sprintf('_wf%d_adc%d.fig',wf,adc)];
      saveas(2,fig2_fn);
      
      if debug_level == 6
        keyboard
      end
      param_collate = param;
      ref_windowed = true;
      ref_window = param.get_heights.ft_wind;
      ref_nonnegative = ref(param.analysis.specular.ref_nonnegative{img});
      ref_negative = ref(param.analysis.specular.ref_negative{img} + end);

      param_collate = param;
      param_analysis = spec.param_analysis;

      records_fn = ct_filename_support(param,'','records');
      records = load(records_fn);
      rec = find(records.gps_time > spec.deconv_gps_time(best_idx),1);
      [~,board_idx] = adc_to_board(param.radar_name,adc);
      file_idx = find(records.relative_rec_num{board_idx} <= rec,1,'last');
      raw_fn = records.relative_filename{board_idx}{file_idx};
      fprintf('Best Raw File: %s\n', raw_fn);
      fprintf('UTC time: %s\n', datestr(epoch_to_datenum(gps_to_utc(spec.deconv_gps_time(best_idx)))))
      
      fn_dir = fileparts(ct_filename_out(param,spec_file_input_type, ''));
      fn = fullfile(fn_dir,sprintf('deconv_wf_%d_adc_%d_%s.mat', wf, adc, param.day_seg));
      fprintf('Saving %s img %d wf_adc %d: %s\n', param.day_seg, img, wf_adc, fn);
      save(fn,'ref_nonnegative','ref_negative','ref_windowed','ref_window','param_collate','best_idx','param_collate','param_analysis');
      continue;
    end
    
    %% Determine which deconv waveforms are good
    % Nyquist zone twtt barriers
    % - For each group of twtt execute the relative metric
    % - For all groups force the absolute metric
    wfs = spec.param_analysis.radar.wfs;
    BW = (wfs.f1-wfs.f0)*wfs.fmult;
    chirp_rate = BW/wfs.Tpd;
    nz_twtt = abs(spec.param_analysis.radar.fs / chirp_rate / 2 / TWTT_GROUPS_PER_NZ);
    twtt_bin_spacing = nz_twtt;
    
    %rel_metric = [1 1 1 0.06 1 1]; % For Debug: Set to 1 for no mask
    %abs_metric = [inf inf inf -23 inf inf]; % For Debug: Set to inf for no mask
    rel_metric = param.analysis.specular.rel_metric;
    abs_metric = param.analysis.specular.abs_metric;
    twtt_zone = 1 + floor(spec.deconv_twtt / nz_twtt);
    twtt_zones = unique(twtt_zone);
    
    mask = logical(ones(1,size(metric,2)));
    rel_thresh = [];
    for metric_idx = 1:size(metric,1)
      % Apply the absolute metric
      mask = mask & metric(metric_idx,:) <= abs_metric(metric_idx);
      % Apply the relative metric (but only compare similar ranges/twtts)
      % This chooses the best X% from among the good waveforms where X% is
      % defined by rel_metric.
      for twtt_zone_idx = 1:length(twtt_zones)
        twtt_mask = twtt_zone == twtt_zones(twtt_zone_idx) & mask;
        if any(twtt_mask)
          sorted = sort(metric(metric_idx,twtt_mask));
          rel_thresh(twtt_zone_idx,metric_idx) = sorted(1 + round(rel_metric(metric_idx) * (length(sorted)-1)));
          if ~isfinite(rel_thresh(twtt_zone_idx,metric_idx))
            % This happens for waveforms that were too close to the fast time gate edge
            rel_thresh(twtt_zone_idx,metric_idx) = inf;
          end
          mask(twtt_mask) = mask(twtt_mask) ...
            & metric(metric_idx,twtt_mask) <= rel_thresh(twtt_zone_idx,metric_idx);
        end
      end
    end
    
    % specular.gps_times: user forced this waveform to be used
    for forced_idxs = 1:length(param.analysis.specular.gps_times)
      [gps_times_offset,gps_times_idx] = min(abs(spec.deconv_gps_time - param.analysis.specular.gps_times(forced_idxs)));
      fprintf('Added %d with offset %.1f sec\n', gps_times_idx, gps_times_offset);
      mask(gps_times_idx) = true;
    end
    
    % spec.bad_gps_times: user forced these waveforms to not be used
    if isfield(param.analysis.specular,'bad_gps_times')
      for forced_idxs = 1:size(param.analysis.specular.bad_gps_times,1)
        bad_mask = spec.deconv_gps_time >= param.analysis.specular.bad_gps_times(forced_idxs,1) ...
          & spec.deconv_gps_time <= param.analysis.specular.bad_gps_times(forced_idxs,2);
        for bad_idx = find(bad_mask)
          fprintf('Removed waveform idx %d\n', bad_idx);
        end
        mask(bad_mask) = 0;
      end
    end
    
    if all(mask==0)
      warning('%d waveforms passed out of %d', sum(mask), length(mask));
      [~,sort_idxs] = sort(spec.metric(4,:));
      spec.metric(:,sort_idxs(1:min(size(spec.metric,2),10))).'
    else
      fprintf('%d waveforms passed out of %d\n', sum(mask), length(mask));
    end
    
    %% Create table of results
    mask_idxs = find(mask);
    table = [];
    if isempty(mask_idxs)
      able.deconv_impulse_response = {};
      table.freq = {};
      table.deconv_H = {};
    else
      for idx = 1:length(mask_idxs)
        table.deconv_impulse_response{idx} = spec.deconv_raw{mask_idxs(idx)};
        table.freq{idx} = spec.freq{mask_idxs(idx)};
        table.deconv_H{idx} = spec.deconv_H{mask_idxs(idx)};
      end
    end
    table.Tpd = spec.Tpd;
    table.deconv_gps_time = spec.deconv_gps_time(mask);
    table.twtt = spec.deconv_twtt(mask);
    table.deconv_DDC_Mt = spec.deconv_DDC_Mt(mask);
    metric = metric(:,mask);
    
    if 0
      %% DEBUG CODE: Aligns data to frames for helping to adjust parameters
      % to find good deconvolution waveforms (e.g. specular.gps_times,
      % specular.metric).
      
      % Debug code which determines which frames each deconvolution waveform
      % comes from.
      records = load(ct_filename_support(spec.param_analysis,'','records'));
      load(ct_filename_support(spec.param_analysis,'','frames'));
      for idx = 1:length(table.deconv_gps_time)
        table.deconv_frame(idx) = find(table.deconv_gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      end
      fig_h = 1; figure(fig_h); clf;
      plot(table.deconv_frame);
      h_axis(fig_h) = gca;
      grid on;
      
      % fig_h = 2; figure(fig_h); clf;
      % imagesc(lp(table.deconv_H));
      % h_axis(fig_h) = gca;
      % grid on;
      
      % fig_h = 3; figure(fig_h); clf;
      % deconv_H_angle = unwrap(angle(table.deconv_H));
      % deconv_H_angle = deconv_H_angle - repmat(deconv_H_angle(round(size(deconv_H_angle,1)/2),:),[size(deconv_H_angle,1) 1]);
      % imagesc(deconv_H_angle);
      % h_axis(fig_h) = gca;
      % grid on;
      
      linkaxes(h_axis,'x');
      
      keyboard
    end
    
    %% Group similar deconvolution waveforms, average them together
    % First group waveforms by number of samples(the DDC filter should be the same in this case,
    % then for those with the same number of samples filter, uses correlation statistics to group
    tmp = [];
    final = [];
    final.metric = [];
    final.num_response = [];
    final.deconv_H = {};
    final.deconv_gps_time = [];
    final.deconv_twtt_min = [];
    final.deconv_twtt_max = [];
    final.deconv_impulse_response = {};
    final.freq = {};
    final.deconv_DDC_Mt = [];
    num_samples = zeros(size(table.deconv_H));
    for idx = 1:length(num_samples)
      num_samples(idx) = length(table.deconv_H{idx});
    end
    num_sam_jumps = find(abs(diff(num_samples))>0);
    num_sam_jumps = [0 num_sam_jumps length(num_samples)];
    for idx = 1:length(num_sam_jumps)-1
      idxs = [num_sam_jumps(idx)+1:num_sam_jumps(idx+1)];
      tmp.deconv_gps_time = table.deconv_gps_time(idxs);
      tmp.twtt = table.twtt(idxs);
      tmp.deconv_DDC_Mt = table.deconv_DDC_Mt(idxs);
      if isempty(idxs) % cases for segments no deconvolution waveforms found
        break;
      else
        tmp.deconv_H = zeros(length(table.deconv_H{idxs(1)}),length(idxs));
        tmp.deconv_impulse_response = zeros(length(table.deconv_impulse_response{idxs(1)}),length(idxs));
        tmp.freq = zeros(length(table.freq{idxs(1)}),length(idxs));
        for idx1 = 1:length(idxs)
          tmp.deconv_impulse_response(:,idx1) = table.deconv_impulse_response{idxs(idx1)};
          tmp.freq(:,idx1) = table.freq{idxs(idx1)};
          tmp.deconv_H(:,idx1) = table.deconv_H{idxs(idx1)};
        end
      end
      done = zeros(size(tmp.twtt));
      next_idx = 1;
      deconv_corr_norm = sqrt(sum(tmp.deconv_H .* conj(tmp.deconv_H)));
      while any(~done)
        test_idx = next_idx;
        
        % Find the correlation of the current deconvolution function with
        % all other deconvolution functions
        deconv_corr = sum(tmp.deconv_H.*repmat(conj(tmp.deconv_H(:,test_idx)), [1 size(tmp.deconv_H,2)]) );
        
        % Normalize the correlation
        corr_metric = abs(deconv_corr ./ deconv_corr_norm);
        corr_metric = corr_metric/max(corr_metric);
        corr_metric = corr_metric .* exp(log(CORR_METRIC_THRESHOLD) ...
          *(tmp.deconv_gps_time - tmp.deconv_gps_time(test_idx)).^2 / CORR_METRIC_TIME_CONSTANT^2);
        
        if 0
          % Debug/analysis plots
          figure(1); clf;
          subplot(3,1,1);
          imagesc(angle(tmp.deconv_H))
          a1 = gca;
          subplot(3,1,2);
          plot(corr_metric)
          grid on;
          a2 = gca;
          subplot(3,1,3);
          plot(tmp.twtt)
          a3 = gca;
          linkaxes([a1 a2 a3],'x');
        end
        
        % Average all results that have similar correlation values...
        % - Still only group results that are close in time though... the system
        %   seems to jump between various common stability regimes and
        %   we need to keep track of the gps time when each regime occurs
        % - So the average is over all waveforms from the segment that have
        %   good correlation, but we still only time tag and group adjacent
        %   waveforms in time.
        
        H_mask = corr_metric >= CORR_METRIC_THRESHOLD;
        mask = logical(zeros(size(corr_metric)));
        next_idx = find(corr_metric(test_idx+1:end) < CORR_METRIC_THRESHOLD,1);
        if isempty(next_idx)
          next_idx = length(corr_metric)+1;
        else
          next_idx = test_idx + next_idx;
        end
        mask(test_idx:next_idx-1) = 1;
        
        % Remove lower performing deconvolution waveforms from the set
        metric_val = -metric(3,H_mask) - 4*metric(4,H_mask);
        metric_val = metric_val - min(metric_val);
        metric_val_idxs = find(H_mask);
        % Lower than 75% of the normalized best score gets dropped
        H_mask(metric_val_idxs) = metric_val >= max(metric_val)*0.75;
        
        % Collect all the results in an output structure "final"
        final.metric = cat(2,final.metric,mean(metric(:,H_mask),2));
        final.num_response = cat(2,final.num_response,sum(H_mask));
        final.deconv_H = cat(2,final.deconv_H,mean(tmp.deconv_H(:,H_mask),2));
        final.deconv_gps_time = cat(2,final.deconv_gps_time,mean(tmp.deconv_gps_time(mask)));
        final.deconv_DDC_Mt = cat(2,final.deconv_DDC_Mt,mean(tmp.deconv_DDC_Mt(mask)));
        % Grab the best deconv impulse response
        [~,best_idx] = max(-metric(3,H_mask) -2*metric(4,H_mask));
        H_mask_idxs = find(H_mask);
        best_idx = H_mask_idxs(best_idx);
        final.deconv_impulse_response ...
          = cat(2,final.deconv_impulse_response,tmp.deconv_impulse_response(:,best_idx));
        final.freq = cat(2,final.freq,tmp.freq(:,best_idx));
        
        twtt_bin = round(mean(tmp.twtt(mask))/twtt_bin_spacing);
        final.deconv_twtt_min = cat(2,final.deconv_twtt_min,twtt_bin-1);
        final.deconv_twtt_max = cat(2,final.deconv_twtt_max,twtt_bin+1);
        
        % Set this group of deconvolution waveforms to done
        done(mask) = true;
      end
    end
    final.match_freq = (spec.wfs(wf).fc + 2*spec.wfs(wf).chirp_rate / spec.wfs(wf).fs_raw * ((0:Nt-1) - floor(Nt/2))).';
    final.Tpd = spec.Tpd;
    final.param_analysis = spec.param_analysis;
    
    if isempty(final.deconv_H)
      warning('This segment has no deconvolution waveforms that pass metrics.');
      final.freq = {};
    end
    
    %% Save output
    final.param_collate = param;
    fn_out = fullfile(fn_dir,sprintf('deconv_tmp_%s.mat',spec.param_analysis.day_seg));
    fprintf('  Saving %s\n', fn_out);
    save(fn_out,'-struct','final');
    
    if 0
      %% DEBUG CODE: Aligns data to frames for helping to adjust parameters
      % to find good deconvolution waveforms (e.g. specular.gps_times,
      % specular.metric).
      
      % Debug code which determines which frames each deconvolution waveform
      % comes from.
      records = load(ct_filename_support(spec.param_analysis,'','records'));
      load(ct_filename_support(spec.param_analysis,'','frames'));
      for idx = 1:length(final.deconv_gps_time)
        final.deconv_frame(idx) = find(final.deconv_gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      end
      for idx = 1:length(spec.deconv_gps_time)
        spec.deconv_frame(idx) = find(spec.deconv_gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      end
      for idx = 1:length(spec.gps_time)
        spec.frame(idx) = find(spec.gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      end
      
      figure(1); clf;
      imagesc(lp(final.deconv_H(:,final.deconv_frame>=1 & final.deconv_frame<=inf)));
      
      figure(2); clf;
      subplot(2,1,1);
      plot(spec.peakiness);
      ylabel('Peakiness');
      a1 = gca;
      grid on;
      subplot(2,1,2);
      plot(spec.frame,'.');
      ylabel('Frame');
      a2 = gca;
      grid on;
      linkaxes([a1 a2],'x');
      
      keyboard
    end
    
  end
end

if stage_two_en
  %% STAGE TWO
  % =========================================================================
  % Loads all the tmp_* files and creates deconv_* files that have
  % waveforms added in for missing elevations.  A segment may have data
  % collected at an altitude or with a DDC filter where no good deconvolution
  % waveform was collected and this code looks in other segments for these missing waveforms.
  % =========================================================================

  [output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
  if any(strcmpi(output_dir,{'accum','rds'}))
    error('Stage two does not support accum or rds (not needed).\n');
  end
  
  deconv_H = [];
  deconv_elev = [];
  deconv_f0 = [];
  deconv_f1 = [];
  deconv_metric = [];
  deconv_gps_time = [];
  overall_min_twtt = inf;
  overall_max_twtt = -inf;
  overall_min_DDC_Mt = inf;
  overall_max_DDC_Mt = -inf;
  fns = {};
  twtt_table = {};
  freq_table = {};
  Tpd_table = {};
  DDC_Mt_table = {};
  for param_idx = 1:length(params)
    param = params(param_idx);
    
    %% Is this segment selected in the param spreadsheet
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
    
    %% Load the specular surface file
    fprintf('Loading %s\n', param.day_seg);
    fn_dir = fileparts(ct_filename_out(param,spec_file_input_type, ''));
    fn = fullfile(fn_dir,sprintf('deconv_tmp_%s.mat', param.day_seg));
    fns{end+1} = fn;
    if ~exist(fn)
      continue
    end
    final = load(fn);
    
    min_twtt = min(final.deconv_twtt_min);
    if overall_min_twtt > min_twtt
      overall_min_twtt = min_twtt;
    end
    max_twtt = max(final.deconv_twtt_max);
    if overall_max_twtt < max_twtt
      overall_max_twtt = max_twtt;
    end
    twtt_table{length(fns)} = [];
    for idx = 1:length(final.deconv_twtt_min)
      twtt_table{length(fns)} = unique([twtt_table{length(fns)}, final.deconv_twtt_min(idx):final.deconv_twtt_max(idx)]);
    end
    freq_table{length(fns)} = [min(cellfun(@(x) min(x(:)),final.freq)) ...
      max(cellfun(@(x) max(x(:)),final.freq))];
    Tpd_table{length(fns)} = final.Tpd;
    
    min_DDC_Mt = min(final.deconv_DDC_Mt);
    if overall_min_DDC_Mt > min_DDC_Mt
      overall_min_DDC_Mt = min_DDC_Mt;
    end
    max_DDC_Mt = max(final.deconv_DDC_Mt);
    if overall_max_DDC_Mt < max_DDC_Mt
      overall_max_DDC_Mt = max_DDC_Mt;
    end
    DDC_Mt_table{length(fns)} = unique([min(final.deconv_DDC_Mt) max(final.deconv_DDC_Mt)]);
    
    
  end
  
  overall_min_twtt
  overall_max_twtt
  overall_min_DDC_Mt
  overall_max_DDC_Mt
  
  for param_idx = 1:length(params)
    param = params(param_idx);
    
    %% Is this segment selected in the param spreadsheet
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) ...
        || ischar(param.cmd.generic) || ~param.cmd.generic ...
        || (~isempty(day_seg_debug) && ~strcmp(day_seg_debug,param.day_seg))
      continue;
    end
    
    %% Load the specular surface file
    fprintf('Loading %s\n', param.day_seg);
    fn_dir = fileparts(ct_filename_out(param,spec_file_input_type, ''));
    fn = fullfile(fn_dir,sprintf('deconv_tmp_%s.mat', param.day_seg));
    if ~exist(fn)
        final.metric = [];
        final.num_response = [];
        final.deconv_H = {};
        final.deconv_gps_time = [];
        final.deconv_twtt_min = [];
        final.deconv_twtt_max = [];
        final.deconv_impulse_response = {};
        final.freq = {};
        final.deconv_DDC_Mt = [];
        final.deconv_frame = [];
        final.param_collate = param;
    else
      final = load(fn);
    end
    
    % Create "twtts" the vector of all missing two way travel times
    if isempty(final.deconv_twtt_min)
      twtts = overall_min_twtt:overall_max_twtt;
    else
      min_twtt = min(final.deconv_twtt_min);
      max_twtt = max(final.deconv_twtt_max);
      twtts = [];
      if max_twtt < overall_max_twtt
        twtts = cat(2,twtts,max_twtt+1:overall_max_twtt);
      end
      if min_twtt > overall_min_twtt
        twtts = cat(2,twtts,overall_min_twtt:min_twtt-1);
      end
    end
    
    fn_idx = strmatch(fn,fns,'exact');
    for twtt = twtts
      % Find the closest file with a twtt that matches the missing time.
      % The start/stop frequencies and pulse duration also need to match.
      found = false;
      for fn_idx_offset = 1:length(fns)
        % Start looking at the closest segments first and then move out from
        % there.
        fn_idx_new = fn_idx+fn_idx_offset;
        if fn_idx_new <= length(fns)
          if any(twtt_table{fn_idx_new} == twtt) ...
              && abs(freq_table{fn_idx_new}(1) - min(final.match_freq)) < 0.1e9 ...
              && abs(freq_table{fn_idx_new}(2) - max(final.match_freq)) < 0.1e9 ...
              && abs(Tpd_table{fn_idx_new} - final.Tpd) < 1e-6
            found = true;
            break
          end
        end
        fn_idx_new = fn_idx-fn_idx_offset;
        if fn_idx_new >= 1
          if any(twtt_table{fn_idx_new} == twtt) ...
              && abs(freq_table{fn_idx_new}(1) - min(final.match_freq)) < 0.1e9 ...
              && abs(freq_table{fn_idx_new}(2) - max(final.match_freq)) < 0.1e9 ...
              && abs(Tpd_table{fn_idx_new} - final.Tpd) < 1e-6
            found = true;
            break
          end
        end
      end
      if ~found
        continue;
      end
      
      %% Add the new twtt to the current file
      final_new = load(fns{fn_idx_new});
      
      % Find the one with the best match
      valid_idxs = find(final_new.deconv_twtt_min <= twtt & final_new.deconv_twtt_max >= twtt);
      [~,best_idx] = min(final_new.metric(4,valid_idxs));
      best_idx = valid_idxs(best_idx);
      
      if any(final_new.deconv_gps_time(:,best_idx) == final.deconv_gps_time)
        continue;
      end
      final.metric = cat(2,final.metric,final_new.metric(:,best_idx));
      final.num_response = cat(2,final.num_response,final_new.num_response(:,best_idx));
      final.deconv_DDC_Mt = cat(2,final.deconv_DDC_Mt,final_new.deconv_DDC_Mt(:,best_idx));
      final.deconv_H = cat(2,final.deconv_H,final_new.deconv_H{best_idx});
      final.deconv_gps_time = cat(2,final.deconv_gps_time,final_new.deconv_gps_time(:,best_idx));
      final.deconv_twtt_min = cat(2,final.deconv_twtt_min,final_new.deconv_twtt_min(best_idx));
      final.deconv_twtt_max = cat(2,final.deconv_twtt_max,final_new.deconv_twtt_max(best_idx));
      final.deconv_impulse_response{end+1} = final_new.deconv_impulse_response{best_idx};
      final.freq{end+1} = final_new.freq{best_idx};
    end
    
    if isempty(final.deconv_DDC_Mt)
      DDC_Mts = overall_min_DDC_Mt:overall_max_DDC_Mt;
      DDC_Mts = intersect(DDC_Mts,[1 2 4 8 16]);
    else
      min_DDC_Mt = min(final.deconv_DDC_Mt);
      max_DDC_Mt = max(final.deconv_DDC_Mt);
      DDC_Mts = [];
      if max_DDC_Mt < overall_max_DDC_Mt
        DDC_Mts = cat(2,DDC_Mts,2.^([log2(max_DDC_Mt)+1:log2(overall_max_DDC_Mt)]));
      end
      if min_DDC_Mt > overall_min_DDC_Mt
        DDC_Mts = cat(2,DDC_Mts,2.^([log2(overall_min_DDC_Mt):log2(min_DDC_Mt)-1]));
      end
    end
    
    for DDC_Mt = DDC_Mts
      % Find the closest file with a missing deconv_DDC_Mt.
      found = false;
      for fn_idx_offset = 1:length(fns)
        % Start looking at the closest segments first and then move out from
        % there.
        fn_idx_new = fn_idx+fn_idx_offset;
        if fn_idx_new <= length(fns)
          if any(DDC_Mt_table{fn_idx_new} == DDC_Mt) ...
              && abs(freq_table{fn_idx_new}(1) - min(final.match_freq)) < 0.1e9 ...
              && abs(freq_table{fn_idx_new}(2) - max(final.match_freq)) < 0.1e9 ...
              && abs(Tpd_table{fn_idx_new} - final.Tpd) < 1e-6
            found = true;
            break
          end
        end
        fn_idx_new = fn_idx-fn_idx_offset;
        if fn_idx_new >= 1
          if any(DDC_Mt_table{fn_idx_new} == DDC_Mt) ...
              && abs(freq_table{fn_idx_new}(1) - min(final.match_freq)) < 0.1e9 ...
              && abs(freq_table{fn_idx_new}(2) - max(final.match_freq)) < 0.1e9 ...
              && abs(Tpd_table{fn_idx_new} - final.Tpd) < 1e-6
            found = true;
            break
          end
        end
      end
      if ~found
        continue;
      end
      
      %% Add the new DDC_Mt to the current file
      final_new = load(fns{fn_idx_new});
      
      % Find the one with the best match
      valid_idxs = find(final_new.deconv_DDC_Mt == DDC_Mt);
      [~,best_idx] = min(final_new.metric(4,valid_idxs));
      best_idx = valid_idxs(best_idx);
      
      final.metric = cat(2,final.metric,final_new.metric(:,best_idx));
      final.num_response = cat(2,final.num_response,final_new.num_response(:,best_idx));
      final.deconv_DDC_Mt = cat(2,final.deconv_DDC_Mt,final_new.deconv_DDC_Mt(:,best_idx));
      final.deconv_H = cat(2,final.deconv_H,final_new.deconv_H{best_idx});
      final.deconv_gps_time = cat(2,final.deconv_gps_time,final_new.deconv_gps_time(:,best_idx));
      final.deconv_twtt_min = cat(2,final.deconv_twtt_min,final_new.deconv_twtt_min(:,best_idx));
      final.deconv_twtt_max = cat(2,final.deconv_twtt_max,final_new.deconv_twtt_max(:,best_idx));
      final.deconv_impulse_response{end+1} = final_new.deconv_impulse_response{best_idx};
      final.freq{end+1} = final_new.freq{best_idx};
    end
    
    % Check to make sure there is at least one good waveform
    if isempty(final.deconv_impulse_response)
      warning('This segment has no deconvolution waveforms');
      keyboard
    end
    
    % Nyquist zone twtt barriers
    % - For each group of twtt execute the relative metric
    % - For all groups force the absolute metric
    wfs = final.param_analysis.radar.wfs;
    BW = (wfs.f1-wfs.f0)*wfs.fmult;
    chirp_rate = BW/wfs.Tpd;
    nz_twtt = abs(final.param_analysis.radar.fs / chirp_rate / 2 / TWTT_GROUPS_PER_NZ);
    final.twtt_bin_spacing = nz_twtt;
    
    %% Aligns data to frames for helping to adjust parameters
    % to find good deconvolution waveforms (e.g. specular.gps_times,
    % specular.metric).
    
    % Debug code which determines which frames each deconvolution waveform
    % comes from.
    records = load(ct_filename_support(param,'','records'));
    load(ct_filename_support(param,'','frames'));
    final.deconv_frame = NaN*size(final.deconv_gps_time);
    for idx = 1:length(final.deconv_gps_time)
      frame = find(final.deconv_gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      if isempty(frame) || final.deconv_gps_time(idx) > records.gps_time(end)
        final.deconv_frame(idx) = NaN; % This waveform is from another segment
      else
        final.deconv_frame(idx) = frame;
      end
    end
    
    if debug_level > 0
      
      clear h_axis;
      fig_h = 1; figure(fig_h); clf;
      plot(final.deconv_frame,'.-');
      h_axis(fig_h) = gca;
      grid on;
      title(param.day_seg,'Interpreter','none')
      xlabel('Deconv wf index');
      ylabel('Frame');
      
      if ~isempty(final.deconv_gps_time)
        % Interpolate cell contents to create a matrix for plotting
        [max_num_sam,max_num_sam_idx] = max(cellfun('length',final.deconv_H));
        deconv_H = zeros(max_num_sam,length(final.deconv_H));
        for rline=1:length(final.deconv_H)
          deconv_H(:,rline) = interp1(final.freq{rline},final.deconv_H{rline},final.freq{max_num_sam_idx});
        end
        fig_h = 2; figure(fig_h); clf;
        imagesc(lp(deconv_H));
        h_axis(fig_h) = gca;
        grid on;
        title(sprintf('%s (dB)', param.day_seg),'Interpreter','none')
        xlabel('Deconv wf index');
        ylabel('Frequency bin');
        
        fig_h = 3; figure(fig_h); clf;
        deconv_H_angle = unwrap(angle(deconv_H));
        deconv_H_angle = deconv_H_angle - repmat(deconv_H_angle(round(size(deconv_H_angle,1)/2),:),[size(deconv_H_angle,1) 1]);
        imagesc(deconv_H_angle);
        h_axis(fig_h) = gca;
        grid on;
        title(sprintf('%s (rad)', param.day_seg),'Interpreter','none')
        xlabel('Deconv wf index');
        ylabel('Frequency bin');
      end
      
      frame_time = gps_time_to_frame(records.gps_time,frames.frame_idxs);
      fig_h = 4; figure(fig_h); clf;
      plot(frame_time, records.elev)
      hold on;
      deconv_frame_time = interp1(records.gps_time,frame_time,final.deconv_gps_time);
      plot(deconv_frame_time, interp1(frame_time,records.elev,deconv_frame_time),'rx','LineWidth',2);
      grid on;
      xlabel('Frame');
      ylabel('Elevation (m)');
      legend('Elev','Waveform','Location','best');
      title(param.day_seg,'Interpreter','none')
      
      linkaxes(h_axis,'x');
      
      keyboard
    end
    
    [fn_dir,fn_name] = fileparts(fn);
    fn_out = fullfile(fn_dir,sprintf('deconv_%s.mat',param.day_seg));
    fprintf('  Saving %s\n', fn_out);
    if preserve_old
      final_old = load(fn_out);
      fprintf('\nPRESERVING OLD COLLATE DECONV MODE\n');
      fprintf('There are %d different waveform gps times\n', ...
        length(setdiff(final.deconv_gps_time,final_old.deconv_gps_time)));
      fprintf('There are %d same waveform gps times\n', ...
        length(intersect(final.deconv_gps_time,final_old.deconv_gps_time)));
      fprintf('EITHER RUN\n\tsave(fn_out,''-struct'',''final'');\n');
      fprintf('OR COPY AND PASTE CODE BELOW FOR "Add new gps times"\n');
      keyboard
      if 0
        % Add new gps times
        new_final = final;
        final = final_old;
        [~,new_idxs] = setdiff(new_final.deconv_gps_time,final.deconv_gps_time);
        
        final.metric = cat(2,final.metric,new_final.metric(:,new_idxs));
        final.num_response = cat(2,final.num_response,new_final.num_response(:,new_idxs));
        final.deconv_H = cat(2,final.deconv_H,new_final.deconv_H(new_idxs));
        final.deconv_gps_time = cat(2,final.deconv_gps_time,new_final.deconv_gps_time(:,new_idxs));
        final.deconv_twtt_min = cat(2,final.deconv_twtt_min,new_final.deconv_twtt_min(:,new_idxs));
        final.deconv_twtt_max = cat(2,final.deconv_twtt_max,new_final.deconv_twtt_max(:,new_idxs));
        final.deconv_impulse_response = cat(2,final.deconv_impulse_response,new_final.deconv_impulse_response(:,new_idxs));
        final.freq = cat(2,final.freq,new_final.freq(:,new_idxs));
        final.deconv_DDC_Mt = cat(2,final.deconv_DDC_Mt,new_final.deconv_DDC_Mt(:,new_idxs));
        final.deconv_frame = cat(2,final.deconv_frame,new_final.deconv_frame(:,new_idxs));
        final.param_collate = new_final.param_collate;
        
        % Update gps times if waveform has changed
        
        save(fn_out,'-struct','final');
      end
    else
      save(fn_out,'-struct','final');
    end
  end
end
