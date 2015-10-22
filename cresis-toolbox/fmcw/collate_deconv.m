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
% Author: John Paden

% param_fn = ct_filename_param('snow_param_2009_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2010_Greenland_DC8.xls');
% param_fn = ct_filename_param('snow_param_2011_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls');
param_fn = ct_filename_param('snow_param_2014_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2015_Greenland_C130.xls');
analysis_sheet = 'analysis_spec';

physical_constants;

stage_one_en = true;
CORR_METRIC_THRESHOLD = 0.996; % Found through experimentation
TWTT_GROUPS_PER_NZ = 5; % Number of two way travel time groups per Nyquist zone
Mt = 8; % Amount to over-sample when estimating peaks of lobes

stage_two_en = true; % It is important to enable all segments in the param sheet at once for this stage

spec_file_input_type = 'noise'; % e.g. set to 'noise' to input from CSARP_noise folder
spec_file_output_type = 'noise'; % e.g. set to 'noise' to output to CSARP_noise folder

debug_level = 0; % Set to zero to run with no plots/outputs/stops

preserve_old = false; % Set to true to not overwrite old deconv file

%% AUTOMATED SECTION
% =========================================================================
params = read_param_xls(param_fn,'',{analysis_sheet 'analysis'});
% For debugging, leave empty otherwise
day_seg_debug = '20140325_01'; % Set to 'YYYYMMDD_SS' to debug one segment
% day_seg_debug = '20100323_01'; % Set to 'YYYYMMDD_SS' to debug one segment
% day_seg_debug = '20150324_04'; % Set to 'YYYYMMDD_SS' to debug one segment

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
  
  deconv_H = [];
  deconv_elev = [];
  deconv_f0 = [];
  deconv_f1 = [];
  deconv_metric = [];
  deconv_gps_time = [];
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
    fn = fullfile(fn_dir,sprintf('specular_%s.mat', param.day_seg));
    spec = load(fn);
    
    %% Create the frequency spectrum axis
    wf = spec.param_analysis.analysis.imgs{1}(1);
    if isempty(spec.deconv_mean)
      Nt = length(spec.wfs(wf).time);
    else
      Nt = mode(cellfun(@length,spec.deconv_mean)); % HACK: Force Nt to be constant... need to handle differently for multiple NZ in same processing block and DDC
    end
    [output_dir,radar_type] = ct_output_dir(param.radar_name);
    if strcmpi(radar_type,'fmcw')
%       spec.freq = (spec.wfs(wf).fc + 2*spec.wfs(wf).chirp_rate / spec.wfs(wf).fs_raw * ((0:Nt-1) - floor(Nt/2))).';
      spec.freq = spec.wf_freq ;
      spec.Tpd = spec.wfs(wf).Tpd;
    end
    
    %% Handle the case where no specular targets were found
    if size(spec.deconv_mean,2) == 0
      warning('This segment has no deconvolution waveforms');
      final.freq = spec.freq;
      final.Tpd = spec.Tpd;
      final.deconv_H = [];
      final.deconv_gps_time = [];
      final.deconv_twtt_min = [];
      final.deconv_twtt_max = [];
      final.deconv_impulse_response = [];
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
    spec.deconv_H = zeros(Nt,num_rlines);
    spec.deconv_raw = zeros(Nt,num_rlines);
    spec.rising_edge_SL = zeros(1,num_rlines);
    spec.falling_edge_SL = zeros(1,num_rlines);
    spec.rising_edge_ISL = zeros(1,num_rlines);
    spec.falling_edge_ISL = zeros(1,num_rlines);
    spec.width_ML = zeros(1,num_rlines);
    spec.peak = zeros(1,num_rlines);
    spec.metric = inf(6,num_rlines);
    mask = zeros(1,num_rlines);
    
    for rline = 1:num_rlines
%       if rline == 2267 | rline == 2337 | rline == 5082 % for 20140317_04
%         continue
%       end
      sig = spec.deconv_mean{rline};
      sig_std = spec.deconv_std{rline};
      sig_sample = spec.deconv_sample{rline};
      
      if length(sig) ~= Nt
        % HACK: Ignore waveforms that don't have the most common length (need to handle for DDC)
        continue;
      end
      
      %% Move time zero to the center to make indexing easier
      sig = fftshift(sig,1);
      
      %% Estimate SNR as a function of range bin
      SNR = abs(sig).^2 ./ sig_std.^2;
      
      %% Apply windowing or time-gating to the specular target
      [~,max_bin] = max(sig);
      good_bins = max_bin + param.analysis.specular.rbins;
      sig_tg = zeros(size(sig));
      sig_tg(good_bins) = sig(good_bins) ...
        .* repmat(tukeywin_trim(length(good_bins),0.1), [1 size(sig,2)]);
      
      if any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
        %% Perform required FFT shifts to bring fc into center and time zero back to start
        sig_tg = ifft(fftshift(fft(ifftshift(sig_tg,1)),1));
      else
        %% Perform FFT shift to bring time zero back to start
        sig_tg = ifftshift(sig_tg,1);
      end
      
      %% Extract the deconvolution information
      sig_tg_fft = fft(sig_tg);
      spec.deconv_H(:,rline) = zeros(size(sig_tg_fft));
      Nt_new = Nt-sum(param.analysis.specular.Nt_shorten);
      
      if any(strcmpi(param.radar_name,{'snow','kuband','snow2','kuband2','snow3','kuband3'}))
        spec.deconv_H(param.analysis.specular.Nt_shorten(1)+(0:Nt_new-1),rline) ...
          = param.get_heights.ft_wind(Nt_new) ./ sig_tg_fft(param.analysis.specular.Nt_shorten(1)+(0:Nt_new-1));
      else
        H_window_shortened = hanning(Nt_new);
        
        H_window = zeros(Nt,1);
        H_window(1:floor(Nt+1)/2 - param.analysis.specular.Nt_shorten(1)) ...
          = H_window_shortened(floor(Nt+1)/2 - param.analysis.specular.Nt_shorten(1):-1:1);
        H_window(end:-1:floor(Nt+1)/2+1 + sum(param.analysis.specular.Nt_shorten)) ...
          = H_window_shortened(floor(Nt+1)/2+1 : end);
        
        deconv_H(:,end) = H_window ./ final2_fft;
        deconv_H(:,end) = deconv_H(:,end) ...
          / max(abs(fft(g_data(:,center_rline) .* deconv_H(:,end), size(g_data,1)*10))) ...
          * max(abs(fft(g_data(:,center_rline), size(g_data,1)*10)));
      end
      
      %% Interpolating deconvolution spectrum where specified
      deconv_H = spec.deconv_H(:,rline);
      good_mask = logical(ones(size(spec.deconv_H(:,rline))));
%       param.analysis.specular.interp_rbins = [4262:4712]; % DEBUG LINE
      good_mask(param.analysis.specular.interp_rbins) = 0;
      interp_mag = abs(deconv_H);
      interp_mag(~good_mask) = interp1(find(good_mask),interp_mag(good_mask),find(~good_mask),'spline');
      interp_angle = angle(deconv_H);
      interp_angle = unwrap(angle(deconv_H));
      interp_angle(~good_mask) = interp1(find(good_mask),interp_angle(good_mask),find(~good_mask),'spline');
      
      spec.deconv_H(:,rline) = interp_mag .* exp(j*interp_angle);
      
      %% Apply deconvolution to sample range line
      sig_deconv = ifft(fftshift(fft(sig_sample),1) .* spec.deconv_H(:,rline));
     
      %% Normalize
      sample_Mt = lp(ifft(fft(sig_sample),Mt*length(sig_deconv)));
      sample_peak = max(sample_Mt);
      sig_deconv_Mt = lp(ifft(fft(sig_deconv),Mt*length(sig_deconv)));
      [sig_deconv_peak,peak_idx] = max(sig_deconv_Mt);
      rising_edge_bins = param.analysis.specular.rbins(1)*Mt : -param.analysis.specular.SL_guard_bins*Mt;
      falling_edge_bins = param.analysis.specular.SL_guard_bins*Mt : param.analysis.specular.rbins(end)*Mt;
            
      % Check to make sure peak is not too close to start/stop of range line
      if peak_idx + rising_edge_bins(1) < 1 || peak_idx+falling_edge_bins(end) > length(sig_deconv_Mt)
        spec.metric(:,rline) = NaN(6,1);
        continue;
      end
      
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
      
      spec.width_ML(rline) = (falling_idx - rising_idx) / Mt;
      
      spec.rising_edge_SL(rline) = sig_deconv_peak - max(sig_deconv_Mt(peak_idx+rising_edge_bins));
      
      spec.falling_edge_SL(rline) = sig_deconv_peak - max(sig_deconv_Mt(peak_idx+falling_edge_bins));
    
      peak_idx = round(peak_idx/Mt);
      rising_idx = round(rising_idx/Mt); 
      falling_idx = round(falling_idx/Mt); 
      if rising_idx + param.analysis.specular.rbins(1) <1 | falling_idx + param.analysis.specular.rbins(end) >= length(sig_deconv)
        warning('waveform %d may not be a good one, skipped',rline);
        continue
      end
      spec.rising_edge_ISL(rline) = sum(abs(sig_deconv(rising_idx + param.analysis.specular.rbins(1):rising_idx)).^2);
      spec.rising_edge_ISL(rline) = lp(spec.rising_edge_ISL (rline))-sig_deconv_peak;
      spec.falling_edge_ISL(rline) = sum(abs(sig_deconv(falling_idx:falling_idx + param.analysis.specular.rbins(end))).^2);
      spec.falling_edge_ISL(rline) = lp(spec.falling_edge_ISL (rline))-sig_deconv_peak;
      normal_factor = sample_peak - sig_deconv_peak;
      deconv_H = deconv_H * 10^(normal_factor/20);
      spec.deconv_H(:,rline) = spec.deconv_H(:,rline) * 10^(normal_factor/20);
      sig_deconv = sig_deconv * 10^(normal_factor/20);
      
      spec.peak(rline) = sample_peak;
      spec.deconv_raw(:,rline) = sig_deconv;
      spec.metric(:,rline) = [-spec.peak(rline) spec.width_ML(rline) -spec.falling_edge_SL(rline) ...
        -spec.rising_edge_SL(rline),spec.falling_edge_ISL(rline), spec.rising_edge_ISL(rline)];

      %% DEBUG CODE
      abs_metric = param.analysis.specular.abs_metric;
      %abs_metric = [-4 5.375 -27 -32]; % DEBUG LINE FOR OVERRIDING PARAM SPREADSHEET
      if spec.peak(rline) >= -abs_metric(1) && spec.width_ML(rline) <= abs_metric(2) && spec.falling_edge_SL(rline) >= -abs_metric(3) ...
          && spec.rising_edge_SL(rline) >= -abs_metric(4) && spec.falling_edge_ISL(rline) <= abs_metric(5) && spec.rising_edge_ISL(rline) <= abs_metric(6)
        mask(rline) = 1;
        if debug_level > 1
          
          rline
          
          figure(1); clf;
          subplot(2,1,1);
          plot(lp(spec.deconv_H(:,rline)))
          hold on;
%           plot(lp(deconv_H),'r')
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
          legend('inverse deconvolution filter','deconvoled sample signal')
          h_axes = gca;
          subplot(2,1,2);
          plot(angle(spec.deconv_H(:,rline)))
          hold on;
%           plot(angle(deconv_H),'r')
          hold off;
          grid on;
          xlabel('freq index')
          ylabel('phase(rad)')
          h_axes(2) = gca;
          linkaxes(h_axes,'x');
          
          %% Extract metrics from deconvolved sample range line
          
          [max_val,max_bin] = max(sig_sample);
          good_bins = max_bin + param.analysis.specular.rbins;
          figure(2); clf;
          plot(lp(sig_deconv));
          hold on
          plot(lp(sig_sample),'r');
          plot(lp(+max_val)+circshift(lp(sig_tg),[max_bin-1 0]),'g')
          hold off;
          grid on;
          xlim(good_bins([1 end]))
          xlabel('range bin')
          ylabel('power(dB)')
          legend('deconvoled ice lead signal','ice lead signal','averaged ice lead signal')
          keyboard
        end
      end
      
    end
    
    if debug_level > 0
      %% DEBUG CODE
      mask = logical(mask);
      fprintf('%12.2f %12.2f %12.2f %12.2f\n', ...
        mean(spec.rising_edge_SL), ...
        mean(spec.falling_edge_SL), ...
        mean(spec.rising_edge_SL(mask)), ...
        mean(spec.falling_edge_SL(mask)));
      
      figure(1); clf;
      plot(spec.peak(mask),'k');
      hold on;
      plot(spec.width_ML(mask),'r');
      plot(spec.rising_edge_SL(mask),'g');
      plot(spec.falling_edge_SL(mask),'b');
      hold off;
      grid on;
      aa = gca;
      
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
      plot(spec.deconv_frame(mask));
      grid on;
      aa(2) = gca;
      linkaxes(aa,'x');
      
      mask_idxs = find(mask);
      % Useful debug fprintf:
      %fprintf('%5.0f %5.0f %6.1f %6.3f %6.1f %6.1f\n',[1:length(spec.deconv_frame); spec.deconv_frame; spec.metric])
      keyboard
    end
    
    metric = spec.metric;
    
    
    %% Determine which deconv waveforms are good
    % Nyquist zone twtt barriers
    % - For each group of twtt execute the relative metric
    % - For all groups force the absolute metric
    wfs = spec.param_analysis.radar.wfs;
    BW = (wfs.f1-wfs.f0)*wfs.fmult;
    chirp_rate = BW/wfs.Tpd;
    nz_twtt = abs(spec.param_analysis.radar.fs / chirp_rate / 2 / TWTT_GROUPS_PER_NZ);
    twtt_bin_spacing = nz_twtt;
    
    rel_metric = [1 1 1 0.06 1 1]; % Set to 1 for no mask
    abs_metric = [inf inf inf -23]; % Set to inf for no mask
    rel_metric = param.analysis.specular.rel_metric;
    abs_metric = param.analysis.specular.abs_metric;
    twtt_zone = 1 + floor(spec.deconv_twtt / nz_twtt);
    twtt_zones = unique(twtt_zone);
    
    mask = logical(ones(1,size(metric,2)));
    rel_thresh = [];
    for metric_idx = 1:size(metric,1)
      mask = mask & metric(metric_idx,:) <= abs_metric(metric_idx);
      for twtt_zone_idx = 1:length(twtt_zones)
        twtt_mask = twtt_zone == twtt_zones(twtt_zone_idx);
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
    
    % spec.gps_times: user forced this waveform to be used
    for forced_idxs = 1:length(param.analysis.specular.gps_times)
      [gps_times_offset,gps_times_idx] = min(abs(spec.deconv_gps_time - param.analysis.specular.gps_times(forced_idxs)));
      fprintf('Added %d with offset %.1f sec\n', gps_times_idx, gps_times_offset);
      mask(gps_times_idx) = true;
    end
    
    if all(mask==0)
      warning('%d waveforms passed out of %d', sum(mask), length(mask));
      [~,sort_idxs] = sort(spec.metric(4,:));
      spec.metric(:,sort_idxs(1:min(end,10))).'
    else
      fprintf('%d waveforms passed out of %d\n', sum(mask), length(mask));
    end
    
    %% Hacks for specific segments
    mask_idxs = find(mask);
%     if strcmpi(spec.param_analysis.radar_name,'snow2') && strcmpi(spec.param_analysis.day_seg,'20120315_03')
%       mask(mask_idxs([73 75 78])) = 0;
%     end
    
    table = [];
    table.freq = spec.freq(:,mask);
    table.Tpd = spec.Tpd;
    table.deconv_gps_time = spec.deconv_gps_time(mask);
    table.deconv_H = spec.deconv_H(:,mask);
    table.twtt = spec.deconv_twtt(mask);
    table.deconv_impulse_response = spec.deconv_raw(:,mask);
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
      
      fig_h = 2; figure(fig_h); clf;
      imagesc(lp(table.deconv_H));
      h_axis(fig_h) = gca;
      grid on;
      
      fig_h = 3; figure(fig_h); clf;
      deconv_H_angle = unwrap(angle(table.deconv_H));
      deconv_H_angle = deconv_H_angle - repmat(deconv_H_angle(round(size(deconv_H_angle,1)/2),:),[size(deconv_H_angle,1) 1]);
      imagesc(deconv_H_angle);
      h_axis(fig_h) = gca;
      grid on;
      
      linkaxes(h_axis,'x');
      
      keyboard
    end

    %% Group similar deconvolution waveforms, average them together
    % Uses correlation statistics to group
    done = zeros(size(table.twtt));
    tmp = table;
    final = [];
    final.metric = [];
    final.num_response = [];
    final.deconv_H = [];
    final.deconv_gps_time = [];
    final.deconv_twtt_min = [];
    final.deconv_twtt_max = [];
    final.deconv_impulse_response = [];
    final.freq = [];
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
      metric_val = -metric(3,H_mask) - 2*metric(4,H_mask);
      metric_val = metric_val - min(metric_val);
      metric_val_idxs = find(H_mask);
      % Lower than 75% of the normalized best score gets dropped
      H_mask(metric_val_idxs) = metric_val >= max(metric_val)*0.75;
      
      % Collect all the results in an output structure "final"
      final.metric = cat(2,final.metric,mean(metric(:,H_mask),2));
      final.num_response = cat(2,final.num_response,sum(H_mask));
      final.deconv_H = cat(2,final.deconv_H,mean(tmp.deconv_H(:,H_mask),2));
      final.deconv_gps_time = cat(2,final.deconv_gps_time,mean(tmp.deconv_gps_time(mask)));
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
%     final.freq = spec.freq;
    final.Tpd = spec.Tpd;
    final.param_analysis = spec.param_analysis;

    %% HACK: Required to remove bad waveforms that the code did not automatically find
%     if strcmpi(final.param_analysis.radar_name,'snow') && strcmpi(final.param_analysis.day_seg,'20110317_01')
%       warning('Use this if you ever needed to hand remove certain results, changes to code or parameters will break this removal!');
%       keyboard
%       bad_mask = logical(zeros(size(final.num_response)));
%       bad_mask(end-4:end-3) = 1; % <--- SET BAD DECONV WAVEFORMS HERE
%       final.metric = final.metric(:,~bad_mask);
%       final.num_response = final.num_response(:,~bad_mask);
%       final.deconv_H = final.deconv_H(:,~bad_mask);
%       final.deconv_gps_time = final.deconv_gps_time(:,~bad_mask);
%       final.deconv_twtt_min = final.deconv_twtt_min(:,~bad_mask);
%       final.deconv_twtt_max = final.deconv_twtt_max(:,~bad_mask);
%     end
%     if strcmpi(final.param_analysis.radar_name,'snow') && strcmpi(final.param_analysis.day_seg,'20110322_01')
%       warning('Use this if you ever needed to hand remove certain results, changes to code or parameters will break this removal!');
%       keyboard
%       bad_mask = logical(zeros(size(final.num_response)));
%       bad_mask([2:4,6:8]) = 1; % <--- SET BAD DECONV WAVEFORMS HERE
%       final.metric = final.metric(:,~bad_mask);
%       final.num_response = final.num_response(:,~bad_mask);
%       final.deconv_H = final.deconv_H(:,~bad_mask);
%       final.deconv_gps_time = final.deconv_gps_time(:,~bad_mask);
%       final.deconv_twtt_min = final.deconv_twtt_min(:,~bad_mask);
%       final.deconv_twtt_max = final.deconv_twtt_max(:,~bad_mask);
%     end
%     if strcmpi(final.param_analysis.radar_name,'snow2') && strcmpi(final.param_analysis.day_seg,'20120327_01')
%       warning('Use this if you ever needed to hand remove certain results, changes to code or parameters will break this removal!');
%       keyboard
%       bad_mask = logical(zeros(size(final.num_response)));
%       bad_mask([3:14]) = 1; % <--- SET BAD DECONV WAVEFORMS HERE
%       final.metric = final.metric(:,~bad_mask);
%       final.num_response = final.num_response(:,~bad_mask);
%       final.deconv_H = final.deconv_H(:,~bad_mask);
%       final.deconv_gps_time = final.deconv_gps_time(:,~bad_mask);
%       final.deconv_twtt_min = final.deconv_twtt_min(:,~bad_mask);
%       final.deconv_twtt_max = final.deconv_twtt_max(:,~bad_mask);
%     end
    
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
  % collected at an altitude where no good deconvolution waveform was
  % collected and this code looks in other segments for these missing waveforms.
  % =========================================================================
  
  deconv_H = [];
  deconv_elev = [];
  deconv_f0 = [];
  deconv_f1 = [];
  deconv_metric = [];
  deconv_gps_time = [];
  overall_min_twtt = inf;
  overall_max_twtt = -inf;
  fns = {};
  twtt_table = {};
  freq_table = {};
  Tpd_table = {};
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
    final = load(fn);
    
    min_twtt = min(final.deconv_twtt_min);
    if overall_min_twtt > min_twtt
      overall_min_twtt = min_twtt;
    end
    max_twtt = max(final.deconv_twtt_max);
    if overall_max_twtt < max_twtt
      overall_max_twtt = max_twtt;
    end
    
    twtt_table{length(fns)} = unique([final.deconv_twtt_min final.deconv_twtt_max]);
    freq_table{length(fns)} = [min(final.freq) max(final.freq)];
    Tpd_table{length(fns)} = final.Tpd;
    
    if min_twtt > 4
    end
    if max_twtt < 21
    end
  end
  
  overall_min_twtt
  overall_max_twtt
  
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
    final = load(fn);
    % Change to cell formatting to allow different lengths for deconv_impulse_response
    final.deconv_impulse_response = {final.deconv_impulse_response};

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
      % The start/stop frequencies also need to match.
      found = false;
      for fn_idx_offset = 1:length(fns)
        % Start looking at the closest segments first and then move out from
        % there.
        fn_idx_new = fn_idx+fn_idx_offset;
        if fn_idx_new <= length(fns)
          if any(twtt_table{fn_idx_new} == twtt) ...
              && abs(freq_table{fn_idx_new}(1) - min(final.freq)) < 0.1e9 ...
              && abs(freq_table{fn_idx_new}(2) - max(final.freq)) < 0.1e9 ...
              && abs(Tpd_table{fn_idx_new} - final.Tpd) < 1e-6
            found = true;
            break
          end
        end
        fn_idx_new = fn_idx-fn_idx_offset;
        if fn_idx_new >= 1
          if any(twtt_table{fn_idx_new} == twtt) ...
              && abs(freq_table{fn_idx_new}(1) - min(final.freq)) < 0.1e9 ...
              && abs(freq_table{fn_idx_new}(2) - max(final.freq)) < 0.1e9 ...
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
      
      final.metric = cat(2,final.metric,final_new.metric(:,best_idx));
      final.num_response = cat(2,final.num_response,final_new.num_response(:,best_idx));
      final.deconv_H = cat(2,final.deconv_H,interp1(final_new.freq,final_new.deconv_H(:,best_idx),final.freq));
      final.deconv_H(:,end) = interp_finite(final.deconv_H(:,end));
      final.deconv_gps_time = cat(2,final.deconv_gps_time,final_new.deconv_gps_time(:,best_idx));
      final.deconv_twtt_min = cat(2,final.deconv_twtt_min,twtt);
      final.deconv_twtt_max = cat(2,final.deconv_twtt_max,twtt);
      final.deconv_impulse_response{end+1} = final_new.deconv_impulse_response(:,best_idx);
    end
    final.twtt_bin_spacing = twtt_bin_spacing;
    
    %% Aligns data to frames for helping to adjust parameters
    % to find good deconvolution waveforms (e.g. specular.gps_times,
    % specular.metric).
    
    % Debug code which determines which frames each deconvolution waveform
    % comes from.
    records = load(ct_filename_support(param,'','records'));
    load(ct_filename_support(param,'','frames'));
    for idx = 1:length(final.deconv_gps_time)
      frame = find(final.deconv_gps_time(idx) >= records.gps_time(frames.frame_idxs),1,'last');
      if isempty(frame) || final.deconv_gps_time(idx) > records.gps_time(end)
        final.deconv_frame(idx) = NaN; % This waveform is from another segment
      else
        final.deconv_frame(idx) = frame;
      end
    end
    
    if debug_level > 0

      fig_h = 1; figure(fig_h); clf;
      plot(final.deconv_frame);
      h_axis(fig_h) = gca;
      grid on;
      
      fig_h = 2; figure(fig_h); clf;
      imagesc(lp(final.deconv_H));
      h_axis(fig_h) = gca;
      grid on;
      
      fig_h = 3; figure(fig_h); clf;
      deconv_H_angle = unwrap(angle(final.deconv_H));
      deconv_H_angle = deconv_H_angle - repmat(deconv_H_angle(round(size(deconv_H_angle,1)/2),:),[size(deconv_H_angle,1) 1]);
      imagesc(deconv_H_angle);
      h_axis(fig_h) = gca;
      grid on;
      
      linkaxes(h_axis,'x');
      
      keyboard
    end
    
    [fn_dir,fn_name] = fileparts(fn);
    fn_out = fullfile(fn_dir,sprintf('deconv_%s.mat',param.day_seg));
    fprintf('  Saving %s\n', fn_out);
    if preserve_old
      final_old = load(fn_out);
      keyboard
    else
      save(fn_out,'-struct','final');
    end
  end
end
