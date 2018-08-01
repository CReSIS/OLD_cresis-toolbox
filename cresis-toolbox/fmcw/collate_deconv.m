%function collate_deconv(param,param_override)
% collate_deconv(param,param_override)
%
% This scripts takes the results from analysis cmd spectral
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

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

physical_constants;

%% Input checks
% =====================================================================

cmd = param.analysis.cmd{param.collate_deconv.cmd_idx};
debug_level = param.collate_deconv.debug_level;

if ~isfield(param.collate_deconv,'imgs') || isempty(param.collate_deconv.imgs)
  param.collate_deconv.imgs = 1:length(param.analysis.imgs);
end

if ~isfield(param.collate_deconv,'wf_adcs') || isempty(param.collate_deconv.wf_adcs)
  param.collate_deconv.wf_adcs = [];
end

if ~isfield(param.collate_deconv,'gps_time_penalty') || isempty(param.collate_deconv.gps_time_penalty)
  param.collate_deconv.gps_time_penalty = 1/7200;
end

if ~isfield(param.collate_deconv,'in_dir') || isempty(param.collate_deconv.in_dir)
  param.collate_deconv.in_dir = 'analysis';
end

if ~isfield(param.collate_deconv,'min_score') || isempty(param.collate_deconv.min_score)
  param.collate_deconv.min_score = -10;
end

if ~isfield(param.collate_deconv,'Mt') || isempty(param.collate_deconv.Mt)
  param.collate_deconv.Mt = 10;
end
Mt = param.collate_deconv.Mt;

if ~isfield(param.collate_deconv,'out_dir') || isempty(param.collate_deconv.out_dir)
  param.collate_deconv.out_dir = 'analysis';
end

if ~isfield(param.collate_deconv,'preserve_old') || isempty(param.collate_deconv.preserve_old)
  param.collate_deconv.preserve_old = false;
end

if ~isfield(param.collate_deconv,'stage_one_en') || isempty(param.collate_deconv.stage_one_en)
  param.collate_deconv.stage_one_en = true;
end

if ~isfield(param.collate_deconv,'stage_two_en') || isempty(param.collate_deconv.stage_two_en)
  param.collate_deconv.stage_two_en = true;
end

if ~isfield(param.collate_deconv,'twtt_penalty') || isempty(param.collate_deconv.twtt_penalty)
  param.collate_deconv.twtt_penalty = 1e6;
end

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);
if ~strcmpi(radar_type,'deramp')
  % Only deramp radars need to run stage 2
  cmd.stage_two_en = false;
end

if ~isfield(cmd,'abs_metric') || isempty(cmd.abs_metric)
  error('The "abs_metric" field must be set in the cmd.');
end

if ~isfield(cmd,'bad_gps_times') || isempty(cmd.bad_gps_times)
  cmd.bad_gps_times = [];
end

if ~isfield(cmd,'f0') || isempty(cmd.f0)
  % Default is no limits on lower frequency
  cmd.f0 = -inf;
end

if ~isfield(cmd,'f1') || isempty(cmd.f1)
  % Default is no limits on upper frequency
  cmd.f1 = inf;
end

if ~isfield(cmd,'gps_times') || isempty(cmd.gps_times)
  cmd.gps_times = [];
end

if ~isfield(cmd,'interp_rbins') || isempty(cmd.interp_rbins)
  cmd.interp_rbins = [];
end

if ~isfield(cmd,'metric_weights') || isempty(cmd.metric_weights)
  cmd.metric_weights = [0.5 0 3 5 0 0];
end

if ~isfield(cmd,'ML_threshold') || isempty(cmd.ML_threshold)
  cmd.ML_threshold = 15;
end

if ~isfield(cmd,'rbins') || isempty(cmd.rbins)
  error('The "rbins" field must be set in the cmd to a range of indices about the peak to use in the deconvolution waveform, e.g. cmd.rbins = [-110 40] to use 110 bins before the peak and 40 bins after the peak. rbins should be a cell array with each element corresponding to the param.collate_deconv.imgs array.');
end

if ~isfield(cmd,'SL_guard_bins') || isempty(cmd.SL_guard_bins)
  cmd.SL_guard_bins = 3;
end

%% Stage 1
% =========================================================================
% Loads all the specular_* files and creates deconv_tmp_* files that have the
% poorly performing deconvolution waveforms removed and the good ones
% grouped and averaged.
% Outputs:
%  deconv.gps_time: gps-time of the impulse reponse, 1 by Nx
%  deconv.lat = 1 by Nx vector containing the latitude of this response
%  deconv.lon = 1 by Nx vector containing the longitude of this response
%  deconv.elev = 1 by Nx vector containing the elevation of this response
%  deconv.roll = 1 by Nx vector containing the roll of this response
%  deconv.pitch = 1 by Nx vector containing the pitch of this response
%  deconv.heading = 1 by Nx vector containing the heading of this response
%  deconv.frm: frame containing impulse response, 1 by Nx
%  deconv.rec: record containing impulse response, 1 by Nx
%  deconv.ref_windowed: logical indicating if reference impulse responses
%    are windowed
%  deconv.ref_window: reference window function handle
%  deconv.ref_nonnegative: impulse responses for non-negative time, Nt_nonneg by Nx
%  deconv.ref_negative: impulse responses for negative time, Nt_neg by Nx
%  deconv.ref_mult_factor: multiplication factor for deconv wf, 1 by Nx
%  deconv.impulse_reponse: deconvolved impulse reponse: Nt by Nx
%  deconv.metric: metrics for each impulse response: 6 by Nx
%  deconv.peakiness: peakiness score of each response: 1 by Nx
%  deconv.fc: center frequency (Hz): 1 by Nx
%  deconv.dt: fast-time sample spacing (sec)
%  deconv.twtt = 1 by Nx vector containing the exact TWTT of this response
%  deconv.twtt_min = 1 by Nx vector containing the minimum TWTT
%      bin that this impulse response should be used for; deramp only.
%  deconv.twtt_max = 1 by Nx vector with max TWTT bin; deramp only.
%  deconv.param_collate_deconv = parameter structure during call to collate_deconv.
%  deconv.param_analysis = parameter structure during call to analysis.
%  deconv.param_records = parameter structure during call to create records.
% =========================================================================

if param.collate_deconv.stage_one_en
  for img = param.collate_deconv.imgs
    
    if isempty(param.collate_deconv.wf_adcs)
      % If no wf-adc pairs specified, then do them all.
      wf_adcs = 1:size(param.analysis.imgs{img},1);
    else
      wf_adcs = param.collate_deconv.wf_adcs;
    end
    for wf_adc = wf_adcs
      wf = param.analysis.imgs{img}(wf_adcs,1);
      adc = param.analysis.imgs{img}(wf_adcs,2);
      
      %% Stage 1: Load specular analysis file
      % ===================================================================
      fn_dir = fileparts(ct_filename_out(param,param.collate_deconv.in_dir, ''));
      fn = fullfile(fn_dir,sprintf('specular_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
      fprintf('Loading %s img %d wf %d adc %d\n  %s\n', param.day_seg, img, wf, adc, fn);
      spec = load(fn);
      
      %% Preallocate Outputs
      deconv = [];
      deconv.gps_time = [];
      deconv.lat = [];
      deconv.lon = [];
      deconv.elev = [];
      deconv.roll = [];
      deconv.pitch = [];
      deconv.heading = [];
      [~,deconv.frm,deconv.rec] = get_frame_id(param,spec.deconv_gps_time);
      if ~isempty(spec.param_analysis.radar.wfs.ft_wind)
        deconv.ref_windowed = true;
        deconv.ref_window = spec.param_analysis.radar.wfs.ft_wind;
      else
        error('collate_deconv requires a window to have been used.');
      end
      deconv.ref_nonnegative = [];
      deconv.ref_negative = [];
      deconv.ref_mult_factor = [];
      deconv.impulse_response = {};
      deconv.metric = [];
      deconv.peakiness = [];
      deconv.fc = [];
      deconv.dt = spec.dt;
      if strcmpi(radar_type,'deramp')
        deconv.twtt = spec.deconv_twtt;
      end
      param.analysis.cmd{param.collate_deconv.cmd_idx} = cmd;
      deconv.param_collate_deconv = param;
      deconv.param_analysis = spec.param_analysis;
      deconv.param_records = spec.param_records;
      deconv.file_version = '1';
      
      %% Handle the case where no specular targets were found
      % ===================================================================
      if size(spec.deconv_mean,2) == 0
        warning('This segment has no deconvolution waveforms');
        continue;
      end
      
      % ===================================================================
      % ===================================================================
      %% Analyze each waveform, extract impulse response and metrics
      % ===================================================================
      % ===================================================================
      for rline = 1:length(spec.deconv_gps_time)
        
        %% Create impulse response
        
        % h: impulse response
        h = spec.deconv_mean{rline};
        
        % Estimate SNR as a function of range bin
        SNR = lp(abs(h).^2 ./ spec.deconv_std{rline}.^2 * cmd.rlines,1);
        
        % Time gate signal according to cmd.rlines
        Ntg = cmd.rbins(end)-cmd.rbins(1)+1;
        Htg = tukeywin(Ntg);
        h_nonnegative = h(1:1+cmd.rbins(end)) .* Htg(end-cmd.rbins(end):end);
        h_negative = h(end+1+cmd.rbins(1):end) .* Htg(1:-cmd.rbins(1));
        
        %% Check cmd.rlines
        if 0
          h_fig = 1;
          figure(h_fig); clf;
          h_axes = axes('parent',h_fig);
          plot(h_axes(1), lp(spec.deconv_mean{rline}))
          hold(h_axes(1),'on');
          [max_val,max_idx] = max(lp(spec.deconv_sample{rline}));
          plot(h_axes(1), circshift(lp(spec.deconv_sample{rline}) - max_val,[-max_idx 1]))
          plot(h_axes(1), lp(spec.deconv_std{rline},2) - lp(cmd.rlines))
          h_filled = [h_nonnegative; zeros(length(spec.deconv_sample{rline})-length(h_nonnegative)-length(h_negative),1); h_negative];
          plot(h_axes(1), lp(h_filled))
          xlabel(h_axes(1), 'Range bin');
          ylabel(h_axes(1), 'Relative power (dB)');
          title(h_axes(1), 'Impulse response');
          legend(h_axes(1), 'mean','sample','std','h','location','best');
          xlim(h_axes(1), [1 2*cmd.rbins(2)]);
          ylim(h_axes(1), [-50 0]);
          grid(h_axes(1), 'on');
          
          h_fig = 2;
          figure(h_fig); clf;
          h_axes(2) = axes('parent',h_fig);
          plot(h_axes(2), lp(spec.deconv_mean{rline}))
          hold(h_axes(2),'on');
          [max_val,max_idx] = max(lp(spec.deconv_sample{rline}));
          plot(h_axes(2), circshift(lp(spec.deconv_sample{rline}) - max_val,[-max_idx 1]))
          plot(h_axes(2), lp(spec.deconv_std{rline},2) - lp(cmd.rlines))
          h_filled = [h_nonnegative; zeros(length(spec.deconv_sample{rline})-length(h_nonnegative)-length(h_negative),1); h_negative];
          plot(h_axes(2), lp(h_filled))
          xlabel(h_axes(2), 'Range bin');
          ylabel(h_axes(2), 'Relative power (dB)');
          title(h_axes(2), 'Impulse response');
          legend(h_axes(2), 'mean','sample','std','h','location','best');
          xlim(h_axes(2), [Nt+2*cmd.rbins(1) Nt]);
          ylim(h_axes(2), [-50 0]);
          grid(h_axes(2), 'on');
          
          h_fig = 3;
          figure(h_fig); clf;
          h_axes(3) = axes('parent',h_fig);
          plot(h_axes(3),SNR)
          hold(h_axes(3),'on');
          plot(h_axes(3),[1 length(SNR)], 20*[1 1],'k--');
          title(h_axes(3),'SNR');
          xlabel(h_axes(3),'Range bin');
          ylabel(h_axes(3),'SNR (dB)');
          ylim(h_axes(3), [0 50]);
          grid(h_axes(3), 'on');
          
          h_fig = 4;
          figure(h_fig); clf;
          h_axes(4) = axes('parent',h_fig);
          plot(h_axes(4),SNR)
          hold(h_axes(4),'on');
          plot(h_axes(4),[1 length(SNR)], 20*[1 1],'k--');
          title(h_axes(4),'SNR');
          xlabel(h_axes(4),'Range bin');
          ylabel(h_axes(4),'SNR (dB)');
          ylim(h_axes(4), [0 50]);
          grid(h_axes(4), 'on');
          
          linkaxes(h_axes([1 3]),'x');
          linkaxes(h_axes([2 4]),'x');
        end
        
        %% Test deconvolution function
        
        % Create frequency axis
        Nt = length(spec.deconv_sample{rline});
        
        % Get sample rline parameters
        fc = spec.deconv_fc(rline);
        t0 = spec.deconv_t0(rline);
        dt = spec.dt;
        time = t0 + dt*(0:Nt-1).';
        
        % Adjust deconvolution signal to match sample rline
        h_filled = [h_nonnegative; zeros(length(spec.deconv_sample{rline})-length(h_nonnegative)-length(h_negative),1); h_negative];
        deconv_Nt = length(h_filled);
        % Is dt different? Error
        if dt ~= spec.dt
          error('The fast time sample spacing of the data (%g) does not match the deconvolution waveform sampling (%g).',dt,spec.dt);
        end
        % Is fc different? Multiply time domain by exp(1i*2*pi*dfc*deconv_time)
        dfc = fc - spec.deconv_fc(rline);
        if dfc/fc > 1e-6
          deconv_time = t0 + dt*(0:Nt-1).';
          h_filled = h_filled .* exp(1i*2*pi*dfc*deconv_time);
        end
        % Take FFT of deconvolution impulse response
        h_filled = fft(h_filled);
        
        % Create inverse filter relative to window
        df = 1/(Nt*dt);
        freq = spec.deconv_fc(rline) + df * ifftshift(-floor(Nt/2) : floor((Nt-1)/2)).';       
        freq = fftshift(freq);
        Nt_shorten = find(deconv.cmd.f0 <= freq,1);
        Nt_shorten(2) = length(freq) - find(deconv.cmd.f1 >= freq,1,'last');
        Nt_Hwind = Nt - sum(Nt_shorten);
        Hwind = deconv.ref_window(Nt_Hwind);
        Hwind_filled = ifftshift([zeros(Nt_shorten(1),1); Hwind; zeros(Nt_shorten(end),1)]);
        h_filled_inverse = Hwind_filled ./ h_filled;
        
        % Normalize so that reflection is 0 dB (i.e. we assume this
        % specular target is a perfect reflector) at this range, R when
        % voltage is scaled R.
        R = interp1(spec.gps_time,spec.surface,spec.deconv_gps_time(rline)) * c/2;
        % Apply deconvolution with unnormalized filter
        h_deconvolved = ifft(fft(spec.deconv_sample{rline}) .* h_filled_inverse);
        % Oversample to get a good measurement of the peak value
        h_deconvolved = interpft(h_deconvolved,Mt*Nt);
        % Scale so that the peak value is 0dB at param.collate_deconv.R_norm range
        h_mult_factor = param.collate_deconv.R_norm / (R*max(abs(h_deconvolved)));
        h_filled_inverse = h_filled_inverse * h_mult_factor;
        h_deconvolved = h_deconvolved * h_mult_factor;

        % Oversample sample signal by the same amount as the deconvolution
        h_sample = interpft(spec.deconv_sample{rline},Mt*Nt);
        
        % Find maximum values and indices for the deconvolved and
        % undeconvolved signals.
        [deconv_max_val,deconv_max_idx] = max(h_deconvolved);
        [max_val,max_idx] = max(h_sample);
        
        
        %% Check h_deconvolved
        if 0
          
          comp_bins = cmd.rbins(2)+1:2*cmd.rbins(2);
          comp_bins = Nt+2*cmd.rbins(1) : Nt+cmd.rbins;
          dnoise = lp(mean(abs(h_sample(comp_bins)).^2) ./ mean(abs(h_deconvolved(comp_bins)).^2));
          dsignal = lp(max(abs(h_sample).^2) ./ max(abs(h_deconvolved).^2));
          dSNR = dsignal - dnoise;
          
          fprintf('SNR loss: %.1f dB\n', dSNR)
          fprintf('Peak magnitude: %.1f (dB)\n', lp(max_val./deconv_max_val));
          fprintf('Peak angle: %.1f (deg)\n', angle(max_val./deconv_max_val) * 180/pi);
          fprintf('Index offset: %d\n', mod(max_idx - deconv_max_idx + Nt*Mt/2, Nt*Mt)-Nt*Mt/2);
          
          bins_Mt = 0:1/Mt:Nt-1/Mt;
          
          h_fig = 1;
          figure(h_fig); clf;
          h_axes = axes('parent',h_fig);
          [max_val,max_idx] = max(lp(h_sample));
          plot(h_axes(1), bins_Mt, circshift(lp(h_sample) - max_val,[-max_idx 1]))
          hold(h_axes(1),'on');
          plot(h_axes(1), bins_Mt, circshift(lp(h_deconvolved) - max_val + dsignal,[-max_idx 1]))
          xlabel(h_axes(1), 'Range bin');
          ylabel(h_axes(1), 'Relative power (dB)');
          title(h_axes(1), 'Impulse response falling edge');
          legend(h_axes(1), 'sample','deconvolved','location','best');
          xlim(h_axes(1), [0 2*cmd.rbins(2)]);
          ylim(h_axes(1), [-50 0]);
          grid(h_axes(1), 'on');
          
          h_fig = 2;
          figure(h_fig); clf;
          h_axes(2) = axes('parent',h_fig);
          [max_val,max_idx] = max(lp(h_sample));
          plot(h_axes(2), fliplr(bins_Mt+1), circshift(lp(h_sample) - max_val,[-max_idx 1]))
          hold(h_axes(2),'on');
          plot(h_axes(2), fliplr(bins_Mt+1), circshift(lp(h_deconvolved) - max_val + dsignal,[-max_idx 1]))
          xlabel(h_axes(2), 'Range bin');
          ylabel(h_axes(2), 'Relative power (dB)');
          title(h_axes(2), 'Impulse response rising edge');
          legend(h_axes(2), 'sample','deconvolved','location','best');
          xlim(h_axes(2), [0 -2*cmd.rbins(1)]);
          ylim(h_axes(2), [-50 0]);
          grid(h_axes(2), 'on');
          
          h_fig = 3;
          figure(h_fig); clf;
          h_axes(3) = axes('parent',h_fig);
          plot(h_axes(3), lp(h_filled))
          hold(h_axes(3),'on');
          plot(h_axes(3), lp(Hwind_filled,1) - max(lp(Hwind_filled,1)) + max(lp(h_filled)))
          plot(h_axes(3), lp(h_filled_inverse) - max(lp(h_filled_inverse)) + max(lp(h_filled)))
          xlabel(h_axes(3), 'Range bin');
          ylabel(h_axes(3), 'Relative power (dB)');
          title(h_axes(3), 'Transfer function');
          legend(h_axes(3), 'sample','window','inverse','location','best');
          grid(h_axes(3), 'on');
        end
        
        
        %% Create metrics
        h_metric = abs(h_deconvolved).^2; % Convert to power
        [max_val,max_idx] = max(h_metric); % Find maximum value and index
        h_metric = circshift(h_metric,[-max_idx 1]); % Shift peak to first bin
        
        % Compute negated peak power
        peak = -max(lp(max_val,1));
        
        % Compute main lobe width
        main_lobe = 2;
        for bin = 2:cmd.rbins(2)*Mt
          if lp(h_metric(bin),1) < lp(max_val,1)-cmd.ML_threshold
            break;
          end
          main_lobe = main_lobe + 1;
        end
        for bin = 0:-cmd.rbins(1)*Mt
          if lp(h_metric(end-bin),1) < lp(max_val,1)-cmd.ML_threshold
            break;
          end
          main_lobe = main_lobe + 1;
        end
        main_lobe = main_lobe / Mt;
        
        % Compute falling edge peak sidelobe
        peak_sidelobe_falling_edge = -inf;
        for bin = round(cmd.SL_guard_bins*Mt) : cmd.rbins(2)*Mt
          if lp(h_metric(bin),1) > peak_sidelobe_falling_edge
            peak_sidelobe_falling_edge = lp(h_metric(bin),1);
          end
        end
        peak_sidelobe_falling_edge = peak_sidelobe_falling_edge + peak;
        
        % Compute rising edge peak sidelobe
        peak_sidelobe_rising_edge = -inf;
        for bin = round(cmd.SL_guard_bins*Mt)-1:-cmd.rbins(1)*Mt
          if lp(h_metric(end-bin),1) > peak_sidelobe_rising_edge
            peak_sidelobe_rising_edge = lp(h_metric(end-bin),1);
          end
        end
        peak_sidelobe_rising_edge = peak_sidelobe_rising_edge + peak;
        
        % Compute falling edge integrated sidelobe
        bins = round(cmd.SL_guard_bins*Mt) : cmd.rbins(2)*Mt;
        integrated_sidelobe_falling_edge = sum(h_metric(bins))/Mt;
        integrated_sidelobe_falling_edge = lp(integrated_sidelobe_falling_edge,1);
        integrated_sidelobe_falling_edge = integrated_sidelobe_falling_edge + peak;
        
        % Compute falling edge integrated sidelobe
        bins = round(cmd.SL_guard_bins*Mt)-1 : -cmd.rbins(1)*Mt;
        integrated_sidelobe_rising_edge = sum(h_metric(end-bins))/Mt;
        integrated_sidelobe_rising_edge = lp(integrated_sidelobe_rising_edge,1);
        integrated_sidelobe_rising_edge = integrated_sidelobe_rising_edge + peak;
        
        deconv.metric(:,rline) = [peak, main_lobe, peak_sidelobe_falling_edge, peak_sidelobe_rising_edge, ...
          integrated_sidelobe_falling_edge, integrated_sidelobe_rising_edge];
        
        
        %% Prepare deconvolution filter for storage
        [~,match_idx] = min(abs(spec.gps_time - spec.deconv_gps_time(rline)));
        deconv.gps_time(rline) = spec.gps_time(match_idx);
        deconv.lat(rline) = spec.lat(match_idx);
        deconv.lon(rline) = spec.lon(match_idx);
        deconv.elev(rline) = spec.elev(match_idx);
        deconv.roll(rline) = spec.roll(match_idx);
        deconv.pitch(rline) = spec.pitch(match_idx);
        deconv.heading(rline) = spec.heading(match_idx);
        deconv.peakiness(rline) = spec.peakiness(match_idx);
        deconv.fc(rline) = spec.deconv_fc(rline);
        deconv.ref_nonnegative{rline} = h_nonnegative;
        deconv.ref_negative{rline} = h_negative;
        deconv.ref_mult_factor(rline) = h_mult_factor;
        deconv.impulse_response{rline} = h_deconvolved;
        
      end
      
      %% Print metric
      if 0
        % Compare results to metric
        pass = bsxfun(@lt,deconv.metric,cmd.abs_metric(:));
        
        score = nansum(bsxfun(@times, cmd.metric_weights(:), bsxfun(@minus, cmd.abs_metric(:), deconv.metric)));
        
        % Find the highest score in each bin
        twtts = unique(round(deconv.twtt*param.collate_deconv.twtt_penalty*10)/10/param.collate_deconv.twtt_penalty);
        
        max_score_rlines = zeros(size(twtts));
        for twtt_idx = 1:length(twtts)
          twtt = twtts(twtt_idx);
          [~,max_score_rlines(twtt_idx)] = max(score(abs(deconv.twtt-twtt) < 1/param.collate_deconv.twtt_penalty/10*2));
        end
        
        % Plot metrics
        
        h_fig = 1;
        figure(h_fig); clf;
        h_axes = axes('parent',h_fig);
        plot(h_axes(1), find(all(pass)), deconv.metric(1,all(pass)), '.-')
        hold(h_axes(1),'on');
        plot(h_axes(1), find(all(pass)), deconv.metric(2,all(pass)), '.-')
        plot(h_axes(1), find(all(pass)), deconv.metric(3,all(pass)), '.-')
        plot(h_axes(1), find(all(pass)), deconv.metric(4,all(pass)), '.-')
        plot(h_axes(1), find(all(pass)), deconv.metric(5,all(pass)), '.-')
        plot(h_axes(1), find(all(pass)), deconv.metric(6,all(pass)), '.-')
        xlabel(h_axes(1), 'Deconv waveform index');
        ylabel(h_axes(1), 'Metric (lower is better)');
        title(h_axes(1), sprintf('%s: pass', param.day_seg), 'interpreter', 'none');
        legend(h_axes(1), 'P','ML','FSL','RSL','IFSL','IRSL','location','best');
        grid(h_axes(1), 'on');
        
        h_fig = 2;
        figure(h_fig); clf;
        h_axes(2) = axes('parent',h_fig);
        h_plot = plot(h_axes(2), deconv.metric(1,:));
        hold(h_axes(2),'on');
        h_plot(2) = plot(h_axes(2), deconv.metric(2,:));
        h_plot(3) = plot(h_axes(2), deconv.metric(3,:));
        h_plot(4) = plot(h_axes(2), deconv.metric(4,:));
        h_plot(5) = plot(h_axes(2), deconv.metric(5,:));
        h_plot(6) = plot(h_axes(2), deconv.metric(6,:));
        plot(h_axes(2), find(all(pass)), deconv.metric(1,all(pass)), '.','Color',get(h_plot(1),'Color'))
        plot(h_axes(2), find(all(pass)), deconv.metric(2,all(pass)), '.','Color',get(h_plot(2),'Color'))
        plot(h_axes(2), find(all(pass)), deconv.metric(3,all(pass)), '.','Color',get(h_plot(3),'Color'))
        plot(h_axes(2), find(all(pass)), deconv.metric(4,all(pass)), '.','Color',get(h_plot(4),'Color'))
        plot(h_axes(2), find(all(pass)), deconv.metric(5,all(pass)), '.','Color',get(h_plot(5),'Color'))
        plot(h_axes(2), find(all(pass)), deconv.metric(6,all(pass)), '.','Color',get(h_plot(6),'Color'))
        xlabel(h_axes(2), 'Deconv waveform index');
        ylabel(h_axes(2), 'Metric (lower is better)');
        title(h_axes(2), param.day_seg, 'interpreter', 'none');
        legend(h_axes(2), 'P','ML','FSL','RSL','IFSL','IRSL','location','best');
        grid(h_axes(2), 'on');
        
        linkaxes(h_axes);
        
        twtt_bin = deconv.twtt * param.collate_deconv.twtt_penalty;
        
        h_fig = 3;
        figure(h_fig); clf;
        h_axes(3) = axes('parent',h_fig);
        plot(h_axes(3), twtt_bin)
        hold(h_axes(3),'on');
        %plot(h_axes(3), find(all(pass)), twtt_bin(all(pass)), '.','markersize',14)
        scatter(h_axes(3), find(all(pass)),twtt_bin(all(pass)), 200,score(all(pass)).','.')
        xlabel(h_axes(3), 'Deconv waveform index');
        ylabel(h_axes(3), 'TWTT bin');
        title(h_axes(3), sprintf('%s', param.day_seg), 'interpreter', 'none');
        legend(h_axes(3), 'TWTT','Passed','location','best');
        grid(h_axes(3), 'on');
        h_colorbar = colorbar;
        set(get(h_colorbar,'YLabel'),'String','Score');
        
        h_fig = 4;
        figure(h_fig); clf;
        h_axes(4) = axes('parent',h_fig);
        plot(h_axes(4), spec.gps_time, spec.surface * param.collate_deconv.twtt_penalty);
        hold(h_axes(4),'on');
        scatter(h_axes(4), deconv.gps_time(all(pass)),twtt_bin(all(pass)), 200,score(all(pass)).','.')
        scatter(h_axes(4), deconv.gps_time(all(pass)),twtt_bin(all(pass)), 200,score(all(pass)).','x')
        h_colorbar = colorbar;
        set(get(h_colorbar,'YLabel'),'String','Score');
        xlabel(h_axes(4), 'GPS time (sec)');
        ylabel(h_axes(4), 'TWTT bin');
        title(h_axes(4), sprintf('%s', param.day_seg), 'interpreter', 'none');
        grid(h_axes(4), 'on');
        xlim(h_axes(4), spec.gps_time([1 end]))
        ylim(h_axes(4), [min(5,min(spec.surface * param.collate_deconv.twtt_penalty)) max(15,max(spec.surface * param.collate_deconv.twtt_penalty))]);
        
        % Print table
        fprintf('Metric Threshold (Must be below this to pass)\n');
        fprintf('Peak\tML\tPSL FE\tPSL RS\tISL FE\tISL RE\n');
        fprintf('%.1f\t', cmd.abs_metric); fprintf('\n');
        fprintf('Metric ( '); fprintf(2,'Red Failed '); fprintf(')\n');
        fprintf('INDEX\tFRM\tREC\tPeak\tML\tPSL FE\tPSL RS\tISL FE\tISL RE\tPASS\tSCORE\tTWTT\n');
        for rline = 1:length(deconv.gps_time)
          fprintf('%d\t',rline);
          fprintf('%d\t',deconv.frm(rline));
          fprintf('%d\t',deconv.rec(rline));
          for metric = 1:6
            if pass(metric,rline)
              %fprintf('<strong>%.1f</strong>\t', deconv.metric(metric,rline));
              fprintf('%.1f\t', deconv.metric(metric,rline));
            else
              fprintf(2,'%.1f\t', deconv.metric(metric,rline));
            end
          end
          if all(pass(:,rline))
            fprintf('*\t');
          else
            fprintf('\t');
          end
          max_score_idx = find(rline == max_score_rlines);
          if ~isempty(max_score_idx)
            fprintf('<strong>%.1f\t', score(rline));
            fprintf('%.3g ', twtts(max_score_idx)); fprintf('</strong>');
          else
            fprintf('%.1f\t', score(rline));
          end
          fprintf('\n');
        end
      end
      
      %% Store result
      fn_dir = fileparts(ct_filename_out(param,param.collate_deconv.out_dir, ''));
      if ~exist(fn_dir,'dir')
        mkdir(fn_dir);
      end
      out_fn = fullfile(fn_dir,sprintf('deconv_lib_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
      fprintf('Saving %s img %d wf %d adc %d\n  %s\n', param.day_seg, img, wf, adc, out_fn);
      file_locked = false;
      if exist(out_fn,'file')
        tmp = load(out_fn,'file_version');
        if isfield(tmp,'file_version')
          file_locked = ~isempty(tmp.file_version=='L');
        end
      end
      if file_locked
        error('  File is locked.');
      end
      save(out_fn,'-v7.3','-struct','deconv');
      
    end
  end
end

if param.collate_deconv.stage_two_en
  for img = param.collate_deconv.imgs
    
    if isempty(param.collate_deconv.wf_adcs)
      % If no wf-adc pairs specified, then do them all.
      wf_adcs = 1:size(param.analysis.imgs{img},1);
    else
      wf_adcs = param.collate_deconv.wf_adcs;
    end
    for wf_adc = wf_adcs
      wf = param.analysis.imgs{img}(wf_adcs,1);
      adc = param.analysis.imgs{img}(wf_adcs,2);
      
      %% Stage 2: Determine the best waveform for each record
      % ===================================================================
      % 1. Load all segments that are specified (default is just this segment)
      cmd.day_segs = {param.day_seg};
      for day_seg_idx = 1:length(cmd.day_segs)
        day_seg = cmd.day_segs{day_seg_idx};
        fn = fullfile(fn_dir,sprintf('deconv_lib_%s_wf_%d_adc_%d.mat', day_seg, wf, adc));
        fprintf('Loading %s img %d wf %d adc %d\n  %s\n', day_seg, img, wf, adc, fn);
        deconv = load(fn);
        
        deconv.day_seg = repmat({day_seg},[length(deconv.gps_time) 1]);
      end
      if isempty(deconv.gps_time)
        error('There are no deconvolution waveforms in the files loaded. Specify other cmd.day_seg to load or remake current day_seg files with lower metric thresholds.');
      end
      
      % 2. Load surface using opsLoadLayers to determine which waveforms
      %    are needed
      layer_params.name = 'surface';
      layer_params.source = 'records';
      layer = opsLoadLayers(param,layer_params);
      
      % 3. Decimate layer
      decim_idxs = get_equal_alongtrack_spacing_idxs(layer.gps_time,10);
      layer.gps_time = layer.gps_time(decim_idxs);
      layer.twtt = layer.twtt(decim_idxs);
      layer.lat = layer.lat(decim_idxs);
      layer.lon = layer.lon(decim_idxs);
      layer.elev = layer.elev(decim_idxs);
      
      % 4. Compare results to metric
      pass = bsxfun(@lt,deconv.metric,cmd.abs_metric(:));
      score = nansum(bsxfun(@times, cmd.metric_weights(:), bsxfun(@minus, cmd.abs_metric(:), deconv.metric)));
      
      % 5. Find best score at each point along the flight track
      min_score = min(score);
      score = score-min_score;
      for rline = 1:length(layer.twtt)
        % Score with twtt penalty and time constant penalty term
        d_twtt = layer.twtt(rline) - deconv.twtt;
        d_gps_time = layer.gps_time(rline) - deconv.gps_time;
        adjusted_score = min_score + score .* exp(-abs(param.collate_deconv.twtt_penalty*d_twtt).^2) ...
          .* exp(-abs(param.collate_deconv.gps_time_penalty*d_gps_time).^2);
        
        [max_score(rline),max_idx(rline)] = max(adjusted_score);
      end
      if any(max_score < param.collate_deconv.min_score)
        warning('Score is too low for %d blocks of range lines.', sum(max_score < param.collate_deconv.min_score));
      end
      [max_idxs,~,max_idxs_mapping] = unique(max_idx);
      
      % 6. Create the final output structure
      final = [];
      final.gps_time = deconv.gps_time(max_idxs);
      final.lat = deconv.lat(max_idxs);
      final.lon = deconv.lon(max_idxs);
      final.elev = deconv.elev(max_idxs);
      final.roll = deconv.roll(max_idxs);
      final.pitch = deconv.pitch(max_idxs);
      final.heading = deconv.heading(max_idxs);
      final.frm = deconv.frm(max_idxs);
      final.rec = deconv.rec(max_idxs);
      final.ref_windowed = deconv.ref_windowed;
      final.ref_window = deconv.ref_window;
      final.ref_nonnegative = deconv.ref_nonnegative(max_idxs);
      final.ref_negative = deconv.ref_negative(max_idxs);
      final.ref_mult_factor = deconv.ref_mult_factor(max_idxs);
      final.impulse_response = deconv.impulse_response(max_idxs);
      final.metric = deconv.metric(:,max_idxs);
      final.peakiness = deconv.peakiness(max_idxs);
      final.fc = deconv.fc(max_idxs);
      final.dt = deconv.dt;
      final.twtt = deconv.twtt(max_idxs);
      final.param_collate_deconv_final = param;
      final.param_collate_deconv = deconv.param_collate_deconv;
      final.param_analysis = deconv.param_analysis;
      final.param_records = deconv.param_records;
      final.map_day_seg = deconv.day_seg(max_idxs);
      final.map_gps_time = layer.gps_time;
      final.map_idxs = max_idxs_mapping(:).';
      final.max_score = max_score;
      final.file_version = '1';
      
      % 7. Store final output file
      fn_dir = fileparts(ct_filename_out(param,param.collate_deconv.out_dir, ''));
      if ~exist(fn_dir,'dir')
        mkdir(fn_dir);
      end
      out_fn = fullfile(fn_dir,sprintf('deconv_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
      fprintf('Saving %s img %d wf %d adc %d\n  %s\n', param.day_seg, img, wf, adc, out_fn);
      ct_file_lock_check(out_fn);
      save(out_fn,'-v7.3','-struct','final');
      
    end
  end
end
