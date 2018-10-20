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
  param.collate_deconv.out_dir = param.collate_deconv.in_dir;
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

if ~isfield(cmd,'day_segs') || isempty(cmd.day_segs)
  % Default is to use this day_seg only to find deconvolution waveforms
  cmd.day_segs = {param.day_seg};
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

if param.collate_deconv.stage_one_en && (~isfield(cmd,'rbins') || isempty(cmd.rbins))
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
      if isempty(spec.deconv_gps_time)
        warning('No specular waveforms found.');
        continue;
      end
      
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
      if param.ct_file_lock
        deconv.file_version = '1L';
      else
        deconv.file_version = '1';
      end
      
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
      
      %% HACK: To deal with broken data_pulse_compress
      hack_records = load(ct_filename_support(param,'','records'));
%       hack_trajectory_param = struct('gps_source',hack_records.gps_source, ...
%   'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
%   'tx_weights', [], 'lever_arm_fh', param.radar.lever_arm_fh);
% hack_records = trajectory_with_leverarm(hack_records,hack_trajectory_param);

      % Break records in segment into blocks
      hack_breaks = 1:spec.param_analysis.analysis.block_size:length(hack_records.gps_time);
      
      % If the last block is less than half the desired block size, then combine
      % with earlier block if possible
      if length(hack_records.gps_time)-hack_breaks(end) < spec.param_analysis.analysis.block_size/2 ...
          && length(hack_breaks) > 1
        hack_breaks = hack_breaks(1:end-1);
      end
      hack_breaks(end+1) = length(hack_records.gps_time)+1;
      hack_elev = interp1(hack_records.gps_time,hack_records.elev,spec.deconv_gps_time,'nearest','extrap');
      for rline = 1:length(spec.deconv_gps_time)
        hack_block_idx = find(deconv.rec(rline) >= hack_breaks,1,'last');
        hack_elev_block = hack_records.elev(hack_breaks(hack_block_idx):hack_breaks(hack_block_idx+1)-1);
        hack_twtt(rline) = spec.deconv_twtt(rline) + (hack_elev(rline) - hack_elev_block(1))/(c/2);
      end
      figure(1); clf;
      plot(spec.deconv_twtt);
      hold on
      plot(hack_twtt)
      hack_end = spec.deconv_t0+spec.param_analysis.radar.fs/2*spec.param_analysis.radar.wfs.Tpd/(spec.param_analysis.radar.wfs.f1-spec.param_analysis.radar.wfs.f0);
      plot(hack_end)
      hack_bad_guard = 2e-9;
      hack_bad_mask = abs(mod(hack_twtt+hack_bad_guard-min(spec.deconv_t0),spec.param_analysis.radar.fs/2*spec.param_analysis.radar.wfs.Tpd/(spec.param_analysis.radar.wfs.f1-spec.param_analysis.radar.wfs.f0))) < 2*hack_bad_guard;
      plot(find(~hack_bad_mask),hack_twtt(~hack_bad_mask),'.');
      %% END HACK

      for rline = 1:length(spec.deconv_gps_time)
        
        %% Create impulse response
        
        % HACK: To deal with broken data_pulse_compress
        if hack_bad_mask(rline)
          fprintf('HACK: skipping %d: %g\n', rline, hack_twtt(rline));
          deconv.metric(:,rline) = NaN(6,1);
          [~,match_idx] = min(abs(spec.gps_time - spec.deconv_gps_time(rline)));
          deconv.gps_time(rline) = spec.gps_time(match_idx);
          continue
        end
        % END HACK
        
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
          
          keyboard
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
        Nt_shorten = find(cmd.f0 <= freq,1);
        Nt_shorten(2) = length(freq) - find(cmd.f1 >= freq,1,'last');
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

        % Oversample sample signal by the same amount as the deconvolved
        % signal
        h_sample = interpft(spec.deconv_sample{rline},Mt*Nt);
        
        % Scale so that the peak value is 1/R
        h_mult_factor = max(abs(h_sample)) / max(abs(h_deconvolved));
        h_filled_inverse = h_filled_inverse * h_mult_factor;
        h_deconvolved = h_deconvolved * h_mult_factor;
        
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
          xlabel(h_axes(3), 'Frequency bin');
          ylabel(h_axes(3), 'Relative power (dB)');
          title(h_axes(3), 'Transfer function');
          legend(h_axes(3), 'sample','window','inverse','location','best');
          grid(h_axes(3), 'on');
          
          keyboard
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
      if 1
        % Compare results to metric
        pass = bsxfun(@lt,deconv.metric,cmd.abs_metric(:));
        
        score = nansum(bsxfun(@times, cmd.metric_weights(:), bsxfun(@minus, cmd.abs_metric(:), deconv.metric)));
        score(:,any(isnan(deconv.metric))) = -inf;
        
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
        
        try; linkaxes(h_axes); end;
        
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
        fprintf('Peak\tML\tPSL FE\tPSL RE\tISL FE\tISL RE\n');
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
      ct_file_lock_check(out_fn,2);
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

      % Load deconv
      fn_dir = fileparts(ct_filename_out(param,param.collate_deconv.out_dir, ''));
      fn = fullfile(fn_dir,sprintf('deconv_lib_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
      fprintf('Loading %s img %d wf %d adc %d\n  %s\n', param.day_seg, img, wf, adc, fn);
      if exist(fn)
        deconv = load(fn,'param_collate_deconv','param_analysis','param_records');
        new_file = false;
      else
        fprintf('  Does not exist.\n');
        new_file = true;
      end
      
      % 1. Load all segments that are specified (default is just this segment)
      for day_seg_idx = 1:length(cmd.day_segs)
        day_seg = cmd.day_segs{day_seg_idx};
        fn_dir = fileparts(ct_filename_out(param,param.collate_deconv.out_dir, ''));
        fn = fullfile(fn_dir,sprintf('deconv_lib_%s_wf_%d_adc_%d.mat', day_seg, wf, adc));
        fprintf('Loading lib %s img %d wf %d adc %d\n  %s\n', day_seg, img, wf, adc, fn);
        if day_seg_idx == 1
          deconv_lib = load(fn);
          deconv_lib.day_seg = repmat({day_seg},[length(deconv_lib.gps_time) 1]);
        else
          tmp = load(fn);
          deconv_lib.gps_time(end+(1:length(tmp.gps_time))) = tmp.gps_time;
          deconv_lib.lat(end+(1:length(tmp.lat))) = tmp.lat;
          deconv_lib.lon(end+(1:length(tmp.lon))) = tmp.lon;
          deconv_lib.elev(end+(1:length(tmp.elev))) = tmp.elev;
          deconv_lib.roll(end+(1:length(tmp.roll))) = tmp.roll;
          deconv_lib.pitch(end+(1:length(tmp.pitch))) = tmp.pitch;
          deconv_lib.heading(end+(1:length(tmp.heading))) = tmp.heading;
          deconv_lib.frm(end+(1:length(tmp.frm))) = tmp.frm;
          deconv_lib.rec(end+(1:length(tmp.rec))) = tmp.rec;
          deconv_lib.ref_nonnegative(end+(1:length(tmp.ref_nonnegative))) = tmp.ref_nonnegative;
          deconv_lib.ref_negative(end+(1:length(tmp.ref_negative))) = tmp.ref_negative;
          deconv_lib.ref_mult_factor(end+(1:length(tmp.ref_mult_factor))) = tmp.ref_mult_factor;
          deconv_lib.impulse_response(end+(1:length(tmp.impulse_response))) = tmp.impulse_response;
          deconv_lib.metric(:,end+(1:length(tmp.metric))) = tmp.metric;
          deconv_lib.peakiness(end+(1:length(tmp.peakiness))) = tmp.peakiness;
          deconv_lib.fc(end+(1:length(tmp.fc))) = tmp.fc;
          deconv_lib.twtt(end+(1:length(tmp.twtt))) = tmp.twtt;
          deconv_lib.day_seg(end+(1:length(tmp.gps_time))) = repmat({day_seg},[length(tmp.gps_time) 1]);
        end
        if new_file
          deconv.param_collate_deconv = deconv_lib.param_collate_deconv;
          deconv.param_analysis = deconv_lib.param_analysis;
          deconv.param_records = deconv_lib.param_records;
        end
        
      end
      if isempty(deconv_lib.gps_time)
        error('There are no deconvolution waveforms in the files loaded. Specify other cmd.day_seg to load or remake current day_seg files with lower metric thresholds.');
      end
      
      % 2. Load surface using opsLoadLayers to determine which waveforms
      %    are needed
      layer_params = [];
      layer_params.name = 'surface';
      layer_params.source = 'layerData';
      layer = opsLoadLayers(param,layer_params);
      
      % 3. Decimate layer
      decim_idxs = get_equal_alongtrack_spacing_idxs(layer.gps_time,10);
      layer.gps_time = layer.gps_time(decim_idxs);
      layer.twtt = layer.twtt(decim_idxs);
      layer.lat = layer.lat(decim_idxs);
      layer.lon = layer.lon(decim_idxs);
      layer.elev = layer.elev(decim_idxs);
      
      % 4. Compare results to metric
      pass = bsxfun(@lt,deconv_lib.metric,cmd.abs_metric(:));
      score = nansum(bsxfun(@times, cmd.metric_weights(:), bsxfun(@minus, cmd.abs_metric(:), deconv_lib.metric)));
      score(:,any(isnan(deconv_lib.metric))) = nan;
      
      % 5. Find best score at each point along the flight track
      min_score = nanmin(score);
      score = score-min_score;
      clear max_score unadjusted_score max_idx;
      for rline = 1:length(layer.twtt)
        % Score with twtt penalty and time constant penalty term
        d_twtt = layer.twtt(rline) - deconv_lib.twtt;
        d_gps_time = layer.gps_time(rline) - deconv_lib.gps_time;
        %adjusted_score = min_score + score .* exp(-abs(param.collate_deconv.twtt_penalty*d_twtt).^2) ...
        %  .* exp(-abs(param.collate_deconv.gps_time_penalty*d_gps_time).^2);
        adjusted_score = min_score + score - (100-100*exp(-abs(param.collate_deconv.twtt_penalty*d_twtt).^2)) ...
          - (50-50*exp(-abs(param.collate_deconv.gps_time_penalty*d_gps_time).^2));
        
        [max_score(rline),max_idx(rline)] = max(adjusted_score);
        unadjusted_score(rline) = min_score + score(max_idx(rline));
      end
      if any(max_score < param.collate_deconv.min_score)
        warning('Score is too low for %d of %d blocks of range lines.', sum(max_score < param.collate_deconv.min_score), length(max_score));
      end
      [max_idxs,~,max_idxs_mapping] = unique(max_idx);
      
      [~,sort_idxs] = sort(deconv_lib.twtt(max_idxs));
      unsort_idxs(sort_idxs) = 1:length(sort_idxs);
      max_idxs = max_idxs(sort_idxs);
      max_idxs_mapping = unsort_idxs(max_idxs_mapping);
      
      % 6. Create the final output structure
      final = [];
      final.gps_time = deconv_lib.gps_time(max_idxs);
      final.lat = deconv_lib.lat(max_idxs);
      final.lon = deconv_lib.lon(max_idxs);
      final.elev = deconv_lib.elev(max_idxs);
      final.roll = deconv_lib.roll(max_idxs);
      final.pitch = deconv_lib.pitch(max_idxs);
      final.heading = deconv_lib.heading(max_idxs);
      final.frm = deconv_lib.frm(max_idxs);
      final.rec = deconv_lib.rec(max_idxs);
      final.ref_windowed = deconv_lib.ref_windowed;
      final.ref_window = deconv_lib.ref_window;
      final.ref_nonnegative = deconv_lib.ref_nonnegative(max_idxs);
      final.ref_negative = deconv_lib.ref_negative(max_idxs);
      final.ref_mult_factor = deconv_lib.ref_mult_factor(max_idxs);
      final.impulse_response = deconv_lib.impulse_response(max_idxs);
      final.metric = deconv_lib.metric(:,max_idxs);
      final.peakiness = deconv_lib.peakiness(max_idxs);
      final.fc = deconv_lib.fc(max_idxs);
      final.dt = deconv_lib.dt;
      final.twtt = deconv_lib.twtt(max_idxs);
      final.param_collate_deconv_final = param;
      final.param_collate_deconv = deconv.param_collate_deconv;
      final.param_analysis = deconv.param_analysis;
      final.param_records = deconv.param_records;
      final.map_day_seg = deconv_lib.day_seg(max_idxs);
      final.map_gps_time = layer.gps_time;
      final.map_twtt = layer.twtt;
      final.map_idxs = max_idxs_mapping(:).';
      final.max_score = max_score;
      final.unadjusted_score = unadjusted_score;
      if param.ct_file_lock
        final.file_version = '1L';
      else
        final.file_version = '1';
      end

      % 7. Plot final deconv waveforms
      
      % TWTT Figure
      % ===================================================================
      h_fig = 1;
      figure(h_fig); clf;
      h_axes = axes('parent',h_fig);
      
      legend_str = {};
      h_plot = [];
      for idx = 1:length(final.gps_time)
        h_plot(idx+1) = plot(h_axes(1), find(final.map_idxs==idx), final.twtt(final.map_idxs(final.map_idxs==idx)),'.');
        hold(h_axes(1),'on');
        legend_str{idx+1} = sprintf('%d',idx);
      end
      
      h_plot(1) = plot(h_axes(1), final.map_twtt, 'k', 'LineWidth',2);
      legend_str{1} = 'TWTT';
      
      xlabel(h_axes(1), 'Block');
      ylabel(h_axes(1), 'Two way travel time (\mus)');
      title(h_axes(1), param.day_seg, 'interpreter', 'none');
      legend(h_axes(1), h_plot, legend_str,'location','best');
      grid(h_axes(1), 'on');
      
      fig_fn = ct_filename_ct_tmp(param,'','collate_deconv',sprintf('%s_twtt',param.collate_deconv.out_dir));
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      saveas(h_fig,[fig_fn '.fig']);
      saveas(h_fig,[fig_fn '.jpg']);
      
      % Score Figure
      % ===================================================================
      h_fig = 2;
      figure(h_fig); clf;
      h_axes(2) = axes('parent',h_fig);
      
      legend_str = {};
      h_plot = [];
      for idx = 1:length(final.gps_time)
        h_plot(idx) = plot(h_axes(2), find(final.map_idxs==idx), final.max_score(find(final.map_idxs==idx)),'.');
        hold(h_axes(2),'on');
        legend_str{idx} = sprintf('%d',idx);
      end
      for idx = 1:length(final.gps_time)
        h_new_plot = plot(h_axes(2), find(final.map_idxs==idx), final.unadjusted_score(find(final.map_idxs==idx)),'.');
        set(h_new_plot, 'Color', get(h_plot(idx),'Color'))
      end
      
      xlabel(h_axes(2), 'Block');
      ylabel(h_axes(2), 'Score');
      title(h_axes(2), param.day_seg, 'interpreter', 'none');
      legend(h_axes(2), h_plot, legend_str,'location','best');
      grid(h_axes(2), 'on');
      
      fig_fn = ct_filename_ct_tmp(param,'','collate_deconv',sprintf('%s_score',param.collate_deconv.out_dir));
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      saveas(h_fig,[fig_fn '.fig']);
      saveas(h_fig,[fig_fn '.jpg']);
      
      % Transfer Function Figure
      % ===================================================================
      h_fig = 3;
      figure(h_fig); clf;
      pos = get(h_fig,'Position');
      set(h_fig,'Position',[pos(1:2) 1000 600]);
      h_axes(3) = subplot(5,1,1:2,'parent',h_fig);
      h_axes(4) = subplot(5,1,3:5,'parent',h_fig);
      
      legend_str = {};
      for idx = 1:length(final.gps_time)
        % Get the reference function
        h_nonnegative = final.ref_nonnegative{idx};
        h_negative = final.ref_negative{idx};
        h_mult_factor = final.ref_mult_factor(idx);
        
        % Adjust deconvolution signal to match sample rline
        h_filled = [h_nonnegative; h_negative];
        
        % Take FFT of deconvolution impulse response
        h_filled = fft(h_filled);
        
        Nt = numel(h_filled);
        df = 1/(Nt*final.dt);
        freq = final.fc(idx) + df * ifftshift(-floor(Nt/2) : floor((Nt-1)/2)).';
        freq = fftshift(freq);
        fc_idx = find(freq==final.fc(idx));
        
        h_filled_lp = fftshift(lp(h_filled));
        h_filled_phase = fftshift(angle(h_filled));
        h_filled_phase = unwrap(h_filled_phase)*180/pi;
        h_filled_phase = h_filled_phase - h_filled_phase(fc_idx);
        
        if max(freq) > 2e9
          freq_scale = 1e9;
        else
          freq_scale = 1e6;
        end
        plot(h_axes(3), freq/freq_scale, h_filled_lp);
        hold(h_axes(3),'on');
        plot(h_axes(4), freq/freq_scale, h_filled_phase);
        hold(h_axes(4),'on');
        legend_str{idx} = sprintf('%d %s_%03d %4.0f %4.1fus',idx, ...
          final.map_day_seg{idx},final.frm(idx),round(lp(final.ref_nonnegative{idx}(1),2)), ...
          round(final.twtt(idx)*1e7)/10);
      end
      
      if freq_scale == 1e9
        xlabel(h_axes(4), 'Frequency (GHz)');
      else
        xlabel(h_axes(4), 'Frequency (MHz)');
      end
      ylabel(h_axes(3), 'Relative power (dB)');
      ylabel(h_axes(4), 'Relative angle (deg)');
      title(h_axes(3), sprintf('%s (Legend idx:frm:peak:twtt)', param.day_seg), 'interpreter', 'none');
      grid(h_axes(3), 'on');
      grid(h_axes(4), 'on');
      h_legend = legend(h_axes(3), legend_str, 'location', 'northeastoutside', 'interpreter','none');
      drawnow;
      pos3 = get(h_axes(3),'Position');
      pos4 = get(h_axes(4),'Position');
      set(h_axes(4),'Position',[pos4(1:2) pos3(3) pos4(4)]);
      
      fig_fn = ct_filename_ct_tmp(param,'','collate_deconv',sprintf('%s_transfer_func',param.collate_deconv.out_dir));
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      saveas(h_fig,[fig_fn '.fig']);
      saveas(h_fig,[fig_fn '.jpg']);
      
      
      % 8. Store final output file
      fn_dir = fileparts(ct_filename_out(param,param.collate_deconv.out_dir, ''));
      if ~exist(fn_dir,'dir')
        mkdir(fn_dir);
      end
      out_fn = fullfile(fn_dir,sprintf('deconv_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
      fprintf('Saving %s img %d wf %d adc %d\n  %s\n', param.day_seg, img, wf, adc, out_fn);
      ct_file_lock_check(out_fn,2);
      save(out_fn,'-v7.3','-struct','final');
      
    end
  end
end
