% function collate_coh_noise(param,param_override)
% collate_coh_noise(param,param_override)
%
% Collects analysis.m results from coherent noise tracking (coh_noise
% command) and creates files for removing the coherent noise during data
% loading.
% Loads all the coh_noise_* files and creates coh_noise_simp_* files that
% are pre-filtered for speed and saves as netcdf so that subsets of the
% files can be loaded efficiently.
%
% Example:
%  See run_collate_coh_noise for how to run.
%
% Authors: John Paden

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =====================================================================

if ~isfield(param.collate_coh_noise,'cmd_idx') || isempty(param.collate_coh_noise.cmd_idx)
  param.collate_coh_noise.cmd_idx = 1;
end

if ~isfield(param.collate_coh_noise,'debug_plots') || isempty(param.collate_coh_noise.debug_plots)
  param.collate_coh_noise.debug_plots = {};
end
enable_visible_plot = any(strcmp('visible',param.collate_coh_noise.debug_plots));
enable_threshold_plot = any(strcmp('threshold_plot',param.collate_coh_noise.debug_plots));
enable_cn_plot = any(strcmp('cn_plot',param.collate_coh_noise.debug_plots));
if ~isempty(param.collate_coh_noise.debug_plots)
  h_fig = get_figures(5,enable_visible_plot);
end

if ~isfield(param.collate_coh_noise,'dft_corr_time') || isempty(param.collate_coh_noise.dft_corr_time)
  param.collate_coh_noise.dft_corr_time = inf;
end

if ~isfield(param.collate_coh_noise,'firdec_fs') || isempty(param.collate_coh_noise.firdec_fs)
  param.collate_coh_noise.firdec_fs = 0;
end

if ~isfield(param.collate_coh_noise,'firdec_fcutoff') || isempty(param.collate_coh_noise.firdec_fcutoff)
  param.collate_coh_noise.firdec_fcutoff = 0;
end

if ~isfield(param.analysis,'imgs') || isempty(param.analysis.imgs)
  param.analysis.imgs = {[1 1]};
end
if ~isfield(param.collate_coh_noise,'imgs') || isempty(param.collate_coh_noise.imgs)
  param.collate_coh_noise.imgs = 1:length(param.analysis.imgs);
end

if ~isfield(param.collate_coh_noise,'in_path') || isempty(param.collate_coh_noise.in_path)
  param.collate_coh_noise.in_path = 'analysis';
end

if ~isfield(param.collate_coh_noise,'method') || isempty(param.collate_coh_noise.method)
  param.collate_coh_noise.method = 'dft';
end

if ~isfield(param.collate_coh_noise,'out_path') || isempty(param.collate_coh_noise.out_path)
  param.collate_coh_noise.out_path = 'analysis';
end

if ~isfield(param.collate_coh_noise,'threshold_en') || isempty(param.collate_coh_noise.threshold_en)
  param.collate_coh_noise.threshold_en = false;
end
% Enable threshold if threshold_plot is enabled
if enable_threshold_plot
  param.collate_coh_noise.threshold_en = true;
end
enable_threshold = param.collate_coh_noise.threshold_en;

if ~isfield(param.collate_coh_noise,'threshold_eval') || isempty(param.collate_coh_noise.threshold_eval)
  % Default is no evaluation
  % Example of how to use evaluation:
  %   params(param_idx).collate_coh_noise.threshold_eval{wf} = 'tmp=threshold(max(1,round(1.1*Tpd_bin)+60):end); tmp(tmp>-110)=-110; threshold(max(1,round(1.1*Tpd_bin)+60):end)=tmp; threshold=threshold+10;';
  param.collate_coh_noise.threshold_eval = {};
end

if ~isfield(param.collate_coh_noise,'threshold_ylims') || isempty(param.collate_coh_noise.threshold_ylims)
  param.collate_coh_noise.threshold_ylims = [];
end

if ~isfield(param.collate_coh_noise,'wf_adcs') || isempty(param.collate_coh_noise.wf_adcs)
  param.collate_coh_noise.wf_adcs = [];
end
if ~isempty(param.collate_coh_noise.wf_adcs) && ~iscell(param.collate_coh_noise,'wf_adcs')
  param.collate_coh_noise.wf_adcs = {param.collate_coh_noise.wf_adcs};
end

for img = param.collate_coh_noise.imgs
  
  if isempty(param.collate_coh_noise.wf_adcs)
    wf_adcs = 1:size(param.analysis.imgs{img},1);
  else
    wf_adcs = param.collate_coh_noise.wf_adcs{img};
  end
  for wf_adc = wf_adcs
    wf = param.analysis.imgs{img}(wf_adc,1);
    adc = param.analysis.imgs{img}(wf_adc,2);
    
    %% Load the coherent noise file
    % =====================================================================
    fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.in_path));
    fn = fullfile(fn_dir,sprintf('coh_noise_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
    fprintf('Loading %s (%s)\n', fn, datestr(now));
    noise = load(fn);
    
    cmd = noise.param_analysis.analysis.cmd{param.collate_coh_noise.cmd_idx};
    if ~isfield(cmd,'min_samples') || isempty(cmd.min_samples)
      cmd.min_samples = 0.5*cmd.block_ave;
    end
    
    % Each processed block has a constant start and stop bin, but between
    % blocks the start/stop range bin may change.
    start_bins = round(noise.t0/noise.dt);
    stop_bins = round(noise.t0/noise.dt) + noise.Nt - 1;
    % Determine the first start bin and the last stop bin
    start_bin = min(start_bins);
    stop_bin = max(stop_bins);
    Nt = stop_bin-start_bin+1;
      
    % Get the size of each processing block
    block_size = cellfun(@(t) size(t,2),noise.coh_ave);
    block_start = 1 + [0 cumsum(block_size(1:end-1))];
    
    % Remove bad blocks
    % =====================================================================
    % NOT DONE YET
    
    %% Fourier analysis of each bin
    % =====================================================================
    Nx = length(noise.gps_time);
    
    Nx_dft = round(Nx / param.collate_coh_noise.dft_corr_time);
    if Nx_dft<1
      Nx_dft = 1;
    end
    dft_freqs = ifftshift(-floor(Nx_dft/2) : floor((Nx_dft-1)/2));
    [~,dft_freqs_idxs] = sort(abs(dft_freqs));
    dft_freqs = dft_freqs(dft_freqs_idxs);
    noise.dft = zeros(Nt, Nx_dft, 'single');
    
    dgps_time = median(diff(noise.gps_time));
    dx = max(1,round(1/(dgps_time * param.collate_coh_noise.firdec_fs)));
    dec_idxs = 1:dx:Nx;
    noise.coh_noise_gps_time = noise.gps_time(dec_idxs);
    Nx_coh_noise = length(dec_idxs);
    noise.coh_noise = zeros(Nt, Nx_coh_noise, 'single');
    
    if enable_threshold
      threshold = zeros(Nt,1);
      cn_before_mag = zeros(Nx,Nt);
    end
    
    if enable_cn_plot
      % Transposed to speed up writes
      cn_before = zeros(Nx,Nt);
      cn_after = zeros(Nx,Nt);
    end
    for bin = start_bin:stop_bin
      bin_idx = bin-start_bin+1;
      if ~mod(bin_idx-1,10^floor(log10(Nt)-1))
        fprintf('  Estimating rbin index %d of %d (%s)\n', bin_idx, Nt, datestr(now));
      end
      coh_bin = nan(size(noise.gps_time));
      if enable_threshold
        coh_bin_mag = nan(size(noise.gps_time));
      end
      for block_idx = 1:length(noise.coh_ave)
        if bin >= start_bins(block_idx) && bin <= stop_bins(block_idx)
          tmp = noise.coh_ave{block_idx}(bin-start_bins(block_idx)+1,:);
          tmp(noise.coh_ave_samples{block_idx}(bin-start_bins(block_idx)+1,:) < cmd.min_samples) = NaN;
          coh_bin(block_start(block_idx)+(0:block_size(block_idx)-1)) = tmp;
          if enable_threshold
            coh_bin_mag(block_start(block_idx)+(0:block_size(block_idx)-1)) = noise.coh_ave_mag{block_idx}(bin-start_bins(block_idx)+1,:);
          end
        end
      end
      if 0
        % Debug plots
        figure(1); clf;
        plot(real(coh_bin)); title(sprintf('%d',bin));
        if enable_threshold
          plot(coh_bin_mag); title(sprintf('%d',bin));
        end
        figure(2); clf;
        plot(imag(coh_bin)); title(sprintf('%d',bin));
      end
      if enable_cn_plot
        cn_before(:,bin_idx) = coh_bin;
      end
      if enable_threshold
        cn_before_mag(:,bin_idx) = coh_bin_mag;
        if size(coh_bin_mag,2) < 10
          threshold(bin_idx) = lp(mean(abs(coh_bin_mag).^2,2));
        else
          threshold(bin_idx) = min(lp(fir_dec(abs(coh_bin_mag).^2,10)),[],2);
        end
      end
      if strcmpi(param.collate_coh_noise.method,'dft')
        for dft_idx = 1:length(dft_freqs)
          mf = exp(1i*2*pi/Nx * dft_freqs(dft_idx) .* (0:Nx-1));
          noise.dft(bin_idx,dft_idx) = nanmean(conj(mf).*coh_bin);
          coh_bin = coh_bin - noise.dft(bin_idx,dft_idx) * mf;
        end
      elseif strcmpi(param.collate_coh_noise.method,'firdec')
        noise.dft(bin_idx,1) = nanmean(coh_bin);
        fcutoff = param.collate_coh_noise.firdec_fcutoff(noise.dt*bin);
        if fcutoff == 0
          % Zero cutoff frequency: Take mean over all values
          noise.coh_noise(bin_idx,:) = nanmean(coh_bin);
        elseif fcutoff < 0
          % Negative cutoff frequency: Disable coherent noise removal
          noise.coh_noise(bin_idx,:) = 0;
        else
          % Positive cutoff frequency: Regular FIR filter
          B = tukeywin(round(1/(fcutoff*dgps_time)/2)*2+1,0.5).';
          B = B / sum(B);
          noise.coh_noise(bin_idx,:) = nan_fir_dec(coh_bin,B,dx);
        end
        %noise.coh_noise(bin_idx,isnan(noise.coh_noise(bin_idx,:))) = 0;
        noise_est = interp_finite(interp1(noise.coh_noise_gps_time,noise.coh_noise(bin_idx,:),noise.gps_time),0);
        coh_bin = coh_bin - noise_est;
      end
      if enable_cn_plot
        cn_after(:,bin_idx) = coh_bin;
      end
      if 0
        % Debug plots
        noise.dft(bin_idx,dft_idx)
        keyboard
      end
    end
    
    time = start_bin*noise.dt + noise.dt*(0:Nt-1).';
    Tpd = param.radar.wfs(wf).Tpd;
    if enable_threshold
      orig_threshold = threshold;
      if ~isempty(param.collate_coh_noise.threshold_eval)
        if numel(param.collate_coh_noise.threshold_eval) < wf
          error('If param.collate_coh_noise.threshold_eval is specified, there must be a cell entry for each waveform that is used. Waveform %d cannot be found since numel(param.collate_coh_noise.threshold_eval)=%d',wf,numel(param.collate_coh_noise.threshold_eval));
        end
        
        %figure(100); plot(threshold); hold on;
        % Example:
        % param.collate_coh_noise.threshold_eval{wf} = 'threshold(time>Tpd+0.85e-6 & threshold>-110) = -100; threshold(time<=Tpd+0.85e-6) = inf;'
        % param.collate_coh_noise.threshold_eval{wf} = 'threshold(time>Tpd+2.3e-6 & threshold>-130) = -110; threshold(time<=Tpd+2.3e-6) = threshold(time<=Tpd+2.3e-6)+20;';
        eval(param.collate_coh_noise.threshold_eval{wf});
      end
    end
    
    if enable_cn_plot
      cn_before = cn_before.';
      cn_after = cn_after.';
      
      dxt = mean_without_outliers(diff(noise.gps_time));
      Xt = Nx*dxt;
      dfx = 1/Xt;
      fx = dfx * (-floor(Nx/2) : floor((Nx-1)/2));
      
      mask = isnan(cn_before);
      cn_before(mask) = 0;
      
      clf(h_fig(1));
      set(h_fig(1), 'name', 'collate_coh_noise FFT');
      h_axes(1) = axes('parent',h_fig(1));
      imagesc(fx, [], fftshift(lp( fft(cn_before,[],2) ),2), 'parent', h_axes(1));
      cn_before(mask) = NaN;
      title(h_axes(1), sprintf('%s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      xlabel(h_axes(1), 'Frequency (1/m)');
      ylabel(h_axes(1), 'Range bin');
      
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_fft_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      saveas(h_fig(1),fig_fn);
      size_fig = whos('cn_before');
      if size_fig.bytes < 1e9
        fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_fft_wf_%02d_adc_%02d',wf,adc)) '.fig'];
        fprintf('Saving %s\n', fig_fn);
        saveas(h_fig(1),fig_fn);
      end
      
      %cn_before(bsxfun(@gt,lp(cn_before),threshold)) = NaN;
      clf(h_fig(2));
      set(h_fig(2), 'name', 'collate_coh_noise Before');
      h_axes(2) = axes('parent',h_fig(2));
      imagesc(lp(cn_before),'parent',h_axes(2));
      cc=caxis(h_axes(2));
      title(h_axes(2), sprintf('Before coherent noise removal %s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      ylabel(h_axes(2), 'Range bin');
      xlabel(h_axes(2), 'Block');
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      saveas(h_fig(2),fig_fn);
      if size_fig.bytes < 1e9
        fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_wf_%02d_adc_%02d',wf,adc)) '.fig'];
        fprintf('Saving %s\n', fig_fn);
        saveas(h_fig(2),fig_fn);
      end
      
      %cn_before(bsxfun(@gt,lp(cn_before),threshold)) = NaN;
      clf(h_fig(3));
      set(h_fig(3), 'name', 'collate_coh_noise Before Phase');
      h_axes(3) = axes('parent',h_fig(3));
      imagesc(angle(cn_before),'parent',h_axes(3));
      title(h_axes(3), sprintf('Before coherent noise removal %s wf %d adc %d (phase)',regexprep(param.day_seg,'_','\\_'), wf, adc));
      ylabel(h_axes(3), 'Range bin');
      xlabel(h_axes(3), 'Block');
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_phase_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      saveas(h_fig(3),fig_fn);
      if size_fig.bytes < 1e9
        fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_phase_wf_%02d_adc_%02d',wf,adc)) '.fig'];
        fprintf('Saving %s\n', fig_fn);
        saveas(h_fig(3),fig_fn);
      end
      
      clf(h_fig(4));
      set(h_fig(4), 'name', 'collate_coh_noise After');
      h_axes(4) = axes('parent',h_fig(4));
      imagesc(lp(cn_after),'parent',h_axes(4));
      caxis(h_axes(4), cc);
      title(h_axes(4), sprintf('After coherent noise removal %s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      xlabel(h_axes(4), 'Block');
      ylabel(h_axes(4), 'Range bin');
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_after_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      saveas(h_fig(4),fig_fn);
      if size_fig.bytes < 1e9
        fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_after_wf_%02d_adc_%02d',wf,adc)) '.fig'];
        fprintf('Saving %s\n', fig_fn);
        saveas(h_fig(4),fig_fn);
      end
      
      linkaxes(h_axes(2:4));
    end
    
    if enable_threshold_plot
      clf(h_fig(5));
      set(h_fig(5), 'name', 'Threshold');
      h_axes(5) = axes('parent',h_fig(5));
      plot(h_axes(5), orig_threshold)
      hold(h_axes(5), 'on');
      grid(h_axes(5), 'on');
      plot(h_axes(5), threshold)
      plot(h_axes(5), lp(abs(noise.dft(:,1)).^2))
      legend(h_axes(5), 'Original', 'Modified', 'DC Noise', 'location', 'best')
      xlabel(h_axes(5), 'Range bin');
      ylabel(h_axes(5), 'Relative power (dB)');
      title(h_axes(5), sprintf('Threshold %s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      if ~isempty(param.collate_coh_noise.threshold_ylims)
        ylim(h_axes(5),param.collate_coh_noise.threshold_ylims);
      end
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('threshold_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      fig_fn_dir = fileparts(fig_fn);
      if ~exist(fig_fn_dir,'dir')
        mkdir(fig_fn_dir);
      end
      saveas(h_fig(5),fig_fn);
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('threshold_wf_%02d_adc_%02d',wf,adc)) '.fig'];
      fprintf('Saving %s\n', fig_fn);
      saveas(h_fig(5),fig_fn);
    end
    
    if enable_visible_plot
      % Bring plots to front
      for h_fig_idx = 1:length(h_fig)
        figure(h_fig(h_fig_idx));
      end
      % Enter debug mode
      keyboard
    end
    
    %% Create the simplified output
    % =====================================================================
    noise_simp = struct('gps_time',noise.gps_time);
    noise_simp.start_bin = start_bin;
    noise_simp.dt = noise.dt;
    noise_simp.fc = noise.fc;
    if strcmpi(param.collate_coh_noise.method,'dft')
      noise_simp.dft_freqs = dft_freqs;
      noise_simp.dft = noise.dft;
    elseif strcmpi(param.collate_coh_noise.method,'firdec')
      noise_simp.coh_noise_gps_time = noise.coh_noise_gps_time;
      noise_simp.coh_noise = noise.coh_noise;
    end
    noise_simp.param_records = noise.param_records;
    noise_simp.param_analysis = noise.param_analysis;
    if enable_threshold
      noise_simp.threshold = threshold;
    end
    noise_simp.param_collate_coh_noise = param;
    noise_simp.datestr = datestr(now);
    noise_simp.recs = noise.param_analysis.analysis.block_size/2 + noise.param_analysis.analysis.block_size * (0:Nx-1);
    if param.ct_file_lock
      noise_simp.file_version = '1L';
    else
      noise_simp.file_version = '1';
    end
    
    %% Save the result
    % =====================================================================
    out_fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.out_path, ''));
    out_fn = fullfile(out_fn_dir,sprintf('coh_noise_simp_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
    fprintf('Saving %s (%s)\n', out_fn, datestr(now));
    save(out_fn,'-v7.3','-struct','noise_simp');

%         case 'custom2'
%           
%           coh_ave_RFI = noise.coh_ave(regime_mask,:).';
%           
%           % Remove RFI
%           B = [ones(1,11), 0, ones(1,11)]; A = 1; B = B / sum(B);
%           threshold = fir_dec(abs(coh_ave_RFI).^2, B, 1);
%           mask = lp(coh_ave_RFI) > lp(threshold) + 20;
%           coh_ave_RFI(mask) = 0;
%           for rbin = 1:size(coh_ave_RFI,1)
%             coh_ave_RFI(rbin,:) = interp_finite(coh_ave_RFI(rbin,:).',0).';
%           end
%           
%           coh_ave = coh_ave_RFI;
%           Mx = 1;
%           Nx = size(coh_ave,2);
%           Nx_cutoff = round(cmd.Wn*Nx);
%           
%           median_coh_ave = median(abs(coh_ave),2);
%           mask = bsxfun(@gt,abs(coh_ave),abs(median_coh_ave*cmd.threshold));
%           
%           coh_ave(mask) = NaN;
%           
%           coh_ave_hpf = coh_ave_RFI;
%           for fc = -round(cmd.Wn*Nx):round(cmd.Wn*Nx)
%             coh_ave = coh_ave_hpf;
%             coh_ave(mask) = NaN;
%             coh_ave_mod = bsxfun(@times,coh_ave,exp(-1i*2*pi*fc/Nx*(0:Nx-1)));
%             mean_coh_ave = nanmean(coh_ave_mod,2);
%             mean_coh_ave = kron(mean_coh_ave, exp(1i*2*pi*fc/Nx*(0:Nx-1)));
%             coh_ave_hpf = coh_ave_hpf - mean_coh_ave;
%           end
%           
%           if 0
%             figure(1); clf;
%             imagesc(lp(coh_ave_hpf))
%             title('Should have coherent noise removed if working well.')
%             caxis([0 100]);
%             
%             figure(2); clf;
%             imagesc(lp(noise.coh_ave(regime_mask,:).'))
%             caxis([0 100]);
%             
%             figure(3); clf;
%             imagesc(mask)
%             
%             keyboard
%           end
%           
%           % Store data (1 - HPF = LPF)
%           noise.coh_ave(regime_mask,:) = coh_ave_RFI.'-coh_ave_hpf.';
%           
    
  end
end

if ~enable_visible_plot
  try
    delete(h_fig);
  end
end
