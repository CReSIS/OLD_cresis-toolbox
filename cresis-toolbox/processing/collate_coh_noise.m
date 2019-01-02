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
cmd = param.analysis.cmd{param.collate_coh_noise.cmd_idx};

if ~isfield(param.collate_coh_noise,'imgs') || isempty(param.collate_coh_noise.imgs)
  param.collate_coh_noise.imgs = 1:length(param.analysis.imgs);
end

if ~isfield(param.collate_coh_noise,'in_dir') || isempty(param.collate_coh_noise.in_dir)
  param.collate_coh_noise.in_dir = 'analysis';
end

if ~isfield(param.collate_coh_noise,'out_dir') || isempty(param.collate_coh_noise.out_dir)
  param.collate_coh_noise.out_dir = 'analysis';
end

if ~isfield(param.collate_coh_noise,'plot_en') || isempty(param.collate_coh_noise.plot_en)
  param.collate_coh_noise.plot_en = false;
end

if ~isfield(param.collate_coh_noise,'threshold_eval') || isempty(param.collate_coh_noise.threshold_eval)
  % Default is no evaluation
  % Example of how to use evaluation:
  %   params(param_idx).collate_coh_noise.threshold_eval{wf} = 'tmp=threshold(max(1,round(1.1*Tpd_bin)+60):end); tmp(tmp>-110)=-110; threshold(max(1,round(1.1*Tpd_bin)+60):end)=tmp; threshold=threshold+10;';
  param.collate_coh_noise.threshold_eval = {};
end

if ~isfield(param.collate_coh_noise,'wf_adcs') || isempty(param.collate_coh_noise.wf_adcs)
  param.collate_coh_noise.wf_adcs = [];
end

if ~isfield(cmd,'dft_corr_length') || isempty(cmd.dft_corr_length)
  cmd.dft_corr_length = 1e6;
end

if ~isfield(cmd,'min_samples') || isempty(cmd.min_samples)
  cmd.min_samples = 0.5*cmd.block_ave;
end

if ~isfield(cmd,'power_grow') || isempty(cmd.power_grow)
  cmd.power_grow = 1;
end

if param.collate_coh_noise.plot_en
  fig_visible = 'off';
  h_fig = figure('Visible',fig_visible);
  h_axes = axes('parent',h_fig);
  for idx=2:4
    h_fig(idx) = figure('Visible',fig_visible);
    h_axes(idx) = axes('parent',h_fig(idx));
  end
end

for img = param.collate_coh_noise.imgs
  
  if isempty(param.collate_coh_noise.wf_adcs)
    wf_adcs = 1:size(param.analysis.imgs{img},1);
  else
    wf_adcs = param.collate_coh_noise.wf_adcs;
  end
  for wf_adc = wf_adcs
    wf = param.analysis.imgs{img}(wf_adc,1);
    adc = param.analysis.imgs{img}(wf_adc,2);
    
    %% Load the coherent noise file
    % =====================================================================
    fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.in_dir));
    fn = fullfile(fn_dir,sprintf('coh_noise_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
    fprintf('Loading %s (%s)\n', fn, datestr(now));
    noise = load(fn);
    
    % Each processed block has a constant start and stop bin, but between
    % blocks the start/stop range bin may change.
    start_bins = round(noise.t0/noise.dt);
    stop_bins = round(noise.t0/noise.dt) + noise.Nt - 1;
    % Determine the first start bin and the last stop bin
    start_bin = min(start_bins);
    stop_bin = max(stop_bins);
      
    % Get the size of each processing block
    block_size = cellfun(@(t) size(t,2),noise.coh_ave);
    block_start = 1 + [0 cumsum(block_size(1:end-1))];
    
    % Remove bad blocks
    % =====================================================================
    % NOT DONE YET
    
    %% Fourier analysis of each bin
    % =====================================================================
    Nx = length(noise.gps_time);
    Nx_dft = round(Nx / cmd.dft_corr_length);
    if Nx_dft<1
      Nx_dft = 1;
    end
    dft_freqs = ifftshift(-floor(Nx_dft/2) : floor((Nx_dft-1)/2));
    [~,dft_freqs_idxs] = sort(abs(dft_freqs));
    dft_freqs = dft_freqs(dft_freqs_idxs);
    noise.dft = zeros(stop_bin-start_bin+1, Nx_dft, 'single');
    threshold = zeros(stop_bin-start_bin+1,1);
    if param.collate_coh_noise.plot_en
      cn_before = zeros(stop_bin-start_bin+1,Nx);
      cn_after = zeros(stop_bin-start_bin+1,Nx);
      cn_before_mag = zeros(stop_bin-start_bin+1,Nx);
    end
    for bin = start_bin:stop_bin
      bin_idx = bin-start_bin+1;
      coh_bin = nan(size(noise.gps_time));
      coh_bin_mag = nan(size(noise.gps_time));
      for block_idx = 1:length(noise.coh_ave)
        if bin >= start_bins(block_idx) && bin <= stop_bins(block_idx)
          tmp = noise.coh_ave{block_idx}(bin-start_bins(block_idx)+1,:);
          tmp(noise.coh_ave_samples{block_idx}(bin-start_bins(block_idx)+1,:) < cmd.min_samples) = NaN;
          coh_bin(block_start(block_idx)+(0:block_size(block_idx)-1)) = tmp;
          coh_bin_mag(block_start(block_idx)+(0:block_size(block_idx)-1)) = noise.coh_ave_mag{block_idx}(bin-start_bins(block_idx)+1,:);
        end
      end
      if 0
        % Debug plots
        figure(1); clf;
        plot(real(coh_bin)); title(sprintf('%d',bin));
        plot(real(coh_bin_mag)); title(sprintf('%d',bin));
        figure(2); clf;
        plot(imag(coh_bin)); title(sprintf('%d',bin));
        plot(imag(coh_bin_mag)); title(sprintf('%d',bin));
      end
      if param.collate_coh_noise.plot_en
        cn_before(bin_idx,:) = coh_bin;
        cn_before_mag(bin_idx,:) = coh_bin_mag;
      end
      if size(coh_bin_mag,2) < 10
        threshold(bin_idx) = lp(mean(abs(coh_bin_mag).^2,2));
      else
        threshold(bin_idx) = min(lp(fir_dec(abs(coh_bin_mag).^2,10)),[],2);
      end
      for dft_idx = 1:length(dft_freqs)
        mf = exp(1i*2*pi/Nx * dft_freqs(dft_idx) .* (0:Nx-1));
        noise.dft(bin_idx,dft_idx) = nanmean(conj(mf).*coh_bin);
        coh_bin = coh_bin - noise.dft(bin_idx,dft_idx) * mf;
      end
      if param.collate_coh_noise.plot_en
        cn_after(bin_idx,:) = coh_bin;
      end
      if 0
        % Debug plots
        noise.dft(bin_idx,dft_idx)
        keyboard
      end
    end
    orig_threshold = threshold;
    if ~isempty(param.collate_coh_noise.threshold_eval)
      if numel(param.collate_coh_noise.threshold_eval) < wf
        error('If param.collate_coh_noise.threshold_eval is specified, there must be a cell entry for each waveform that is used. Waveform %d cannot be found since numel(param.collate_coh_noise.threshold_eval)=%d',wf,numel(param.collate_coh_noise.threshold_eval));
      end

      Tpd_bin  = round(param.radar.wfs(wf).Tpd/noise.dt) - start_bin;
      %figure(100); plot(threshold); hold on;
      % Example:
      %param.collate_coh_noise.threshold_eval{wf} = 'tmp=threshold(max(1,round(1.1*Tpd_bin)+60):end); tmp(tmp>-110)=-110; threshold(max(1,round(1.1*Tpd_bin)+60):end)=tmp; threshold=threshold+10;';
      eval(param.collate_coh_noise.threshold_eval{wf});
    end
    
    if param.collate_coh_noise.plot_en
      dxt = mean_without_outliers(diff(noise.gps_time));
      Xt = Nx*dxt;
      dfx = 1/Xt;
      fx = dfx * (-floor(Nx/2) : floor((Nx-1)/2));
      
      for idx=1:4
        delete(get(h_axes(idx),'Children'));
      end
      
      mask = isnan(cn_before);
      cn_before(mask) = 0;
      imagesc(fx, [], fftshift(lp( fft(cn_before,[],2) ),2), 'parent', h_axes(1));
      cn_before(mask) = NaN;
      title(h_axes(1), sprintf('%s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
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
      
      cn_before(bsxfun(@gt,lp(cn_before),threshold)) = NaN;
      imagesc(lp(cn_before),'parent',h_axes(2));
      title(h_axes(2), sprintf('Before %s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
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
      
      imagesc(lp(cn_after),'parent',h_axes(3));
      cc=caxis(h_axes(2));
      caxis(h_axes(3), cc);
      title(h_axes(3), sprintf('After %s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      xlabel(h_axes(3), 'Block');
      ylabel(h_axes(3), 'Range bin');
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_after_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      saveas(h_fig(3),fig_fn);
      if size_fig.bytes < 1e9
        fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('coh_after_wf_%02d_adc_%02d',wf,adc)) '.fig'];
        fprintf('Saving %s\n', fig_fn);
        saveas(h_fig(3),fig_fn);
      end
      
      linkaxes(h_axes(2:3));
      
      plot(h_axes(4), orig_threshold)
      hold(h_axes(4), 'on');
      grid(h_axes(4), 'on');
      plot(h_axes(4), threshold)
      plot(h_axes(4), lp(abs(noise.dft(:,1)).^2))
      legend(h_axes(4), 'Original', 'Modified', 'DC Noise', 'location', 'best')
      xlabel(h_axes(4), 'Range bin');
      ylabel(h_axes(4), 'Relative power (dB)');
      title(h_axes(4), sprintf('Threshold %s wf %d adc %d',regexprep(param.day_seg,'_','\\_'), wf, adc));
      param.collate_coh_noise.threshold_ylims = [-180 -20];
      ylim(h_axes(4),param.collate_coh_noise.threshold_ylims);
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('threshold_wf_%02d_adc_%02d',wf,adc)) '.jpg'];
      fprintf('Saving %s\n', fig_fn);
      saveas(h_fig(4),fig_fn);
      fig_fn = [ct_filename_ct_tmp(param,'','collate_coh_noise',sprintf('threshold_wf_%02d_adc_%02d',wf,adc)) '.fig'];
      fprintf('Saving %s\n', fig_fn);
      saveas(h_fig(4),fig_fn);

      %keyboard
    end
    
    %% Create the simplified output
    % =====================================================================
    noise_simp = struct('gps_time',noise.gps_time);
    noise_simp.start_bin = start_bin;
    noise_simp.dt = noise.dt;
    noise_simp.fc = noise.fc;
    noise_simp.dft_freqs = dft_freqs;
    noise_simp.dft = noise.dft;
    noise_simp.param_analysis = noise.param_analysis;
    noise_simp.threshold = threshold;
    noise_simp.param_collate = param;
    noise_simp.datestr = datestr(now);
    noise_simp.recs = noise.param_analysis.analysis.block_size/2 + noise.param_analysis.analysis.block_size * (0:Nx-1);
    if param.ct_file_lock
      noise_simp.file_version = '1L';
    else
      noise_simp.file_version = '1';
    end
    
    %% Save the result
    % =====================================================================
    out_fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.out_dir, ''));
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
if param.collate_coh_noise.plot_en
  delete(h_fig);
end
