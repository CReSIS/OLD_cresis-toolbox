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

cmd = param.analysis.cmd{param.collate_coh_noise.cmd_idx};
debug_level = param.collate_coh_noise.debug_level;

if ~isfield(param.collate_coh_noise,'imgs') || isempty(param.collate_coh_noise.imgs)
  param.collate_coh_noise.imgs = 1:length(param.analysis.imgs);
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

for img = param.collate_coh_noise.imgs
  
  if isempty(param.collate_coh_noise.wf_adcs)
    wf_adcs = 1:size(param.analysis.imgs{img},1);
  else
    wf_adcs = param.collate_coh_noise.wf_adcs;
  end
  for wf_adc = wf_adcs
    wf = param.analysis.imgs{img}(wf_adcs,1);
    adc = param.analysis.imgs{img}(wf_adcs,2);
    
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
    for bin = start_bin:stop_bin
      bin_idx = bin-start_bin+1;
      data_bin = nan(size(noise.gps_time));
      for block_idx = 1:length(noise.coh_ave)
        if bin >= start_bins(block_idx) && bin <= stop_bins(block_idx)
          data_bin(block_start(block_idx)+(0:block_size(block_idx)-1)) = noise.coh_ave{block_idx}(bin-start_bins(block_idx)+1,:);
        end
      end
      %figure(1); clf;
      %plot(real(data_bin)); title(sprintf('%d',bin));
      %figure(2); clf;
      %plot(imag(data_bin)); title(sprintf('%d',bin));
      for dft_idx = 1:length(dft_freqs)
        mf = exp(1i*2*pi/Nx * dft_freqs(dft_idx) .* (0:Nx-1));
        noise.dft(bin_idx,dft_idx) = nanmean(mf.*data_bin);
        data_bin = data_bin - noise.dft(bin_idx,dft_idx) * mf;
      end
      %noise.dft(bin,dft_idx)
      %keyboard
    end
    
    %% Create the simplified output
    % =====================================================================
    noise_simp = struct('gps_time',noise.gps_time);
    noise_simp.start_bin = start_bin;
    noise_simp.dt = noise.dt;
    noise_simp.dft_freqs = dft_freqs;
    noise_simp.dftI = real(noise.dft);
    noise_simp.dftQ = imag(noise.dft);
    noise_simp.param_analysis = noise.param_analysis;
    noise_simp.param_collate = param;
    noise_simp.datestr = datestr(now);
    noise_simp.recs = noise.param_analysis.analysis.block_size/2 + noise.param_analysis.analysis.block_size * (0:Nx-1);
    noise_simp.file_version = '1';
    
    %% Store the simplified output in netcdf file
    % =====================================================================
    out_fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.out_dir, ''));
    out_fn = fullfile(out_fn_dir,sprintf('coh_noise_simp_%s_wf_%d_adc_%d.nc', param.day_seg, wf, adc));
    fprintf('Saving %s (%s)\n', out_fn, datestr(now));
    netcdf_from_mat(out_fn,noise_simp);

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
