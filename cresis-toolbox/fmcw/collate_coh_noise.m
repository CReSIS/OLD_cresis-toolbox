function collate_coh_noise(param,param_override)
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

% HACK TO SUPPORT OLD FORMAT... WILL BE REMOVED EVENTUALLY
old_format = 0;
if old_format
  cmd.method = 'sgolayfilt';
  cmd.K = param.get_heights.coh_noise_arg{1};
  cmd.F = param.get_heights.coh_noise_arg{2};
  cmd.W = param.get_heights.coh_noise_arg{3};
end

if ~isfield(cmd,'doppler_window') ...
    || isempty(cmd.doppler_window)
  doppler_window = hanning(61); doppler_window = doppler_window(1:30); % experimental
  cmd.doppler_window = doppler_window;
else
  cmd.doppler_window = cmd.doppler_window(1:floor(length(cmd.doppler_window)/2));
  doppler_window = cmd.doppler_window;
end

if ~isfield(cmd,'doppler_threshold') || isempty(cmd.doppler_threshold)
  cmd.doppler_threshold = inf;
end

if ~isfield(cmd,'doppler_grow') || isempty(cmd.doppler_grow)
  cmd.doppler_grow = 0;
end

if ~isfield(cmd,'min_samples') || isempty(cmd.min_samples)
  cmd.min_samples = 0.5*cmd.block_ave;
end

if ~isfield(cmd,'power_grow') || isempty(cmd.power_grow)
  cmd.power_grow = 1;
end

%% Load the coherent noise file
% =====================================================================
fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.in_dir));
fn = fullfile(fn_dir,sprintf('coh_noise_%s_img_%02d.mat', param.day_seg, param.collate_coh_noise.img));
fprintf('Loading %s\n', fn);
noise = load(fn);

%% Debug
% =====================================================================
if 0
  records = load(ct_filename_support(param,'','records'));
  load(ct_filename_support(param,'','frames'));
  
  aa=noise.coh_ave;
  aa=fftshift(aa,1);
  aa(:,any(lp(aa(100:end-100,:))>40)) = NaN;
  aa=fftshift(aa,1);
  aa(:,sum(isnan(noise.coh_ave))>10) = NaN;
  aa(noise.coh_ave_samples<500)=NaN;
  figure(1); clf;
  imagesc(lp(aa))
  figure(2); clf;
  imagesc(lp(noise.coh_ave))
  figure(1); h_axes=gca; figure(2); h_axes(end+1) = gca; linkaxes(h_axes,'xy');
  return;
end

%% Create the Doppler mask
% =====================================================================
doppler_psd = lp(nanmean(noise.doppler,2));
doppler_psd = interp_finite(doppler_psd,NaN);
doppler_noise_floor = medfilt1(double(doppler_psd),201);
doppler_mask = doppler_psd > doppler_noise_floor + cmd.doppler_threshold;
doppler_mask(1) = 0; % Regular coherent noise removal removes the DC component
doppler_mask = grow(doppler_mask,length(doppler_window) + cmd.doppler_grow);
if debug_level > 1
  % Debug code
  figure(1); clf;
  plot(doppler_psd)
  hold on;
  plot(lp(nanmax(noise.doppler,[],2)),'c')
  plot(doppler_noise_floor + cmd.doppler_threshold,'r')
  plot(find(doppler_mask), doppler_psd(doppler_mask),'k.');
  fprintf('Adjust Doppler mask by running commands like "doppler_mask(1:120) = 0;"\n');
  keyboard
end

% For each masked out section, weight the edges
doppler_weights = zeros(size(doppler_mask));
weight_state = 1-doppler_mask(1);
doppler_weights(1) = weight_state;
idx = 2;
while idx <= length(doppler_weights)
  if ~doppler_mask(idx-1) && doppler_mask(idx)
    % Transition from good to bad
    doppler_weights(idx + (0:length(doppler_window)-1)) = doppler_window(end:-1:1);
    weight_state = 0;
    idx = idx + length(doppler_window);
    
  elseif doppler_mask(idx-1) && ~doppler_mask(idx)
    % Transition from bad to good
    doppler_weights(idx + (-length(doppler_window):-1)) = doppler_window;
    doppler_weights(idx) = 1;
    weight_state = 1;
    idx = idx + 1;
    
  else
    doppler_weights(idx) = weight_state;
    idx = idx + 1;
  end
end
if 0
  figure(1); clf;
  plot(doppler_mask)
  hold on;
  plot(doppler_weights,'r')
  keyboard
end

%% Set the noise regimes
% =====================================================================
noise.regime = ones(size(noise.gps_time));
if strcmpi(noise.param_analysis.radar_name,'snow2') && strcmpi(noise.param_analysis.day_seg,'20120316_03')
  % Segment specific regimes
  bad_mask = noise.gps_time > 1.331923763327688e+09 & noise.gps_time < 1.331923803347247e+09;
  noise.coh_ave_samples(:,bad_mask) = cmd.min_samples;
  noise.regime(find(bad_mask,1):end) = 2;
  
elseif isfield(cmd,'regimes') ...
    && ~isempty(cmd.regimes) ...
    && cmd.regimes.en
  % Cross correlate neighboring range lines to test for changes in statistics which indciate a regime change
  dline = 2;
  dcorr = zeros(1,size(noise.coh_ave,2));
  tmp = interp_finite(noise.coh_ave);
  for rline = 1:size(noise.coh_ave,2)-dline
    dcorr(rline+1) = norm(tmp(:,rline)-tmp(:,rline+2)) ./ norm(tmp(:,rline));
  end
  % Deal with edges
  dcorr(1) = dcorr(2);
  dcorr(end) = dcorr(end-1);
  
  if debug_level > 0
    figure(1); clf;
    imagesc(lp(noise.coh_ave));
    h_axis = gca;
    
    figure(2); clf;
    plot(dcorr)
    h_axis(end+1) = gca;
    
    linkaxes(h_axis,'x')
  end
  
  % Statistic changes vector
  stat_change = lp(dcorr) > cmd.regimes.threshold;
  %noise.coh_ave(:,stat_change == 1) = NaN;
  
  % Create different noise regimes for each statistics change
  cur_regime = 1;
  for rline = 2:length(stat_change)
    if stat_change(rline) && ~stat_change(rline-1)
      cur_regime = cur_regime + 1;
    end
    noise.regime(rline) = cur_regime;
  end
elseif isfield(noise,'nyquist_zone')
  %% Segment into regimes based on nyquist zone
  noise.regime = noise.nyquist_zone;
  
  % Nyquist_zone: bit mask indicating which nyquist zones are used in each
  % block. If only one nyquist zone is used, then only one bit will be
  % one. This means the values will be 1, 2, 4, or 8.  In other words it
  % will be a power of 2 and log2() will produce an integer. If two or
  % more nyquist zones are used then log2() will not produce an integer.
  % If two or more nyquist zones are used, then the mean/average will be
  % contaminated and should not be used and we mark that column with 1 in
  % mixed_nz_blocks and set these columns to zero valid samples.
  mixed_nz_blocks = mod(log2(noise.nyquist_zone),1)~=0;
  noise.coh_ave_samples(:,mixed_nz_blocks) = 0;
  
end
regimes = unique(noise.regime);
if length(regimes) > 1
  warning('There are %d noise regimes\n', length(regimes));
end

%% Apply the filtering from coh_noise_arg across each regime
% =====================================================================
noise.coh_ave_samples = uint32(noise.coh_ave_samples);
noise.coh_ave = single(noise.coh_ave);
old_noise_coh_ave = noise.coh_ave;
if debug_level > 0
  figure(1); clf;
  imagesc(lp(noise.coh_ave));
  aa = gca;
end

noise.coh_ave(noise.coh_ave_samples <= cmd.min_samples) = NaN;

if debug_level > 0
  figure(2); clf;
  imagesc(lp(noise.coh_ave))
  aa(2) = gca;
end

for rline = 1:size(noise.coh_ave,2)
  noise.coh_ave(filter2(cmd.power_grow,double(isnan(noise.coh_ave(:,rline)))) > 0) = NaN;
end

if debug_level > 0
  figure(3); clf;
  imagesc(lp(noise.coh_ave))
  aa(3) = gca;
  linkaxes(aa,'xy')
  keyboard
end

noise.coh_ave = noise.coh_ave.';
for regime = regimes
  regime_mask = find(noise.regime == regime);
  if any(all(isnan(noise.coh_ave(regime_mask,:))))
    regime_fill = find(all(isnan(noise.coh_ave(regime_mask,:)),1));
    noise.coh_ave(regime_mask,regime_fill) = old_noise_coh_ave(regime_fill,regime_mask).';
  end
  for rbin = 1:size(noise.coh_ave,2)
    noise.coh_ave(regime_mask,rbin) = interp_finite(noise.coh_ave(regime_mask,rbin),0);
  end
  
  switch cmd.method
    case 'DC'
      noise.coh_ave(regime_mask,:) = repmat(mean(noise.coh_ave(regime_mask,:)), [length(regime_mask) 1]);
      
    case 'window'
      coh_ave = noise.coh_ave(regime_mask,:).';
      
      % Remove RFI
      B = [ones(1,11), 0, ones(1,11)]; A = 1; B = B / sum(B);
      threshold = fir_dec(abs(coh_ave).^2, B, 1);
      mask = lp(coh_ave) > lp(threshold) + 20;
      coh_ave(mask) = 0;
      for rbin = 1:size(coh_ave,1)
        coh_ave(rbin,:) = interp_finite(coh_ave(rbin,:).',0).';
      end
      
      % HPF
      coh_ave = fft(coh_ave,[],2);
      Nx = size(coh_ave,2);
      Nx_cutoff = max(1,round(Nx*cmd.Wn));
      H = hanning(2*Nx_cutoff-1); H = H(1:(length(H)-1)/2).';
      H = [zeros(1,Nx_cutoff) H];
      coh_ave(:,1:length(H)) = bsxfun(@times, coh_ave(:,1:length(H)), H);
      coh_ave(:,end-length(H)+2:end) = bsxfun(@times, coh_ave(:,end-length(H)+2:end), H(end:-1:2));
      coh_ave = ifft(coh_ave,[],2);
      
      if 0
        figure(1); clf;
        imagesc(lp(coh_ave));
        
        dd = noise.coh_ave(regime_mask,:)-coh_ave.';

        figure(1); clf;
        imagesc(lp(noise.coh_ave(regime_mask,:)-dd).')
        
        keyboard
      end
      
      % Store data (1 - HPF = LPF)
      noise.coh_ave(regime_mask,:) = noise.coh_ave(regime_mask,:)-coh_ave.';
      
    case 'filter'
      if length(cmd.A_coh_noise) > 1
        % Use filtfilt
        noise.coh_ave(regime_mask,:) = single(filtfilt(cmd.B_coh_noise, ...
          cmd.A_coh_noise, double(noise.coh_ave(regime_mask,:))));
      else
        % Use fir_dec (no feedback)
        noise.coh_ave(regime_mask,:) = fir_dec(noise.coh_ave(regime_mask,:).',cmd.B_coh_noise,1).';
      end
      
    case 'threshold'
      Nx = size(noise.coh_ave,1);
      Mx = 10;
      Nx_cutoff = round(cmd.Wn*Mx*Nx);
      Nx_buffer = round(3/cmd.Wn);
      
      coh_ave = noise.coh_ave(regime_mask,:).';
      min_val = mean(abs(coh_ave(:,1:max(1,round(Nx_cutoff/Mx)))),2);
      qq = bsxfun(@times,min_val,exp(-1i*2*angle(bsxfun(@times,coh_ave(:,2:1+Nx_buffer),conj(coh_ave(:,1))))));
      min_val = mean(abs(coh_ave(:,end-max(1,round(Nx_cutoff/Mx))+1:end )),2);
      qq2 = bsxfun(@times,min_val,exp(-1i*2*angle(bsxfun(@times,coh_ave(:,end-Nx_buffer:end-1),conj(coh_ave(:,end))))));
      rr = [fliplr(bsxfun(@times,exp(1i*angle(coh_ave(:,2:1+Nx_buffer))),qq)), coh_ave, fliplr(bsxfun(@times,exp(1i*angle(coh_ave(:,end-Nx_buffer:end-1))),qq2))];
      Nx = size(rr,2);
      
      % Remove RFI
      B = [ones(1,11), 0, ones(1,11)]; A = 1; B = B / sum(B);
      tt = fir_dec(abs(rr).^2, B, 1);
      mask = lp(rr) > lp(tt) + 20;
      rr(mask) = 0;
      for rbin = 1:size(rr,1)
        rr(rbin,:) = interp_finite(rr(rbin,:).',0).';
      end
      
      dt = median(diff(noise.gps_time));
      Nt = size(rr,1);
      T = dt*Nx;
      df = 1/T/Mx;
      freq = df*(0:Nx*Mx-1);
      dd = fft(rr,Mx*Nx,2);
      
      % Create threshold criteria
      ee = lp(mean(abs(dd(:,Nx_cutoff*1:Nx_cutoff*4)).^2,2));
      ee = sgolayfilt(double(ee),3,21);
      
      % Positive frequencies
      mask = zeros(size(dd));
      Nt_cutoff = 1;
      ff=dd(Nt_cutoff:end,1:Nx_cutoff);
      mask = bsxfun(@gt, lp(ff), ee(Nt_cutoff:end) + cmd.threshold);
%       mask(:,1:Mx) = 1;
%       
%       mask = fir_dec(double(~mask).',ones(1,11)/11,1).';
%       mask(mask<0.5) = 0; mask(mask>=0.5) = 1;
%       mask = grow(mask,cmd.grow_shrink(1),8);
%       mask = ~shrink(mask,cmd.grow_shrink(2));
%       mask = fir_dec(double(~mask).',ones(1,11)/11,1).';
%       mask(mask<0.8) = 0; mask(mask>=0.8) = 1;
%       mask = ~shrink(mask,2);
      
      ff(mask) = 0;
      dd(Nt_cutoff:end,1:Nx_cutoff) = ff;
      
      % Negative frequencies
      ff=dd(Nt_cutoff:end,end-Nx_cutoff+1:end);
      mask = bsxfun(@gt, lp(ff), ee(Nt_cutoff:end) + cmd.threshold);
%       mask(:,end-Mx+2:end) = 1;
%       
%       mask = fir_dec(double(~mask).',ones(1,11)/11,1).';
%       mask(mask<0.5) = 0; mask(mask>=0.5) = 1;
%       mask = grow(mask,cmd.grow_shrink(1),8);
%       mask = ~shrink(mask,cmd.grow_shrink(2));
%       mask = fir_dec(double(~mask).',ones(1,11)/11,1).';
%       mask(mask<0.8) = 0; mask(mask>=0.8) = 1;
%       mask = ~shrink(mask,2);
      
      ff(mask) = 0;
      dd(Nt_cutoff:end,end-Nx_cutoff+1:end) = ff;
      
      if 1
        figure(1); clf;
        imagesc(lp(dd));
        fprintf('%s\t%d\n', param.day_seg, Nx-2*Nx_cutoff);
        xlim([1 2*Nx_cutoff]);
        ylim([Nt_cutoff Nt]);
        
        figure(2); clf;
        imagesc(lp(dd));
        fprintf('%s\t%d\n', param.day_seg, Nx-2*Nx_cutoff);
        xlim([size(dd,2)-Nx_cutoff*2 size(dd,2)]);
        ylim([Nt_cutoff Nt]);
      end
      
      dd = ifft(dd,[],2);
      dd = dd(:,1:Nx);
      
      dd = dd(:,1+Nx_buffer:end-Nx_buffer);
      rr = rr(:,1+Nx_buffer:end-Nx_buffer);
      
      coh_ave = rr-dd;
      
      if 0
        debug_fig_offset = 0;
        incoh = noise.coh_ave(regime_mask,:).'-coh_ave;
        figure(3+debug_fig_offset); clf;
        imagesc(lp(incoh));

        
        figure(4+debug_fig_offset); clf;
        imagesc(lp(incoh(Nt_cutoff:Nt,[1:100,end-99:end])));
        
        keyboard
      end
      
      noise.coh_ave(regime_mask,:) = coh_ave.';
      
    case 'sgolayfilt'
      if size(noise.coh_ave(regime_mask,:),1) < cmd.F+2
        sgolayfilt_F = size(noise.coh_ave(regime_mask,:),1);
        if mod(sgolayfilt_F,2)==0
          sgolayfilt_F = sgolayfilt_F - 1;
        end
        sgolayfilt_degree = min(cmd.K, sgolayfilt_F-1);
        noise.coh_ave(regime_mask,:) = single(sgolayfilt(double(noise.coh_ave(regime_mask,:)),sgolayfilt_degree,sgolayfilt_F));
      else
        %    noise.coh_ave(regime_mask,:) = single(sgolayfilt(double(noise.coh_ave(regime_mask,:)),cmd.K,cmd.F,cmd.W));
        regime_mask_tmp = regime_mask(2:end-1);
        noise.coh_ave(regime_mask_tmp,:) = single(sgolayfilt(double(noise.coh_ave(regime_mask_tmp,:)),cmd.K,cmd.F,cmd.W));
      end
      if length(regime_mask) >= 3
        noise.coh_ave(regime_mask(1),:) = noise.coh_ave(regime_mask(2),:);
        noise.coh_ave(regime_mask(end),:) = noise.coh_ave(regime_mask(end-1),:);
      end
      
  end
  
end
clear old_noise_coh_ave;

%% Create the simplified output
% =====================================================================
noise_simp = struct('gps_time',noise.gps_time);
noise_simp.coh_aveI = real(noise.coh_ave);
noise_simp.coh_aveQ = imag(noise.coh_ave);
noise_simp.doppler_weights = doppler_weights;
noise_simp.sw_version = param.sw_version;
noise_simp.param_collate = cmd;
noise_simp.datestr = datestr(now);

%% Store the simplified output in netcdf file
% =====================================================================
out_fn_dir = fileparts(ct_filename_out(param,param.collate_coh_noise.out_dir, ''));
out_fn = fullfile(out_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
fprintf('Saving %s\n', out_fn);
netcdf_from_mat(out_fn,noise_simp);

return


%% Example code to load netcdf file
out_fn_dir = ct_filename_out(param,'', 'CSARP_noise');
out_segment_fn_dir = fileparts(out_fn_dir);
cdf_fn = fullfile(out_segment_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));

finfo = ncinfo(cdf_fn);
% Determine number of records and set recs(1) to this
Nt = finfo.Variables(find(strcmp('coh_aveI',{finfo.Variables.Name}))).Size(2);

noise = [];
noise.gps_time = ncread(cdf_fn,'gps_time');
recs = find(noise.gps_time > records.gps_time(1) - 100 & noise.gps_time < records.gps_time(end) + 100);
noise.gps_time = noise.gps_time(recs);

noise.coh_ave = ncread(cdf_fn,'coh_aveI',[recs(1) 1],[recs(end)-recs(1)+1 Nt]) ...
  + j*ncread(cdf_fn,'coh_aveQ',[recs(1) 1],[recs(end)-recs(1)+1 Nt]);

