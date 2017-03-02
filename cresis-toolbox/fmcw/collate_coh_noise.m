% script collate_coh_noise
%
% Collects coh_noise_tracker.m results from coherent noise tracking
% and creates files for removing the coherent noise.
%
% Example:
%  See run_collate_coh_noise for how to run.
%
% Authors: John Paden

%% AUTOMATED SECTION
% =========================================================================
% Loads all the coh_noise_* files and creates coh_noise_simp_* files that
% are pre-filtered for speed and saves as netcdf so that subsets of the
% files can be loaded efficiently.
% =========================================================================

params = read_param_xls(param_fn,day_seg,analysis_sheet_name);
for param_idx = 1:length(params)
  param = params(param_idx);
  
  %% Is this segment selected in the param spreadsheet
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  %% Segment must have selected this type of coherent noise removal
  if ~any(param.get_heights.coh_noise_method == [7 9 17 19])
    continue;
  end
  
  if ~isfield(param.analysis.coh_ave,'doppler_window') ...
      || isempty(param.analysis.coh_ave.doppler_window)
    doppler_window = hanning(61); doppler_window = doppler_window(1:30); % experimental
    param.analysis.coh_ave.doppler_window = doppler_window;
  else
    param.analysis.coh_ave.doppler_window = param.analysis.coh_ave.doppler_window(1:floor(length(param.analysis.coh_ave.doppler_window)/2));  
    doppler_window = param.analysis.coh_ave.doppler_window;
  end
  
  %% Get the coherent noise removal arguments for this segment
  param.proc.coh_noise_arg = param.get_heights.coh_noise_arg;
  
  %% Load the coherent noise file
  fn_dir = fileparts(ct_filename_out(param,coh_ave_file_input_type, ''));
  fn = fullfile(fn_dir,sprintf('coh_noise_img_%02d_%s.mat', img,param.day_seg));
  fprintf('collate_coh_noise %s: %s\n', param.day_seg, fn);
  noise = load(fn);
  
  %% Create the Doppler mask
  doppler_psd = lp(mean(noise.doppler,2));
  doppler_psd = interp_finite(doppler_psd);
  doppler_noise_floor = medfilt1(double(doppler_psd),201);
  doppler_mask = doppler_psd > doppler_noise_floor + param.analysis.coh_ave.doppler_threshold;
  doppler_mask(1) = 0; % Regular coherent noise removal removes the DC component
  doppler_mask = grow(doppler_mask,length(doppler_window) + param.analysis.coh_ave.doppler_grow);
  if debug_level > 1
    % Debug code
    figure(1); clf;
    plot(doppler_psd)
    hold on;
    plot(lp(max(noise.doppler,[],2)),'c')
    plot(doppler_noise_floor + param.analysis.coh_ave.doppler_threshold,'r')
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
  
  noise.regime = ones(size(noise.gps_time));
  
  %% Segment specific hacks
  if strcmpi(noise.param_analysis.radar_name,'snow2') && strcmpi(noise.param_analysis.day_seg,'20120316_03')
    bad_mask = noise.gps_time > 1.331923763327688e+09 & noise.gps_time < 1.331923803347247e+09;
    noise.coh_ave_samples(:,bad_mask) = param.analysis.coh_ave.min_samples;
    noise.regime(find(bad_mask,1):end) = 2;
    
  elseif isfield(param.analysis.coh_ave,'regimes') ...
      && ~isempty(param.analysis.coh_ave.regimes) ...
      && param.analysis.coh_ave.regimes.en
    %% Cross correlate neighboring range lines to test for changes in statistics
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
    stat_change = lp(dcorr) > param.analysis.coh_ave.regimes.threshold;
    %noise.coh_ave(:,stat_change == 1) = NaN;
    
    % Create different noise regimes for each statistics change
    cur_regime = 1;
    for rline = 2:length(stat_change)
      if stat_change(rline) && ~stat_change(rline-1)
        cur_regime = cur_regime + 1;
      end
      noise.regime(rline) = cur_regime;
    end
  end
  regimes = unique(noise.regime);
  if length(regimes) > 1
    warning('There are %d noise regimes\n', length(regimes));
  end
  
  %% Apply the filtering from coh_noise_arg across each regime
  old_noise = noise;
  if debug_level > 0
    figure(1); clf;
    imagesc(lp(noise.coh_ave));
    aa = gca;
  end
  
  noise.coh_ave(noise.coh_ave_samples <= param.analysis.coh_ave.min_samples) = NaN;
  
  if debug_level > 0
    figure(2); clf;
    imagesc(lp(noise.coh_ave))
    aa(2) = gca;
  end
  
  mask = isnan(noise.coh_ave);
  mask = filter2(param.analysis.coh_ave.power_grow,double(mask));
  noise.coh_ave(mask > 0) = NaN;
  
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
      noise.coh_ave(regime_mask,regime_fill) = old_noise.coh_ave(regime_fill,regime_mask).';
    end
    for rbin = 1:size(noise.coh_ave,2)
      noise.coh_ave(regime_mask,rbin) = interp_finite(noise.coh_ave(regime_mask,rbin),0);
    end
    if size(noise.coh_ave(regime_mask,:),1) < param.proc.coh_noise_arg{2}+2
      %       sgolayfilt_F = size(noise.coh_ave(regime_mask,:),1);
      %       if mod(sgolayfilt_F,2)==0
      %         sgolayfilt_F = sgolayfilt_F - 1;
      %       end
      %       sgolayfilt_degree = min(param.proc.coh_noise_arg{1}, sgolayfilt_F-1);
      %       noise.coh_ave(regime_mask,:) = single(sgolayfilt(double(noise.coh_ave(regime_mask,:)),sgolayfilt_degree,sgolayfilt_F));
    else
      %    noise.coh_ave(regime_mask,:) = single(sgolayfilt(double(noise.coh_ave(regime_mask,:)),param.proc.coh_noise_arg{1},param.proc.coh_noise_arg{2},param.proc.coh_noise_arg{3}));
      regime_mask_tmp = regime_mask(2:end-1);
      noise.coh_ave(regime_mask_tmp,:) = single(sgolayfilt(double(noise.coh_ave(regime_mask_tmp,:)),param.proc.coh_noise_arg{1},param.proc.coh_noise_arg{2},param.proc.coh_noise_arg{3}));
    end
    if length(regime_mask) >= 3
      noise.coh_ave(regime_mask(1),:) = noise.coh_ave(regime_mask(2),:);
      noise.coh_ave(regime_mask(end),:) = noise.coh_ave(regime_mask(end-1),:);
    end
  end
  
  %% Create the simplified output
  noise_simp = struct('gps_time',noise.gps_time);
  noise_simp.coh_aveI = real(noise.coh_ave);
  noise_simp.coh_aveQ = imag(noise.coh_ave);
  noise_simp.doppler_weights = doppler_weights;
  noise_simp.sw_version = param.sw_version;
  noise_simp.param_collate = param.analysis.coh_ave;
  noise_simp.datestr = datestr(now);
  noise_simp.param_collate.coh_noise_arg = param.proc.coh_noise_arg;
  
  %% Store the simplified output in netcdf file
  out_fn_dir = fileparts(ct_filename_out(param,coh_ave_file_output_type, ''));
  out_fn = fullfile(out_fn_dir,sprintf('coh_noise_simp_%s.nc', param.day_seg));
  fprintf('  Saving %s\n', out_fn);
  netcdf_from_mat(out_fn,noise_simp);
end

return

% Example code to load
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

