% function collate_equalization(param)
%
% Function for estimating equalization coefficients. Can be used for
% transmit and receive equalization.

physical_constants
ref_wf_adc = param.analysis.surf.ref_wf_adc;

%% Load data
fn = fullfile(ct_filename_out(param,input_fn_dir,'',1),sprintf('surf_%s_img_01.mat',param.day_seg));
data = load(fn);
wrap_around_window = hanning(10);
wrap_around_window = [wrap_around_window(6:10); 0];
data.surf_vals(end-5:end,:,:) = bsxfun(@times,data.surf_vals(end-5:end,:,:), ...
  wrap_around_window);
zero_surf_bin = 1-data.param_analysis.analysis.surf.bin_rng(1);
Nt = size(data.surf_vals,1);
Nx = size(data.surf_vals,2);
Nc = size(data.surf_vals,3);

rlines = param.analysis.surf.rlines;
if isempty(rlines)
  rlines = 1:size(data.surf_vals,2);
end

% Only supports a single image right now
img = 1;

if debug_level >= 3
  %% DEBUG
  plot_bins = zero_surf_bin;
  test_wf_adc = min(9,Nc); % <== SET DESIRED CHANNEL TO COMPARE
  figure(1); clf;
  subplot(3,1,1);
  plot(lp(data.surf_vals(plot_bins,:,test_wf_adc) ./ data.surf_vals(plot_bins,:,ref_wf_adc)).','.')
  grid on;
  title(sprintf('Compare wf-adc pair %d to %d',test_wf_adc,ref_wf_adc));
  ylabel('Relative power (dB)');
  h_axis = gca;
  subplot(3,1,2);
  plot(180/pi*angle(data.surf_vals(plot_bins,:,test_wf_adc) .* conj(data.surf_vals(plot_bins,:,ref_wf_adc)) ).','.')
  grid on;
  h_axis(end+1) = gca;
  ylabel('Relative angle (deg)');
  subplot(3,1,3);
  plot(180/pi*data.roll.');
  grid on;
  h_axis(end+1) = gca;
  ylabel('Roll angle (deg)');
  xlabel('Range line');
  figure(2); clf;
  subplot(3,1,1:2);
  imagesc(lp(data.surf_vals(:,:,ref_wf_adc)));
  ylabel('Relative range bin');
  h_axis(end+1) = gca;
  subplot(3,1,3);
  plot(data.surf_bins(:,:).');
  h_axis(end+1) = gca;
  xlabel('Range line');
  ylabel('Surface Bin')
  grid on;
  
  figure(3); clf;
  h_plot = zeros(1,Nc);
  legend_str = cell(1,Nc);
  plot_mode = [0 0 0; hsv(7)];
  plot_bins = zero_surf_bin + (-1:1); % <== SET DESIRED RANGE BIN MULTILOOKING
  Nfir_dec = 5; % <== SET DESIRED ALONG TRACK MULTILOOKING
  ref_rline = min(32055,Nx); % <== SET DESIRED REFERENCE RANGE LINE
  for wf_adc = 1:Nc
    unwrapped_angle = angle(mean(fir_dec(data.surf_vals(plot_bins,:,wf_adc) ...
      .* conj(data.surf_vals(plot_bins,:,ref_wf_adc)),ones(1,Nfir_dec)/Nfir_dec,1)));
    unwrapped_angle = unwrap(unwrapped_angle);
    unwrapped_angle = unwrapped_angle - unwrapped_angle(ref_rline);
    h_plot(wf_adc) = ...
      plot(180/pi* unwrapped_angle, ...
      'Color', plot_mode(mod(wf_adc-1,length(plot_mode))+1,:), ...
      'LineStyle','none','Marker', '.');
    hold on;
    legend_str{wf_adc} = sprintf('wf-adc %d', wf_adc);
  end
  grid on;
  legend(h_plot,legend_str);
  xlabel('Range line');
  ylabel('Relative angle (deg)');
  h_axis(end+1) = gca;
  
  linkaxes(h_axis,'x');
  xlim([1 size(data.surf_vals,2)]);
  keyboard;
end

%% Retrack surface
if param.analysis.surf.retrack.en
  surf_bin = NaN*zeros(1,size(data.surf_vals,2));
  ml_data = lp(fir_dec(abs(data.surf_vals(:,:,ref_wf_adc)).^2,ones(1,5)/5,1));
  for rline = 1:size(data.surf_vals,2)
    % Threshold is hard coded to max of 7 dB of noise or 13 below peak
    cur_threshold = max([ml_data(1,rline)+7; ml_data(:,rline)-13]);
    tmp = find(ml_data(:,rline) > cur_threshold,1);
    if ~isempty(tmp)
      [~,max_offset] = max(ml_data(tmp+(0:2),rline));
      tmp = tmp-1 + max_offset;
      surf_bin(rline) = tmp;
    end
  end
  
  if debug_level >= 2
    figure(1); clf;
    imagesc(ml_data);
    hold on;
    plot(surf_bin);
  end
  
  for rline = 1:size(data.surf_vals,2)
    offset = surf_bin(rline) - zero_surf_bin;
    if offset < 0
      data.surf_vals(1-offset:end,rline,:) = data.surf_vals(1:end+offset,rline,:);
      data.surf_vals(1:1-offset) = 0;
    elseif offset > 0
      data.surf_vals(1:end-offset,rline,:) = data.surf_vals(1+offset:end,rline,:);
      data.surf_vals(end-offset:end) = 0;
    end
  end
  
  ml_data = lp(fir_dec(abs(data.surf_vals(:,:,ref_wf_adc)).^2,ones(1,5)/5,1));
  % Check to make sure surface is flat
  if debug_level >= 2
    figure(2); clf;
    imagesc(ml_data);
    keyboard
  end
end

%% Setup code for estimation equalization coefficients
switch(param.analysis.surf.delay.method)
  case 'threshold'
    delay_method = 1;
  case 'xcorr_complex'
    delay_method = 2;
  case 'xcorr_magnitude'
    delay_method = 3;
  otherwise
    error('delay.method %d is not supported', param.analysis.surf.delay.method);
end

Mt = param.analysis.surf.delay.Mt;
search_bins = param.analysis.surf.delay.bin_rng;
ref_bins = param.analysis.surf.delay.bin_rng;

if delay_method == 2
  Hcorr_wind = boxcar(length(ref_bins));
else
  Hcorr_wind = boxcar(Mt*length(ref_bins));
end
zero_padding_offset = length(search_bins) - length(ref_bins);

% Update create settings and wiki to include discussion of rx equal
% settings which capture the same interface in all waveforms
wf = data.param_analysis.analysis.imgs{img}(1,1);
adc = data.param_analysis.analysis.imgs{img}(1,1);

%% 1. Determine time delay and phase correction for position and channel equalization
dtime = zeros(size(data.elev));
if param.analysis.surf.motion_comp.en
  if isempty(data.param_analysis.get_heights.lever_arm_fh)
    error('No leverarm was used during analysis surf, cannot enable motion_comp');
  end
  drange = bsxfun(@minus,data.elev,mean(data.elev));
  dtime = dtime + drange/(c/2);
end
if param.analysis.surf.chan_eq.en
  Tsys = param.radar.wfs(wf).Tsys - data.param_analysis.radar.wfs(wf).Tsys;
  dtime = bsxfun(@plus, dtime, Tsys.');
  
  Tsys = param.radar.wfs(wf).Tsys;
else
  Tsys = data.param_analysis.radar.wfs(wf).Tsys;
end
dtime = permute(dtime,[3 2 1]);

%% 2. Apply time delay, phase, and amplitude correction
df = 1/(Nt*data.wfs(wf).dt);
freq = data.wfs(wf).fc + df * ifftshift( -floor(Nt/2) : floor((Nt-1)/2) ).';
data.surf_vals = ifft(fft(data.surf_vals) .* exp(1i*2*pi*bsxfun(@times,freq,dtime)));

if param.analysis.surf.chan_eq.en
  % Only apply the relative offset between what has already been applied
  % during analysis surf and the new coefficients
  chan_equal_deg = param.radar.wfs(wf).chan_equal_deg - data.param_analysis.radar.wfs(wf).chan_equal_deg;
  chan_equal_dB = param.radar.wfs(wf).chan_equal_dB - data.param_analysis.radar.wfs(wf).chan_equal_dB;
  data.surf_vals = bsxfun(@times, data.surf_vals, ...
    permute(exp(-1i*chan_equal_deg/180*pi) ./ 10.^(chan_equal_dB/20),[1 3 2]));
  
  chan_equal_deg = param.radar.wfs(wf).chan_equal_deg;
  chan_equal_dB = param.radar.wfs(wf).chan_equal_dB;
else
  chan_equal_deg = data.param_analysis.radar.wfs(wf).chan_equal_deg;
  chan_equal_dB = data.param_analysis.radar.wfs(wf).chan_equal_dB;
end

%% 3. Estimate time delay and amplitude and phase
peak_offset = NaN*zeros(Nc,Nx);
peak_val = NaN*zeros(Nc,Nx);
for rline_idx = 1:length(rlines)
  if debug_level >= 1 && ~mod(rline_idx-1,100)
    fprintf('  Estimating rline index %d of %d\n', rline_idx, length(rlines))
  end
  rline = rlines(rline_idx);
  for wf_adc = 1:Nc
    if delay_method == 2
      %% Cross correlation method with complex data
      in = data.surf_vals(zero_surf_bin+search_bins,rline,wf_adc);
      ref_in = data.surf_vals(zero_surf_bin+search_bins,rline,ref_wf_adc);
      [corr_out,lags] = xcorr(in, ref_in .* Hcorr_wind);
      corr_int = interpft(corr_out,Mt*length(corr_out));
      [peak_val(wf_adc,rline) peak_offset(wf_adc,rline)] = max(corr_int);
      peak_offset(wf_adc,rline) = (peak_offset(wf_adc,rline)-1)/Mt+1 ...
        + ref_bins(1) + search_bins(1) - 1 - zero_padding_offset;
    elseif delay_method == 3
      %% Cross correlation method with magnitude data
      error('Not finished');
      in = interpft(data.surf_vals(zero_surf_bin+search_bins,rline,wf_adc), Mt*length(search_bins));
      ref_in = interpft(data.surf_vals(zero_surf_bin+search_bins,rline,ref_wf_adc), Mt*length(search_bins));
      [corr_int,lags] = xcorr(abs(in), abs(ref_in) .* Hcorr_wind);
    elseif delay_method == 1
      %% Threshold method
      threshold = lp(mean(abs(data.surf_vals(param.analysis.surf.noise_bin,:,ref_wf_adc)).^2)) + 30;
      in = interpft(data.surf_vals(:,rline,wf_adc), Mt*Nt);
      ref_in = interpft(data.surf_vals(:,rline,ref_wf_adc), Mt*Nt);
      
      in_bin = find(lp(in)>threshold,1);
      ref_in_bin = find(lp(ref_in)>threshold,1);
      if ~isempty(in_bin) && ~isempty(ref_in_bin)
        % If the threshold was exceeded, then we use this range line
        peak_offset(wf_adc,rline) = (in_bin - ref_in_bin)/Mt;
        if abs(peak_offset(wf_adc,rline)) > 5
          figure(1); clf;
          plot(lp(in));
          hold on;
          plot(lp(ref_in),'r');
          keyboard
        end
        
        [~,peak_idx] = max(data.surf_vals(zero_surf_bin+search_bins,rline,ref_wf_adc));
        peak_val(wf_adc,rline) = data.surf_vals(zero_surf_bin+search_bins(peak_idx),rline,wf_adc) ...
          ./ data.surf_vals(zero_surf_bin+search_bins(peak_idx),rline,ref_wf_adc);
      end
    end
    
    if 0
      % Debug code
      figure(1); clf;
      plot(lp(in));
      hold on;
      plot(lp(ref_in),'r');
      hold off;
      figure(2); clf;
      plot(lp(corr_int))
      angle(peak_val(wf_adc,rline))
      lp(peak_val(wf_adc,rline))
      peak_offset(wf_adc,rline)
      keyboard
    end
    
  end
end

if delay_method == 2
  peak_offset = bsxfun(@minus,peak_offset,peak_offset(ref_wf_adc,:));
end

equal.peak_offset = peak_offset;
equal.peak_val = peak_val;

equal.Tsys_offset = nanmean(peak_offset(:,param.analysis.surf.rlines),2)*data.wfs(wf).dt;
equal.chan_equal_deg_offset = angle(nanmean(peak_val(:,param.analysis.surf.rlines),2)) * 180/pi;
equal.chan_equal_dB_offset = lp(nanmean(peak_val(:,param.analysis.surf.rlines),2),2);

equal.Tsys_offset_std = nanstd(peak_offset(:,param.analysis.surf.rlines),[],2)*data.wfs(wf).dt;
equal.chan_equal_deg_offset_std = angle(nanstd(peak_val(:,param.analysis.surf.rlines),[],2)) * 180/pi;
equal.chan_equal_dB_offset_std = lp(nanstd(peak_val(:,param.analysis.surf.rlines),[],2),2);

equal.Tsys_offset = reshape(equal.Tsys_offset,[1 Nc]);
equal.chan_equal_deg_offset = reshape(equal.chan_equal_deg_offset,[1 Nc]);
equal.chan_equal_dB_offset = reshape(equal.chan_equal_dB_offset,[1 Nc]);

equal.Tsys_offset_std = reshape(equal.Tsys_offset_std,[1 Nc]);
equal.chan_equal_deg_offset_std = reshape(equal.chan_equal_deg_offset_std,[1 Nc]);
equal.chan_equal_dB_offset_std = reshape(equal.chan_equal_dB_offset_std,[1 Nc]);

equal.Tsys = Tsys + equal.Tsys_offset;
equal.chan_equal_deg = chan_equal_deg + equal.chan_equal_deg_offset;
equal.chan_equal_dB = chan_equal_dB + equal.chan_equal_dB_offset;

equal.Tsys_str = sprintf('[');
equal.Tsys_str = [equal.Tsys_str sprintf('%.1f ',1e9*equal.Tsys(1:end-1))];
equal.Tsys_str = [equal.Tsys_str ...
  sprintf('%.1f]/1e9',1e9*equal.Tsys(end))];

equal.chan_equal_dB_str = sprintf('[');
equal.chan_equal_dB_str = [equal.chan_equal_dB_str sprintf('%.1f ',equal.chan_equal_dB(1:end-1))];
equal.chan_equal_dB_str = [equal.chan_equal_dB_str ...
  sprintf('%.1f]',equal.chan_equal_dB(end))];

equal.chan_equal_deg_str = sprintf('[');
equal.chan_equal_deg_str = [equal.chan_equal_deg_str sprintf('%.1f ',equal.chan_equal_deg(1:end-1))];
equal.chan_equal_deg_str = [equal.chan_equal_deg_str ...
  sprintf('%.1f]',equal.chan_equal_deg(end))];

% If Tsys were inserted, this accounts for the phase shift expected by
% doing this.
equal.chan_equal_deg_with_Tsys = chan_equal_deg + equal.chan_equal_deg_offset + 180*data.wfs(wf).fc*equal.Tsys_offset;

equal.chan_equal_deg_with_Tsys_str = sprintf('[');
equal.chan_equal_deg_with_Tsys_str = [equal.chan_equal_deg_with_Tsys_str sprintf('%.1f ',equal.chan_equal_deg_with_Tsys(1:end-1))];
equal.chan_equal_deg_with_Tsys_str = [equal.chan_equal_deg_with_Tsys_str ...
  sprintf('%.1f]',equal.chan_equal_deg_with_Tsys(end))];

if debug_level >= 1
  fprintf('Offsets from Old Coefficients\n');
  fprintf('%.1f\t', 1e9*equal.Tsys_offset); fprintf('\n');
  fprintf('%.1f\t', equal.chan_equal_deg_offset); fprintf('\n');
  fprintf('%.1f\t', equal.chan_equal_dB_offset); fprintf('\n');
  fprintf('New Coefficients\n');
  fprintf('%.1f\t', 1e9*equal.Tsys); fprintf('\n');
  fprintf('%.1f\t', equal.chan_equal_deg); fprintf('\n');
  fprintf('%.1f\t', equal.chan_equal_dB); fprintf('\n');
  fprintf('New Coefficients if using the old Tsys (in spreadsheet format)\n');
  fprintf('%s\n',equal.chan_equal_dB_str);
  fprintf('%s\n',equal.chan_equal_deg_str);
  fprintf('%s\n',equal.Tsys_str);
  fprintf('New Coefficients if updating this Tsys (in spreadsheet format)\n');
  fprintf('%s\n',equal.chan_equal_dB_str);
  fprintf('%s\n',equal.chan_equal_deg_with_Tsys);
  fprintf('%s\n',equal.Tsys_str);
end

equal_fn_dir = ct_filename_out(param,'equal','',1);
if ~exist(equal_fn_dir)
  mkdir(equal_fn_dir);
end
equal_fn = fullfile(equal_fn_dir,sprintf('Equal_%s.mat',param.day_seg));
equal.param_equal = param;
equal.in_fn = dir(fn);
equal.sw_version = current_software_version;
fprintf('Saving %s\n', equal_fn);
save(equal_fn,'-struct','equal');

return;
