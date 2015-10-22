% script slope_tracker_subtask
%
% This is just a part of slope_tracker_task separated off for easier
% debugging and maintenance.
%
% Author: John Paden

cmd = [];
cmd.create_fcs = 1;
cmd.flatten_elevation = 1;
cmd.check_rx_equalization = 0;
cmd.raw_slope_estimate2D = 0;
cmd.raw_slope_estimate1D = 0;
cmd.SAR_process = 0;

if cmd.create_fcs
  %% Create the reference trajectory
  trajectory_param = struct('gps_source',out_records.gps_source, ...
    'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
    'tx_weights', [], 'lever_arm_fh', param.csarp.lever_arm_fh);
  ref = trajectory_with_leverarm(out_records,trajectory_param);
  along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev,50);
  
  %% Create the flight coordinate system for each channel relative to the reference
  SAR_coord_param.type = 2;
  SAR_coord_param.squint = [0 0 -1].';
  SAR_coord_param.Lsar = 100;
  
  clear fcs y_pos z_pos;
  for wf_adc_idx = 1:size(g_data,3)
    fprintf('Trajectory for %d of %d (%s)\n', wf_adc_idx, size(g_data,3), datestr(now,'HH:MM:SS'));
    
    wf = abs(load_param.load.imgs{1}(wf_adc_idx,1));
    adc = abs(load_param.load.imgs{1}(wf_adc_idx,2));
    trajectory_param = struct('gps_source',out_records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name, ...
      'rx_path', wfs(wf).rx_paths(adc), ...
      'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.csarp.lever_arm_fh);
    
    % Determine the location of this specific phase center
    pc = trajectory_with_leverarm(out_records,trajectory_param);
    
    fcs{wf_adc_idx} = SAR_coord_system(SAR_coord_param,pc,ref,along_track,along_track);
    y_pos(wf_adc_idx,:) = fcs{wf_adc_idx}.pos(2,:);
    z_pos(wf_adc_idx,:) = fcs{wf_adc_idx}.pos(3,:);
  end
  
end

if cmd.flatten_elevation
  Mx = 11;
  
  %% Initialize arrays for quick look image to estimate altimetry
  g_data_coh_ave = [];
  for wf_adc_idx = 1:size(g_data,3)
    g_data_coh_ave(:,:,wf_adc_idx) = fir_dec(g_data(:,:,wf_adc_idx),Mx);
  end
  y_pos_coh_ave = fir_dec(y_pos,Mx);
  z_pos_coh_ave = fir_dec(z_pos,Mx);

  %% Create quick look image
  Nx = size(g_data_coh_ave,2);
  for rline = 1:Nx
    
    [theta array_param.sv] = array_proc_sv(1,new_wfs.fc,y_pos_coh_ave(:,rline),z_pos_coh_ave(:,rline));
    
    g_data_coh_ave(:,rline,1) = mean(g_data_coh_ave(:,rline,:) ...
      .* repmat(permute(array_param.sv,[3 2 1]), [size(g_data_coh_ave,1) 1]), 3);
  end
  g_data_coh_ave = fir_dec(abs(g_data_coh_ave(:,:,1)).^2,3);
  [surf_vals surf_bins] = max(lp(g_data_coh_ave));
  surf_bins = medfilt1(surf_bins,11);
  
  surf_bins = refine_layer(g_data_coh_ave,surf_bins,struct('over_sample',10,'method','spline','bin_rng',-5:5));
  surf_bins = fir_dec(surf_bins,hanning(21)'/sum(hanning(21)));
  
  if 0
    figure(1); clf;
    imagesc(lp(g_data_coh_ave));
    hold on;
    plot(surf_bins);
    hold off;
  end
  
  surf_gps_time = fir_dec(fir_dec(fcs{1}.gps_time,Mx),3);
  
  fcs{1}.surf_bins = interp1(surf_gps_time, surf_bins, fcs{1}.gps_time);
  fcs{1}.surf_bins(isnan(fcs{1}.surf_bins)) = interp1(surf_gps_time, surf_bins, ...
    fcs{1}.gps_time(isnan(fcs{1}.surf_bins)), 'nearest', 'extrap');
  
  [Bfilt,Afilt] = butter(2,[0.001 0.04]);
  elev_hpf = filtfilt(Bfilt,Afilt,out_records.elev);
  
  %% Flatten surface changes by time delaying to a constant elevation
  fprintf('Flatten trajectory (%s)\n', datestr(now,'HH:MM:SS'));
  dt = new_wfs(1).time(2) - new_wfs(1).time(1);
  base_range_corr = c/2*dt * (fcs{1}.surf_bins - min(fcs{1}.surf_bins)) - elev_hpf;
  % To prevent wrap around, apply window roll-off to start of data
  Hwin = hanning(30);
  for rbin = 1:15
    g_data(rbin,:,:) = g_data(rbin,:,:) * Hwin(rbin);
  end
  % Flatten
  for wf_adc_idx = 1:size(g_data,3)
    g_data(:,:,wf_adc_idx) =  ifft(fft( g_data(:,:,wf_adc_idx) ) .* exp(j*2*pi*new_wfs.freq*base_range_corr/(c/2)));
  end
end

if cmd.check_rx_equalization
  % Apply motion compensation
  for wf_adc_idx = 1:size(g_data,3)
    g_data(:,:,wf_adc_idx) =  ifft(fft( g_data(:,:,wf_adc_idx) ) .* exp(j*2*pi*new_wfs.freq*z_pos(wf_adc_idx,:)/(c/2)));
  end

  orig_correction = [-70.4 36.7 0.0 -124.2 79.5 55.2 136.5 -188.5 -159.9 -81.6 99.1 75.1 -180.6 -183.7 20.6];
  for wf_adc_idx = 1:size(g_data,3)
    complex_data = filter2(ones(3,11)/33, g_data(:,:,wf_adc_idx) .* conj(g_data(:,:,3)));
    angle_correction(wf_adc_idx) = angle(sum(complex_data(750:1200))) * 180/pi;
  end
  angle_correction = angle_correction + orig_correction;
  for wf_adc_idx = 1:size(g_data,3)
    fprintf('%.1f ', angle_correction(wf_adc_idx));
  end
  fprintf('\n');
  
  imagesc(hsv_plot(complex_data,-90));
  colormap(hsv(256))
  h_colorbar = colorbar;
  caxis([-pi pi])
  set(get(h_colorbar,'ylabel'),'string','angle (rad)')
  
  return;
  
  % Undo motion compensation
  for wf_adc_idx = 1:size(g_data,3)
    g_data(:,:,wf_adc_idx) =  ifft(fft( g_data(:,:,wf_adc_idx) ) .* exp(-j*2*pi*new_wfs.freq*z_pos(wf_adc_idx,:)/(c/2)));
  end
  
end

if cmd.raw_slope_estimate2D
  % In the slope estimate, consider all cross track channels and
  % range line range
  
  Nsv = 256;
  Nx = 256;
  param.get_heights.B_filter = hanning(Nx);
  param.get_heights.decimate_factor = 10;
  
  rline_rng = -length(param.get_heights.B_filter)/2+1 : length(param.get_heights.B_filter)/2;
  slope_map = [];
  rline_out_idx = 0;
  rlines_out = length(param.get_heights.B_filter)/2:param.get_heights.decimate_factor:size(g_data,2)-length(param.get_heights.B_filter)/2;
  rline_vec = length(param.get_heights.B_filter)/2:param.get_heights.decimate_factor:size(g_data,2)-length(param.get_heights.B_filter)/2;
  fprintf('  rlines %d: processing ', length(rline_vec));
  old_rline_string = '';
  for rline_idx = 1:length(rline_vec)
    rline_string = sprintf('%d', rline_idx);
    fprintf('%s%s', char(8*ones(1,length(old_rline_string))), rline_string);
    old_rline_string = rline_string;
    rline = rline_vec(rline_idx);
    rlines = rline + rline_rng;
    
    % Form kx vector
    kx = fftshift(gen_kx(along_track(rlines), Nx));
    k = 4*pi*new_wfs.fc / c;
    kz = sqrt(k.^2 - kx.^2);
    theta = atan2(kz,kx) - pi/2;
    % plot(theta*180/pi)
    
    % 2D Slope detection
    
    % 1. Create the array's DFT (steering vector, sv) matrix for this section of data
    ypos_ave = mean(y_pos(:,rlines),2);
    zpos_ave = mean(z_pos(:,rlines),2);
    [array_theta array_sv] = array_proc_sv(Nsv,new_wfs.fc,ypos_ave,zpos_ave);
    array_theta = fftshift(array_theta);
    array_sv = fftshift(array_sv,2);
    theta_good = abs(array_theta*180/pi) < 10;
    array_sv = array_sv(:,theta_good);
    array_theta = array_theta(theta_good);

    tmp_fft = fftshift(fft(g_data(:,rlines,:),Nx,2),2);
    tmp_fft = permute(tmp_fft, [2 3 1]);
    tmp_fft_os = zeros(length(theta), length(array_theta), size(tmp_fft,3));
    for bin = 1:size(tmp_fft,3)
      tmp_fft_os(:,:,bin) = tmp_fft(:,:,bin) * array_sv;
    end
    tmp_fft_os = permute(tmp_fft_os, [3 1 2]);

    for array_idx = 1:size(tmp_fft_os,3)
      tmp_fft_os(:,:,array_idx) = wiener2(abs(tmp_fft_os(:,:,array_idx)).^2);
    end
    
    % Select max cross-track angle
    [array_peak_max array_peak_idxs] = max(tmp_fft_os,[],3);
    [array_peak_means] = mean(tmp_fft_os,3);
    
    % Select max along-track angle
    [peak_max along_peak_idx] = max(array_peak_max,[],2);
    array_peak_idx = array_peak_idxs((1:size(array_peak_idxs,1)).' + along_peak_idx*size(array_peak_idxs,1));
    array_peak_mean = array_peak_means((1:size(array_peak_means,1)).' + along_peak_idx*size(array_peak_means,1));
    [peak_mean] = mean(array_peak_max,2);

    % Select max cross-track angle after filtering/culling
    good_mask = lp(peak_max) > -120;
    surface_guard = 100;
    bottom_guard = -300;
    surface = mean(out_records.surface(rlines));
    surface = interp1(new_wfs.time,1:length(new_wfs.time),surface);
    bottom = mean(out_records.bottom(rlines));
    bottom = interp1(new_wfs.time,1:length(new_wfs.time),bottom);
    good_mask(1 : round(surface + surface_guard)) = 0;
    good_mask(round(bottom + bottom_guard) : end) = 0;
    
    array_med_filt_length = 201;
    num_truncate = (array_med_filt_length-1)/2;
    array_peak_idx_good = medfilt1(array_peak_idx(good_mask),array_med_filt_length);
    good_mask(find(good_mask,num_truncate)) = 0;
    good_mask(find(good_mask,num_truncate,'last')) = 0;
    array_peak_idx_good = array_peak_idx_good(1+num_truncate : end-num_truncate);
    
    array_peak_idx_final = interp1(find(good_mask),array_peak_idx_good,1:size(tmp_fft_os,1),'nearest','extrap');
    
    %plot(array_theta(array_peak_idx_final)*180/pi)
    
    
    % Select max cross-track angle after filtering/culling
    good_mask = abs(theta(along_peak_idx)*180/pi) < 5;
    surface_guard = 100;
    bottom_guard = -100;
    surface = mean(out_records.surface(rlines));
    surface = interp1(new_wfs.time,1:length(new_wfs.time),surface);
    bottom = mean(out_records.bottom(rlines));
    bottom = interp1(new_wfs.time,1:length(new_wfs.time),bottom);
    good_mask(1 : round(surface + surface_guard)) = 0;
    good_mask(round(bottom + bottom_guard) : end) = 0;
    %plot(theta(along_peak_idx(good_mask))*180/pi,'x')
    
    [B,A] = butter(2,1/50);
    along_peak_idx_good = filtfilt(B,A,along_peak_idx(good_mask));
    
    along_peak_idx_final = interp1(find(good_mask),along_peak_idx_good,1:size(tmp_fft_os,1),'nearest','extrap');
    %plot(theta(round(along_peak_idx_final))*180/pi)
    
    rline_out_idx = rline_out_idx + 1;
    slope_map_x(:,rline_out_idx) = along_peak_idx_final;
    slope_map_y(:,rline_out_idx) = array_peak_idx_final;
  end
  fprintf('\n');
  slope_map = slope_map_x;
  
  figure(1); clf;
  imagesc(slope_map_x);
  figure(2); clf;
  imagesc(slope_map_y);
  
end


if cmd.raw_slope_estimate1D
  % In the slope estimate, consider all cross track channels and
  % range line range
  
  % Apply motion compensation
  tmp = zeros(size(g_data,1),size(g_data,2));
%   for wf_adc_idx = 1:size(g_data,3)
  for wf_adc_idx = 1:7
    tmp = tmp + ifft(fft( g_data(:,:,wf_adc_idx) ) .* exp(j*2*pi*new_wfs.freq*z_pos(wf_adc_idx,:)/(c/2)));
  end
  
  param.get_heights.B_filter = hanning(250);
  rline_rng = -length(param.get_heights.B_filter)/2+1 : length(param.get_heights.B_filter)/2;
  param.get_heights.decimate_factor = 10;
  slope_map = [];
  rline_out_idx = 0;
  rlines_out = length(param.get_heights.B_filter)/2:param.get_heights.decimate_factor:size(g_data,2)-length(param.get_heights.B_filter)/2;
  rline_vec = length(param.get_heights.B_filter)/2:param.get_heights.decimate_factor:size(g_data,2)-length(param.get_heights.B_filter)/2;
  fprintf('  rlines %d: processing ', length(rline_vec));
  old_rline_string = '';
  for rline_idx = 1:length(rline_vec)
    rline_string = sprintf('%d', rline_idx);
    fprintf('%s%s', char(8*ones(1,length(old_rline_string))), rline_string);
    old_rline_string = rline_string;
    rline = rline_vec(rline_idx);
    rlines = rline + rline_rng;
    
    % Form kx vector
    kx = fftshift(gen_kx(along_track(rlines), Nx));
    k = 4*pi*new_wfs.fc / c;
    kz = sqrt(k.^2 - kx.^2);
    theta = atan2(kz,kx) - pi/2;
    % plot(theta*180/pi)
    
    tmp_fft = fftshift(fft(tmp(:,rlines),Nx,2),2);
    tmp_fft = wiener2(abs(tmp_fft).^2, [3 3]);
    
    %     figure(1); clf;
    %     % imagesc(theta*180/pi,[],lp( tmp_fft ))
    %     imagesc([],[],lp( tmp_fft ))
    %     title(sprintf('rline %d',rline));
    %     hold on;
    
    search_idxs = find(abs(theta)<5/180*pi);
    [val idx] = max(tmp_fft(:,search_idxs),[],2);
    mean_val = mean(abs(tmp_fft(:,search_idxs)).^2,2);
    surface_guard = round(wfs(1).Tpd * wfs(1).fs);
    
    surface = round(interp1(new_wfs.time,1:size(g_data,1),out_records.surface(rline)));
    % surface_multiple_mask = round(interp1(wfs(2).time,1:size(g_data,1),out_records.surface(rline)*2));
    bottom_guard = 300;
    bottom = round(interp1(new_wfs.time,1:size(g_data,1),out_records.bottom(rline)));
    quality = abs(val ./ mean_val);
    good_mask = lp(val,1) > lp(mean_val,1) + 11;
    good_mask(1:surface+surface_guard) = 0;
    good_mask(bottom-bottom_guard:end) = 0;
    
    idx = search_idxs(1) + idx - 1;
    
    if 0
      idx = medfilt1(idx(good_mask),21);
      p = polyfit(find(good_mask),idx,2);
      polyfit_idx = polyval(p,1:length(good_mask));
    else
      % idx = idx(good_mask);
      idx = medfilt1(idx(good_mask),11);
      polyfit_idx = interp1(find(good_mask), idx, 1:length(good_mask),'nearest','extrap');
    end
    
    % %     % plot(theta(idx)*180/pi, 1:size(tmp_fft,1),'k.');
    %     plot(idx, find(good_mask),'k.');
    %     plot(polyfit_idx, 1:length(good_mask),'k-');
    %     hold off;
    % %     xlim([2300 2700]);
    % %     ylim([300 1600]);
    %     keyboard;
    
    rline_out_idx = rline_out_idx + 1;
    slope_map(:,rline_out_idx) = kx(polyfit_idx);
    % Constrain above surface and below bed to be zero slope
    slope_map(1:surface,rline_out_idx) = kx(Nx/2);
    slope_map(bottom:end,rline_out_idx) = kx(Nx/2);
  end
  fprintf('\n');
  
  figure(1); clf;
  imagesc(slope_map);
  
end

if cmd.SAR_process
  
  % Convert slope map to kx
  slope_map_kx = reshape(interp1(1:length(kx), kx, reshape(slope_map,[numel(slope_map) 1])),size(slope_map));
  slope_map_kx = medfilt2(slope_map_kx,[11 11]);
  
  %% Apply slope corrections to phase and magnitude
  fprintf('Flatten slope\n', datestr(now,'HH:MM:SS'));
  slope_map_filt = fir_dec(fir_dec(slope_map_kx,ones(1,11)/11).',ones(1,51)/51).';
  slope_map_filt_linear = interp1(rlines_out,slope_map_filt.',1:size(g_data,2),'linear').';
  slope_map_filt_extrap = interp1(rlines_out,slope_map_filt.',1:size(g_data,2),'neighbor','extrap').';
  slope_map_filt = slope_map_filt_linear;
  slope_map_filt(isnan(slope_map_filt_linear)) = slope_map_filt_extrap(isnan(slope_map_filt_linear));
  
  slope_map_filt = sign(slope_map_filt) .* sqrt( slope_map_filt.^2 ./ (k^2 - slope_map_filt.^2) );
%   slope_map_filt = slope_map_filt ./ k;
  
  ypos = cumsum(1.0*slope_map_filt .* repmat([0 diff(along_track,1,2)],[size(g_data,1) 1]),2);
  
  
  range = new_wfs.time * c/2;
  drange = range(2) - range(1);
  param.filt_len = 10 * (drange);
  param.dx = range(2) - range(1);
  param.method = 'sinc';
  
  if 1
    figure(1); clf;
    imagesc(lp(fir_dec(abs(fir_dec(g_data(:,:,1),12)).^2,1),1))
    colorbar;
    hold on
    plot(1006 - ypos(1006,1:12:end)/drange)
    plot(1544 - ypos(1544,1:12:end)/drange)
    plot(2991 - ypos(2991,1:12:end)/drange)
    plot(3144 - ypos(3144,1:12:end)/drange)
    plot(3557 - ypos(3557,1:12:end)/drange)
    hold off;
  end
  
  %tmp2 = zeros(size(tmp) - [20 0]);
  tmp = zeros(size(g_data,1),size(g_data,2));
  cross_angle = interp1(rline_vec,array_theta(slope_map_y).',rline_vec,'nearest','extrap').';
  g_data = permute(g_data,[1 3 2]);
  old_rline_string = '';
  fprintf('  rlines %d: processing ', size(g_data,3));
  for rline = 1:size(g_data,3)
    rline_string = sprintf('%d', rline);
    fprintf('%s%s', char(8*ones(1,length(old_rline_string))), rline_string);
    old_rline_string = rline_string;
    
    [tmp array_sv] = array_proc_sv({'theta',cross_angle(:,rline)'},new_wfs.fc,y_pos(:,rline),z_pos(:,rline));
    tmp(:,rline) = mean(g_data(:,:,rline) * array_sv,1);
  end
  fprintf('\n');
  g_data = permute(g_data,[1 3 2]);
  
  tmp2 = zeros(size(tmp));
  FUDGE_FACTOR = 1.0; % 1.0 is nominal... I think dependent on FFT size and step
  tmp_phase = tmp .* exp(-j*k*ypos * FUDGE_FACTOR);
  for rline = 1:size(g_data,2)
    actual_range = range + ypos(:,rline);
    tmp2(:,rline) = interp1(actual_range, tmp_phase(:,rline).', range);
    %tmp2(:,rline) = interp1(actual_range, tmp_phase(:,rline).', range(11:end-10));
    %   rline
    %   tmp2(:,rline) = arbitrary_resample(tmp(:,rline).', actual_range, range(11:end-10), param).';
  end
  %tmp2 = tmp2 .* exp(-j*k*ypos(11:end-10,:));
  
  % imagesc(lp(tmp2))
  
  
  %% Apply incoherent averaging with decimation
  %   data_incoh = fir_dec(abs(g_data).^2,param.get_heights.inc_ave);
  
  
  %% Undo elevation correction
  fprintf('Undo elevation correction (%s)\n', datestr(now,'HH:MM:SS'));
  range_corr = ref.elev - ref.elev(rline_ref);
  
  
  %% Undo slope correction
  new_out = fir_dec(abs(fir_dec(tmp2,hanning(41).',6)).^2,hanning(11).',6);
  rlines_out3 = 1:6:size(tmp2,2);
  rlines_out2 = 1:6:length(rlines_out3);
  rlines_out3 = rlines_out3(rlines_out2);
  new_out = interp1(rlines_out3,new_out.',1:size(g_data,2),'linear','extrap').';
  warning off
  tmp3 = zeros(size(tmp));
  for rline = 1:size(g_data,2)
    actual_range = range + ypos(:,rline) - range_corr(rline);
    %   tmp3(:,rline) = interp1(range(11:end-10), new_out(:,rline).', actual_range);
    tmp3(:,rline) = interp1(range, new_out(:,rline).', actual_range);
  end
  warning on;
  
  param.get_heights.decimate_factor = 4;
  
  coh_avgs = param.get_heights.inc_ave;
  siz = size(tmp3);
  new_len     = floor(siz(2)/coh_avgs);
  nearest_len = floor(siz(2)/coh_avgs)*coh_avgs;
  data_incoh = tmp3(:,1+4:10:end-5);
  imagesc(lp(data_incoh));
  
end


return;








if cmd.create_fcs
  
  %% Determine the range correction
  % A positive range correction is when the phase center is above the reference
  % (i.e. higher elevation)
  base_range_corr = out_records.elev - out_records.elev(rline_ref);
  range_corr = zeros(size(g_data,3),size(g_data,2));
  % figure(1); clf;
  % plot(base_range_corr);
  
  trajectory_param = struct('gps_source',out_records.gps_source, ...
    'season_name',param.season_name,'radar_name',param.radar_name,'rx_path', 0, ...
    'tx_weights', [], 'lever_arm_fh', param.csarp.lever_arm_fh);
  ref = trajectory_with_leverarm(out_records,trajectory_param);
  along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev,50);
  
  
  SAR_coord_param.type = 2;
  SAR_coord_param.squint = [0 0 -1].';
  SAR_coord_param.Lsar = 100;
  
  clear fcs y_pos z_pos;
  for wf_adc_idx = 1:size(g_data,3)
    fprintf('Trajectory for %d of %d (%s)\n', wf_adc_idx, size(g_data,3), datestr(now,'HH:MM:SS'));
    
    wf = abs(load_param.load.imgs{1}(wf_adc_idx,1));
    adc = abs(load_param.load.imgs{1}(wf_adc_idx,2));
    trajectory_param = struct('gps_source',out_records.gps_source, ...
      'season_name',param.season_name,'radar_name',param.radar_name, ...
      'rx_path', wfs(wf).rx_paths(adc), ...
      'tx_weights', wfs(wf).tx_weights, 'lever_arm_fh', param.csarp.lever_arm_fh);
    
    % Determine the location of this specific phase center
    pc = trajectory_with_leverarm(out_records,trajectory_param);
    
    fcs{wf_adc_idx} = SAR_coord_system(SAR_coord_param,pc,ref,along_track,along_track);
    y_pos(wf_adc_idx,:) = fcs{wf_adc_idx}.pos(2,:);
    z_pos(wf_adc_idx,:) = fcs{wf_adc_idx}.pos(3,:);
  end
  
  %% Estimate cross track surface slope
  g_data_coh_ave = [];
  for wf_adc_idx = 1:size(g_data,3)
    g_data_coh_ave(:,:,wf_adc_idx) = fir_dec(g_data(:,:,wf_adc_idx),21);
  end
  y_pos_coh_ave = fir_dec(y_pos,21);
  z_pos_coh_ave = fir_dec(z_pos,21);
  
  
  
  
  Nsv = 256;
  rline_out_idx = 0;
  SperiodPre = zeros(Nsv,size(g_data_coh_ave,1));
  slope2_map = [];
  Nx = size(g_data_coh_ave,2);
  
  
  fprintf('  rlines %d: processing ', Nx);
  old_rline_string = '';
  for rline = 1:Nx
    rline_string = sprintf('%d', rline);
    fprintf('%s%s', char(8*ones(1,length(old_rline_string))), rline_string);
    old_rline_string = rline_string;
    
    array_param.rline_rng = -5:5;
    %% Check for edge conditions (not enough data to create dataSample
    % properly)
    if rline+array_param.rline_rng(1) < 1
      rline_rng = 1-rline : array_param.rline_rng(end);
    elseif  rline+array_param.rline_rng(end) > Nx
      rline_rng = array_param.rline_rng(1) : Nx-rline;
    else
      rline_rng = array_param.rline_rng;
    end
    
    [theta array_param.sv] = array_proc_sv(Nsv,new_wfs.fc,y_pos_coh_ave(:,rline),z_pos_coh_ave(:,rline));
    Nc = size(g_data,3);
    
    for bin = 1:size(g_data_coh_ave,1)
      dataSample = g_data_coh_ave(bin,rline+rline_rng,:);
      dataSample = reshape(dataSample,[1*length(rline_rng) Nc]);
      SperiodPre(:,bin) = mean(abs(array_param.sv(:,:)'*dataSample.').^2,2);
    end
    Speriod = fftshift(SperiodPre,1).';
    
    
    
    theta_fftshift = fftshift(theta);
    search_idxs = find(abs(theta_fftshift)<10/180*pi);
    [val idx] = max(Speriod(:,search_idxs),[],2);
    mean_val = median(abs(Speriod(:,search_idxs)).^2,2);
    surface_guard = round(wfs(1).Tpd * wfs(1).fs);
    
    surface = round(interp1(new_wfs.time,1:size(g_data,1),out_records.surface(rline)));
    % surface_multiple_mask = round(interp1(wfs(2).time,1:size(g_data,1),out_records.surface(rline)*2));
    bottom_guard = 300;
    bottom = round(interp1(new_wfs.time,1:size(g_data,1),out_records.bottom(rline)));
    quality = abs(val ./ mean_val);
    good_mask = lp(val,2) > lp(mean_val,1) + 9;
    good_mask(1:surface+surface_guard) = 0;
    good_mask(bottom-bottom_guard:end) = 0;
    
    idx = search_idxs(1) + idx - 1;
    
    if 0
      idx = medfilt1(idx(good_mask),21);
      p = polyfit(find(good_mask),idx,2);
      polyfit_idx = polyval(p,1:length(good_mask));
    else
      % idx = idx(good_mask);
      idx = medfilt1(idx(good_mask),11);
      polyfit_idx = interp1(find(good_mask), idx, 1:length(good_mask),'nearest','extrap');
    end
    
    % %     % plot(theta(idx)*180/pi, 1:size(tmp_fft,1),'k.');
    %     plot(idx, find(good_mask),'k.');
    %     plot(polyfit_idx, 1:length(good_mask),'k-');
    %     hold off;
    % %     xlim([2300 2700]);
    % %     ylim([300 1600]);
    %     keyboard;
    
    rline_out_idx = rline_out_idx + 1;
    slope2_map(:,rline_out_idx) = theta_fftshift(polyfit_idx);
    % Constrain above surface and below bed to be zero slope
    slope2_map(1:surface,rline_out_idx) = theta_fftshift(Nsv/2);
    slope2_map(bottom:end,rline_out_idx) = theta_fftshift(Nsv/2);
    
  end
  fprintf('\n');
  
  slope2_map_filt = medfilt2(slope2_map,[5 15]);
  slope2_map_filt = interp1(fir_dec(1:size(g_data,2),21).',slope2_map_filt.',1:size(g_data,2),'nearest','extrap').';
  
  
  %% Correct for trajectory and attitude changes
  fprintf('Correct for trajectory and attitude changes (%s)\n', datestr(now,'HH:MM:SS'));
  for wf_adc_idx = 1:size(g_data,3)
    g_data(:,:,wf_adc_idx) =  ifft(fft( g_data(:,:,wf_adc_idx) ) .* exp(j*2*pi*new_wfs.freq*range_corr(wf_adc_idx,:)/(c/2)));
  end
  
  
  g_data = permute(g_data,[1 3 2]);
  Nx = size(g_data,3);
  for rline = 1:Nx
    [theta sv_theta] = array_proc_sv({'theta',slope2_map_filt(:,rline)'},new_wfs.fc,y_pos(:,rline),z_pos(:,rline));
    [theta sv_theta_nadir] = array_proc_sv({'theta',0},new_wfs.fc,y_pos(:,rline),z_pos(:,rline));
    
    g_data(:,:,rline) = g_data(:,:,rline) .* conj(sv_theta.') .* repmat(sv_theta_nadir.',[size(g_data,1) 1]);
  end
  g_data = permute(g_data,[1 3 2]);
  
  tmp = mean(g_data,3);
  
end

if 0
  %% Correct for trajectory and attitude changes
  fprintf('Correct for trajectory and attitude changes (%s)\n', datestr(now,'HH:MM:SS'));
  for wf_adc_idx = 1:size(g_data,3)
    g_data(:,:,wf_adc_idx) =  ifft(fft( g_data(:,:,wf_adc_idx) ) .* exp(j*2*pi*new_wfs.freq*range_corr(wf_adc_idx,:)/(c/2)));
  end
  
  %% Correct for cross track slope and other cross track aberrations
  if 0
    fprintf('Correct for cross track slope and other cross track aberations (%s)\n', datestr(now,'HH:MM:SS'));
    for wf_adc_idx = 1:size(g_data,3)
      
      % Form interferogram (couple options)
      %   complex_data = mean(g_data{1}(:,:,3:4),3) .* conj(mean(g_data{1}(:,:,5:6),3));
      complex_data = g_data(:,:,wf_adc_idx) .* conj(g_data(:,:,min(3,size(g_data,3))));
      
      coherence = mean(double(complex_data(1500:1900,:)));
      angle_corr(wf_adc_idx,:) = angle(fir_dec(coherence,hanning(1001).'));
      g_data(:,:,wf_adc_idx) = g_data(:,:,wf_adc_idx).*repmat(exp(-j*angle_corr(wf_adc_idx,:)),[size(g_data,1) 1]);
    end
    
    %% Array processing on cross-track (fft will not work since it is not uniform)
    fprintf('Array processing on cross-track (fft will not work since it is not uniform) (%s)\n', datestr(now,'HH:MM:SS'));
    tmp = mean(g_data,3);
    % figure(1); clf;
    % imagesc(lp(tmp));
  end
  
end


%% FFT in along-track to detect layer slope
fprintf('FFT in along-track to detect layer slope (%s)\n', datestr(now,'HH:MM:SS'));
% For each block determine the slopes and process
Nx = 1240;
param.get_heights.B_filter = hanning(248);
rline_rng = -length(param.get_heights.B_filter)/2+1 : length(param.get_heights.B_filter)/2;
param.get_heights.decimate_factor = 10;
slope_map = [];
rline_out_idx = 0;
rlines_out = length(param.get_heights.B_filter)/2:param.get_heights.decimate_factor:size(g_data,2)-length(param.get_heights.B_filter)/2;
rline_vec = length(param.get_heights.B_filter)/2:param.get_heights.decimate_factor:size(g_data,2)-length(param.get_heights.B_filter)/2;
fprintf('  rlines %d: processing ', length(rline_vec));
old_rline_string = '';
for rline_idx = 1:length(rline_vec)
  rline_string = sprintf('%d', rline_idx);
  fprintf('%s%s', char(8*ones(1,length(old_rline_string))), rline_string);
  old_rline_string = rline_string;
  rline = rline_vec(rline_idx);
  rlines = rline + rline_rng;
  
  % Form kx vector
  kx = fftshift(gen_kx(along_track(rlines), Nx));
  k = 4*pi*new_wfs.fc / c;
  kz = sqrt(k.^2 - kx.^2);
  theta = atan2(kz,kx) - pi/2;
  % plot(theta*180/pi)
  
  tmp_fft = fftshift(fft(tmp(:,rlines),Nx,2),2);
  
  %     figure(1); clf;
  %     % imagesc(theta*180/pi,[],lp( tmp_fft ))
  %     imagesc([],[],lp( tmp_fft ))
  %     title(sprintf('rline %d',rline));
  %     hold on;
  
  search_idxs = find(abs(theta)<5/180*pi);
  [val idx] = max(tmp_fft(:,search_idxs),[],2);
  mean_val = mean(abs(tmp_fft(:,search_idxs)).^2,2);
  surface_guard = round(wfs(1).Tpd * wfs(1).fs);
  
  surface = round(interp1(new_wfs.time,1:size(g_data,1),out_records.surface(rline)));
  % surface_multiple_mask = round(interp1(wfs(2).time,1:size(g_data,1),out_records.surface(rline)*2));
  bottom_guard = 300;
  bottom = round(interp1(new_wfs.time,1:size(g_data,1),out_records.bottom(rline)));
  quality = abs(val ./ mean_val);
  good_mask = lp(val,2) > lp(mean_val,1) + 11;
  good_mask(1:surface+surface_guard) = 0;
  good_mask(bottom-bottom_guard:end) = 0;
  
  idx = search_idxs(1) + idx - 1;
  
  if 0
    idx = medfilt1(idx(good_mask),21);
    p = polyfit(find(good_mask),idx,2);
    polyfit_idx = polyval(p,1:length(good_mask));
  else
    % idx = idx(good_mask);
    idx = medfilt1(idx(good_mask),11);
    polyfit_idx = interp1(find(good_mask), idx, 1:length(good_mask),'nearest','extrap');
  end
  
  % %     % plot(theta(idx)*180/pi, 1:size(tmp_fft,1),'k.');
  %     plot(idx, find(good_mask),'k.');
  %     plot(polyfit_idx, 1:length(good_mask),'k-');
  %     hold off;
  % %     xlim([2300 2700]);
  % %     ylim([300 1600]);
  %     keyboard;
  
  rline_out_idx = rline_out_idx + 1;
  slope_map(:,rline_out_idx) = kx(polyfit_idx);
  % Constrain above surface and below bed to be zero slope
  slope_map(1:surface,rline_out_idx) = kx(Nx/2);
  slope_map(bottom:end,rline_out_idx) = kx(Nx/2);
end
fprintf('\n');

%% Apply slope corrections to phase and magnitude
fprintf('Flatten slope\n', datestr(now,'HH:MM:SS'));
%slope_map_filt = slope_map;
% slope_map_filt = medfilt1(slope_map,11,[],2);
slope_map_filt = fir_dec(fir_dec(slope_map,ones(1,11)/11).',ones(1,51)/51).';
slope_map_filt_linear = interp1(rlines_out,slope_map_filt.',1:size(g_data,2),'linear').';
slope_map_filt_extrap = interp1(rlines_out,slope_map_filt.',1:size(g_data,2),'neighbor','extrap').';
slope_map_filt = slope_map_filt_linear;
slope_map_filt(isnan(slope_map_filt_linear)) = slope_map_filt_extrap(isnan(slope_map_filt_linear));

ypos = cumsum(1.0*slope_map_filt/k .* repmat([0 diff(along_track,1,2)],[size(g_data,1) 1]),2);


range = new_wfs.time * c/2;
drange = range(2) - range(1);
param.filt_len = 10 * (drange);
param.dx = range(2) - range(1);
param.method = 'sinc';

if 0
  figure(1); clf;
  imagesc(lp(fir_dec(abs(fir_dec(tmp,12)).^2,1),1))
  colorbar;
  hold on
  plot(1010 - ypos(1010,1:12:end)/drange)
  plot(1468 - ypos(1468,1:12:end)/drange)
  % plot(2308 - ypos(2308,1:12:end)/drange)
  % plot(2672 - ypos(2672,1:12:end)/drange)
  % plot(2764 - ypos(2764,1:12:end)/drange)
  hold off;
end

%tmp2 = zeros(size(tmp) - [20 0]);
tmp2 = zeros(size(tmp));
FUDGE_FACTOR = 1.2; % 1.0 is nominal... I think dependent on FFT size and step
tmp_phase = tmp .* exp(-j*k*ypos * FUDGE_FACTOR);
for rline = 1:size(g_data,2)
  actual_range = range + ypos(:,rline);
  tmp2(:,rline) = interp1(actual_range, tmp_phase(:,rline).', range);
  %tmp2(:,rline) = interp1(actual_range, tmp_phase(:,rline).', range(11:end-10));
  %   rline
  %   tmp2(:,rline) = arbitrary_resample(tmp(:,rline).', actual_range, range(11:end-10), param).';
end
%tmp2 = tmp2 .* exp(-j*k*ypos(11:end-10,:));

% imagesc(lp(tmp2))


%% Apply incoherent averaging with decimation
%   data_incoh = fir_dec(abs(g_data).^2,param.get_heights.inc_ave);


%% Undo elevation correction
fprintf('Undo elevation correction (%s)\n', datestr(now,'HH:MM:SS'));
range_corr = ref.elev - ref.elev(rline_ref);


%% Undo slope correction
new_out = fir_dec(abs(fir_dec(tmp2,hanning(41).',6)).^2,hanning(11).',6);
rlines_out3 = 1:6:size(tmp2,2);
rlines_out2 = 1:6:length(rlines_out3);
rlines_out3 = rlines_out3(rlines_out2);
new_out = interp1(rlines_out3,new_out.',1:size(g_data,2),'linear','extrap').';
warning off
tmp3 = zeros(size(tmp));
for rline = 1:size(g_data,2)
  actual_range = range + ypos(:,rline) - range_corr(rline);
  %   tmp3(:,rline) = interp1(range(11:end-10), new_out(:,rline).', actual_range);
  tmp3(:,rline) = interp1(range, new_out(:,rline).', actual_range);
end
warning on;

param.get_heights.decimate_factor = 4;

coh_avgs = param.get_heights.inc_ave;
siz = size(tmp3);
new_len     = floor(siz(2)/coh_avgs);
nearest_len = floor(siz(2)/coh_avgs)*coh_avgs;
data_incoh = tmp3(:,1+4:10:end-5);
imagesc(lp(data_incoh));

if 0
  figure(1); clf;
  imagesc(lp(tmp3,1))
  figure(1); caxis([-180 -80]+6)
  % axis([ 70       19931         774        1359]);
  colormap(jet(256));
  
  new_out = filter2(ones(1,11),abs(fir_dec(tmp,hanning(11).',6)).^2);
  rlines_out3 = 1:6:size(tmp,2);
  rlines_out2 = 6:6:size(new_out,2)-5;
  new_out = new_out(:,rlines_out2);
  rlines_out3 = rlines_out3(rlines_out2);
  new_out = interp1(rlines_out3,new_out.',1:size(g_data,2),'linear','extrap').';
  
  figure(2); clf;
  imagesc(lp(new_out,1))
  figure(1); aa=axis; figure(2); axis(aa);
  caxis([-180 -80])
  colormap(jet(256));
end

return





























%% Form interferograms
for wf_adc_idx = 1%:size(g_data,3)
  
  % Form interferogram (couple options)
  %   complex_data = mean(g_data{1}(:,:,3:4),3) .* conj(mean(g_data{1}(:,:,5:6),3));
  complex_data = g_data(:,:,wf_adc_idx) .* conj(g_data(:,:,2));
  %   complex_data = ifft(fft( g_data(:,:,wf_adc_idx) ) .* exp(j*2*pi*wfs(wf).freq*range_corr(wf_adc_idx,:)/(c/2)) ) ...
  %     .* conj( ifft(fft( g_data(:,:,3) ) .* exp(j*2*pi*wfs(wf).freq*range_corr(3,:)/(c/2)) ) );
  
  coherence = mean(double(complex_data(1500:1900,:)));
  angle_corr(wf_adc_idx,:) = angle(fir_dec(coherence,hanning(1001).'));
  complex_data = g_data(:,:,wf_adc_idx).*repmat(exp(-j*angle_corr(wf_adc_idx,:)),[size(g_data,1) 1]) .* conj(g_data(:,:,2));
  
  % Plot interferogram
  figure(1); clf;
  imagesc(hsv_plot(complex_data,[-70 -30]));
  colormap(hsv(256))
  h_colorbar = colorbar;
  caxis([-pi pi])
  set(get(h_colorbar,'ylabel'),'string','angle (rad)')
  
  [vals,idx] = max(complex_data([740:790]*2,:));
  figure(2); clf;
  plot(find(idx == medfilt1(idx,11)),medfilt1(double(angle(vals(idx == medfilt1(idx,11)))*180/pi),101),'r')
  hold on;
  plot(out_records.roll*180/pi)
  hold off;
  title(sprintf('%d', wf_adc_idx));
  %   pause;
  fprintf('%d: %f\n', wf_adc_idx, median(double(angle(vals(idx == medfilt1(idx,11)))*180/pi)));
  %   ylim([-100 100]);
  %   saveas(2,sprintf('~/example2_%d.jpg',wf_adc_idx));
  
end




[ -16.8   75.5    0.0  -96.7  126.8  111.9  170.2 -160.6 -133.8  -52.5  120.7   92.5 -101.1 -137.3   35.8]

for wf_adc_idx = 1:size(g_data,3)
  g_data(:,:,wf_adc_idx) = g_data(:,:,wf_adc_idx) .* exp(j*phase_fix(wf_adc_idx)/180*pi);
end


phase_fix = [57.453823
  40.164520
  0.000000
  33.452213
  45.132740
  53.300823
  33.024754
  40.811426
  27.582607
  37.612646
  28.387407
  0.592820
  67.322769
  31.831495
  16.565794]'
new_phase = load_param.radar.wfs(2).chan_equal_deg - phase_fix
fprintf('['); fprintf('%.1f ',new_phase); fprintf('\b]\n');


beam_width = 16;

500/cosd(beam_width/2) + 2000/cosd(beam_width/2/1.78)*1.78 - 500 - 2000*1.78

% 1% increase in range

% Define FCS and correct for phase

% Remove burst noise

% Pick out good SNR cases
