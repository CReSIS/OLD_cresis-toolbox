% function collate_equal(param, param_override)
%
% Function for estimating equalization coefficients. Can be used for
% transmit and receive equalization.
%
% See "Receiver equalization" wiki page for details.
%
% Use debug_plots 'init', 'final' to create plots for status reports.

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =====================================================================

if ~isfield(param.collate_equal,'chan_eq_en') || isempty(param.collate_equal.chan_eq_en)
  param.collate_equal.chan_eq_en = true;
end

if ~isfield(param.collate_equal,'cmd_idx') || isempty(param.collate_equal.cmd_idx)
  param.collate_equal.cmd_idx = 1;
end
cmd = param.analysis.cmd{param.collate_equal.cmd_idx};

if ~isfield(param.collate_equal,'debug_plots')
  param.collate_equal.debug_plots = {'before_comp','after_comp','surf','final','visible','comp_image'};
end
enable_visible_plot = any(strcmp('visible',param.collate_equal.debug_plots));

if ~isfield(param.collate_equal,'delay') || isempty(param.collate_equal.delay)
  param.collate_equal.delay = [];
end
if ~isfield(param.collate_equal.delay,'method') || isempty(param.collate_equal.delay.method)
  param.collate_equal.delay.method = 'xcorr_complex';
end
if ~isfield(param.collate_equal.delay,'Mt') || isempty(param.collate_equal.delay.Mt)
  param.collate_equal.delay.Mt = 64;
end
if ~isfield(param.collate_equal.delay,'ref_bins') || isempty(param.collate_equal.delay.ref_bins)
  param.collate_equal.delay.ref_bins = -20:20;
end
if ~isfield(param.collate_equal.delay,'search_bins') || isempty(param.collate_equal.delay.search_bins)
  param.collate_equal.delay.search_bins = -7:7;
end

if ~isfield(param.collate_equal,'img_lists') || isempty(param.collate_equal.img_lists)
  % Default is to do all of the images and to do each image separately
  % (i.e. no image to image equalization).
  param.collate_equal.img_lists = num2cell(1:length(param.analysis.imgs));
end

if ~isfield(param.collate_equal,'in_dir') || isempty(param.collate_equal.in_dir)
  param.collate_equal.in_dir = 'analysis';
end

if ~isfield(param.collate_equal,'motion_comp_en') || isempty(param.collate_equal.motion_comp_en)
  param.collate_equal.motion_comp_en = true;
end

if ~isfield(param.collate_equal,'out_dir') || isempty(param.collate_equal.out_dir)
  param.collate_equal.out_dir = 'analysis';
end

if ~isfield(param.collate_equal,'ref') || isempty(param.collate_equal.ref)
  param.collate_equal.ref = 1;
end

if ~isfield(param.collate_equal,'retrack_en') || isempty(param.collate_equal.retrack_en)
  param.collate_equal.retrack_en = true;
end

if ~isfield(param.collate_equal,'rlines') || isempty(param.collate_equal.rlines)
  param.collate_equal.rlines = [];
end

if ~isfield(param.collate_equal,'wf_adcs') || isempty(param.collate_equal.wf_adcs)
  param.collate_equal.wf_adcs = [];
end

if ~isfield(param.collate_equal,'zero_surf_bin') || isempty(param.collate_equal.zero_surf_bin)
  param.collate_equal.zero_surf_bin = [];
end

% Other Setup
% =========================================================================
physical_constants;

if ~isempty(param.collate_equal.debug_plots)
  h_fig = get_figures(3,enable_visible_plot);
end

for img_lists_idx = 1:length(param.collate_equal.img_lists)
  %% Load waveform data
  % =====================================================================
  img = param.collate_equal.img_lists{img_lists_idx}(1);
  wf_adc_list = [];
  for sub_img_idx = 1:length(param.collate_equal.img_lists{img_lists_idx})
    sub_img = param.collate_equal.img_lists{img_lists_idx}(sub_img_idx);
    
    if isempty(param.collate_equal.wf_adcs)
      wf_adcs = 1:size(param.analysis.imgs{sub_img},1);
    else
      wf_adcs = param.collate_equal.wf_adcs{img_lists_idx}{sub_img_idx};
    end
    for wf_adc = wf_adcs
      wf = param.analysis.imgs{sub_img}(wf_adc,1);
      adc = param.analysis.imgs{sub_img}(wf_adc,2);
      wf_adc_list = [wf_adc_list; [wf adc]];
      
      % Load the waveform file
      fn_dir = fileparts(ct_filename_out(param,param.collate_equal.in_dir));
      fn = fullfile(fn_dir,sprintf('waveform_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
      fprintf('Loading %s (%s)\n', fn, datestr(now));
      waveform = load(fn);
      if wf_adc == wf_adcs(1)
        gps_time = waveform.gps_time;
        lat = waveform.lat;
        lon = waveform.lon;
        elev = waveform.elev;
        roll = waveform.roll;
        pitch = waveform.pitch;
        heading = waveform.heading;
        time_rng = waveform.time_rng;
        wf_data = waveform.wf_data;
      else
        gps_time(end+1,:) = waveform.gps_time;
        lat(end+1,:) = waveform.lat;
        lon(end+1,:) = waveform.lon;
        elev(end+1,:) = waveform.elev;
        roll(end+1,:) = waveform.roll;
        pitch(end+1,:) = waveform.pitch;
        heading(end+1,:) = waveform.heading;
        time_rng(:,:,end+1) = waveform.time_rng;
        wf_data(:,:,end+1) = waveform.wf_data;
      end
    end
  end
  dt = waveform.dt;
  fc = waveform.fc;
  ref_wf_adc_idx = param.collate_equal.ref;
  
  % Taper off end of record to reduce circular convolution effects that may
  % show up during time delay compensation.
  wrap_around_window = hanning(10);
  wrap_around_window = [wrap_around_window(6:10); 0];
  wf_data(end-5:end,:,:) = bsxfun(@times,wf_data(end-5:end,:,:), ...
    wrap_around_window);
  
  % Determine the zero_surf_binwhere the surface (or desired calibration
  % target) return should be.
  if ~isempty(param.collate_equal.zero_surf_bin);
    % Hard coded zero surf bin
    zero_surf_bin = param.collate_equal.zero_surf_bin;
  else
    if isnumeric(cmd.start_time)
      error('You must specify param.collate_equal.zero_surf_bin because a hard coded cmd.start_time was used.');
    elseif isstruct(cmd.start_time)
      s = 0;
      eval(cmd.start_time.eval.cmd);
    elseif ischar(cmd.start_time)
      s = 0;
      eval(cmd.start_time);
    end
    zero_surf_bin = round(1-s/dt);
  end
  
  % Dimensions of waveform data
  Nt = size(wf_data,1);
  Nx = size(wf_data,2);
  Nc = size(wf_data,3);
  
  % Approximate time axis
  time = waveform.time_rng(1) + dt*(0:Nt-1);
  
  % Determine which range lines will be used for processing
  if isempty(param.collate_equal.rlines)
    % If user passes in empty rlines, then do all rlines
    rlines = 1:size(wf_data,2);
  else
    % Do rlines that user specified; throw out invalid lines using intersect
    rlines = intersect(param.collate_equal.rlines,1:size(wf_data,2));
  end
  
  %% Plot phase/roll of raw data
  % =====================================================================
  if any(strcmp('before_comp',param.collate_equal.debug_plots))
    param.collate_equal.debug_wf_adc_idx = ref_wf_adc_idx-1; % <== SET DESIRED CHANNEL TO TEST
    debug_wf_adc_idx = min(max(1,param.collate_equal.debug_wf_adc_idx),Nc);
    plot_bins = zero_surf_bin;
    
    clf(h_fig(1));
    pos = get(h_fig(1),'Position');
    set(h_fig(1),'Position',[pos(1:2) 700 800]);
    h_axes = subplot(3,1,1,'parent',h_fig(1));
    plot(h_axes(1), lp(wf_data(plot_bins,:,debug_wf_adc_idx) ./ wf_data(plot_bins,:,ref_wf_adc_idx)).','.')
    grid(h_axes(1),'on');
    title(h_axes(1),sprintf('Compare wf-adc pair %d to %d',debug_wf_adc_idx,ref_wf_adc_idx));
    ylabel(h_axes(1),'Relative power (dB)');
    h_axes(2) = subplot(3,1,2,'parent',h_fig(1));
    plot(h_axes(2), 180/pi*angle(wf_data(plot_bins,:,debug_wf_adc_idx) .* conj(wf_data(plot_bins,:,ref_wf_adc_idx)) ).','.')
    grid(h_axes(2),'on');
    ylabel(h_axes(2),'Relative angle (deg)');
    h_axes(3) = subplot(3,1,3,'parent',h_fig(1));
    plot(h_axes(3),180/pi*roll.');
    grid(h_axes(3),'on');
    ylabel(h_axes(3),'Roll angle (deg)');
    xlabel(h_axes(3),'Range line');
    
    clf(h_fig(2));
    h_axes(4) = subplot(3,1,1:2,'parent',h_fig(2));
    imagesc(lp(wf_data(:,:,ref_wf_adc_idx)),'parent',h_axes(4));
    ylabel(h_axes(4), 'Relative range bin');
    h_axes(5) = subplot(3,1,3,'parent',h_fig(2));
    plot(h_axes(5), time_rng(:,:).'*3e8/2);
    xlabel(h_axes(5), 'Range line');
    ylabel(h_axes(5), 'Surface range (m)')
    grid(h_axes(5), 'on');
    
    clf(h_fig(3));
    h_axes(6) = axes('parent',h_fig(3));
    h_plot = zeros(1,Nc);
    legend_str = cell(1,Nc);
    plot_mode = [0 0 0; hsv(7)];
    plot_bins = zero_surf_bin + (-1:1); % <== SET DESIRED RANGE BIN MULTILOOKING
    Nfir_dec = 11; % <== SET DESIRED ALONG TRACK MULTILOOKING
    total_coherence = zeros(1,Nx);
    for wf_adc = 1:Nc
      interferogram = mean(fir_dec(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx)) ./ abs(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx))),ones(1,Nfir_dec)/Nfir_dec,1));
      total_coherence = total_coherence + abs(interferogram);
    end
    total_coherence = ct_smooth(total_coherence,0.01);
    [~,max_coherence_rline] = max(total_coherence);
    ref_rline = min(max_coherence_rline,Nx); % <== SET DESIRED REFERENCE RANGE LINE
    for wf_adc = 1:Nc
      interferogram = mean(fir_dec(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx)) ./ abs(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx))),ones(1,Nfir_dec)/Nfir_dec,1));
      unwrapped_angle = angle(interferogram);
      unwrapped_angle = unwrap(unwrapped_angle);
      unwrapped_angle = unwrapped_angle - unwrapped_angle(ref_rline);
      h_plot(wf_adc) = ...
        plot(h_axes(6), 180/pi* unwrapped_angle, ...
        'Color', plot_mode(mod(wf_adc-1,length(plot_mode))+1,:), ...
        'LineStyle','none','Marker', '.');
      hold(h_axes(6), 'on');
      legend_str{wf_adc} = sprintf('%d-%d', wf_adc_list(wf_adc,1), wf_adc_list(wf_adc,2));
    end
    grid(h_axes(6), 'on');
    legend(h_plot,legend_str);
    xlabel(h_axes(6), 'Range line');
    ylabel(h_axes(6), 'Relative angle (deg)');
    title(h_axes(6), sprintf('Relative phase between channels, ref range line: %d', ref_rline))
    
    linkaxes(h_axes,'x');
    xlim(h_axes(1), [1 size(wf_data,2)]);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_before_single_img_%02d',param.collate_equal.out_dir,img)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    fig_fn_dir = fileparts(fig_fn);
    if ~exist(fig_fn_dir,'dir')
      mkdir(fig_fn_dir);
    end
    saveas(h_fig(1),fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_before_single_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(1),fig_fn);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_before_surface_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(2),fig_fn);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_before_all_img_%02d',param.collate_equal.out_dir,img)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(3),fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_before_all_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(3),fig_fn);
    
    set(h_fig(1),'Position',pos);
    
    if enable_visible_plot
      keyboard
    end
  end
  
  %% Retrack surface
  % =====================================================================
  if param.collate_equal.retrack_en
    ml_data = fir_dec(abs(wf_data(:,:,ref_wf_adc_idx)).^2,ones(1,5)/5,1);
    
    surf_param = param;
    surf_param.cmd.frms = 1;
    surf_param.qlook.surf.min_bin = time(1);
    surf_param.qlook.surf.threshold_noise_rng = [0 (time(1)-time(zero_surf_bin))*2/3 (time(1)-time(zero_surf_bin))*1/3];
    surf_param.qlook.surf.threshold_rel_max = -9;
    surf_param.qlook.surf.max_rng = [0 0];
    surf_param.qlook.surf.en = true;
    surf_param.layer_tracker.echogram_source = struct('Data',ml_data,'Time',time,'GPS_time',gps_time,'Latitude',lat,'Longitude',lon,'Elevation',elev);
    surf_bin = layer_tracker(surf_param,[]);
    surf_bin = round(interp1(time,1:length(time),surf_bin));
    surf_bin = surf_bin + 1;
    
    for rline = 1:size(wf_data,2)
      offset = surf_bin(rline) - zero_surf_bin;
      if offset < 0
        wf_data(1-offset:end,rline,:) = wf_data(1:end+offset,rline,:);
        wf_data(1:1-offset) = 0;
      elseif offset > 0
        wf_data(1:end-offset,rline,:) = wf_data(1+offset:end,rline,:);
        wf_data(end-offset:end) = 0;
      end
    end
  end
  
  if any(strcmp('surf',param.collate_equal.debug_plots))
    clf(h_fig(1));
    h_axes = axes('parent',h_fig(1));
    imagesc(lp(ml_data),'parent',h_axes);
    colormap(h_axes,1-gray(256));
    hold(h_axes,'on');
    plot(h_axes, surf_bin);
    
    if param.collate_equal.retrack_en
      % Check to make sure surface is flat
      clf(h_fig(2));
      h_axes(2) = axes('parent',h_fig(2));
      imagesc(lp(fir_dec(abs(wf_data(:,:,ref_wf_adc_idx)).^2,ones(1,5)/5,1)),'parent',h_axes(2));
      colormap(h_axes(2),1-gray(256));
      
      linkaxes(h_axes);
    end
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_track1_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    fig_fn_dir = fileparts(fig_fn);
    if ~exist(fig_fn_dir,'dir')
      mkdir(fig_fn_dir);
    end
    saveas(h_fig(1),fig_fn);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_track2_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(2),fig_fn);
    
    if enable_visible_plot
      keyboard
    end
  end
  
  %% Setup for estimating equalization coefficients
  switch(param.collate_equal.delay.method)
    case 'threshold'
      delay_method = 1;
    case 'xcorr_complex'
      delay_method = 2;
    case 'xcorr_magnitude'
      delay_method = 3;
    otherwise
      error('delay.method %d is not supported', param.collate_equal.delay.method);
  end
  
  Mt = param.collate_equal.delay.Mt;
  if isfield(param.collate_equal.delay,'bin_rng')
    warning('Deprecated settings bin_rng used. Use search_bins and ref_bins instead. Support will be removed at some point.');
    param.collate_equal.delay.search_bins = param.collate_equal.delay.bin_rng;
    param.collate_equal.delay.ref_bins = param.collate_equal.delay.bin_rng;
  end
  search_bins = param.collate_equal.delay.search_bins(1) : param.collate_equal.delay.search_bins(end);
  ref_bins = param.collate_equal.delay.ref_bins(1) : param.collate_equal.delay.ref_bins(end);
  
  if delay_method == 2
    Hcorr_wind = boxcar(length(ref_bins));
  else
    Hcorr_wind = boxcar(Mt*length(ref_bins));
  end
  zero_padding_offset = length(search_bins) - length(ref_bins);
  
  %% Compensate data for motion/channel equalization
  % =========================================================================
  
  %% Compensate: 1. Create wf-adc to rx map
  rx_paths = zeros(1,size(wf_adc_list,1));
  for wf_adc = 1:Nc
    wf = wf_adc_list(wf_adc,1);
    adc = wf_adc_list(wf_adc,2);
    rx_paths(wf_adc) = param.radar.wfs(wf).rx_paths(adc);
  end
  
  %% Compensate: 2. Determine motion/channel compensation
  % Determine time delay and phase correction for position and channel equalization
  dtime = zeros(size(elev));
  if param.collate_equal.motion_comp_en
    if isempty(waveform.param_analysis.radar.lever_arm_fh)
      error('No leverarm was used during analysis waveform, cannot enable motion_comp');
    end
    drange = bsxfun(@minus,elev,mean(elev));
    dtime = dtime + drange/(c/2);
  end
  if param.collate_equal.chan_eq_en
    if isfield(param.radar.wfs(wf),'Tadc_adjust')
      new_Tadc_adjust = param.radar.wfs(wf).Tadc_adjust;
    else
      new_Tadc_adjust = 0;
    end
    if isfield(waveform.param_analysis.radar.wfs(wf),'Tadc_adjust')
      old_Tadc_adjust = waveform.param_analysis.radar.wfs(wf).Tadc_adjust;
    else
      old_Tadc_adjust = 0;
    end
    Tsys = param.radar.wfs(wf).Tsys(rx_paths) - waveform.param_analysis.radar.wfs(wf).Tsys(rx_paths) ...
      - (new_Tadc_adjust - old_Tadc_adjust);
    dtime = bsxfun(@plus, dtime, Tsys.');
    
    Tsys = param.radar.wfs(wf).Tsys(rx_paths);
  else
    Tsys = waveform.param_analysis.radar.wfs(wf).Tsys(rx_paths);
  end
  dtime = permute(dtime,[3 2 1]);
  
  %% Compensate: 3. Apply motion/channel compensation
  df = 1/(Nt*dt);
  freq = fc + df * ifftshift( -floor(Nt/2) : floor((Nt-1)/2) ).';
  wf_data = ifft(fft(wf_data) .* exp(1i*2*pi*bsxfun(@times,freq,dtime)));
  
  if param.collate_equal.chan_eq_en
    % Only apply the relative offset between what has already been applied
    % during analysis surf and the new coefficients
    chan_equal_deg = param.radar.wfs(wf).chan_equal_deg(rx_paths) ...
      - waveform.param_analysis.radar.wfs(wf).chan_equal_deg(rx_paths);
    chan_equal_dB = param.radar.wfs(wf).chan_equal_dB(rx_paths) ...
      - waveform.param_analysis.radar.wfs(wf).chan_equal_dB(rx_paths);
    wf_data = bsxfun(@times, wf_data, ...
      permute(exp(-1i*chan_equal_deg/180*pi) ./ 10.^(chan_equal_dB/20),[1 3 2]));
    
    chan_equal_deg = param.radar.wfs(wf).chan_equal_deg(rx_paths);
    chan_equal_dB = param.radar.wfs(wf).chan_equal_dB(rx_paths);
  else
    chan_equal_deg = waveform.param_analysis.radar.wfs(wf).chan_equal_deg(rx_paths);
    chan_equal_dB = waveform.param_analysis.radar.wfs(wf).chan_equal_dB(rx_paths);
  end
  
  if any(strcmp('comp_image',param.collate_equal.debug_plots))
    h_comp_fig = get_figures(Nc,true,'comp_image');
    set(h_comp_fig,'WindowStyle','docked')
    for wf_adc = 1:Nc
      clf(h_comp_fig(wf_adc));
      wf = wf_adc_list(wf_adc,1);
      adc = wf_adc_list(wf_adc,2);
      h_axes(wf_adc) = axes('parent',h_comp_fig(wf_adc));
      imagesc(lp(wf_data(:,:,wf_adc)),'Parent',h_axes(end), 'parent', h_axes(wf_adc));
      title(h_axes(wf_adc),sprintf('wf %d adc %d', wf, adc));
      xlabel(h_axes(wf_adc),'Range line');
      ylabel(h_axes(wf_adc),'Range bin');
    end
    linkaxes(h_axes,'xy');
  end
  
  %% Plot phase/roll after compensation
  % =========================================================================
  if any(strcmp('after_comp',param.collate_equal.debug_plots))
    param.collate_equal.debug_wf_adc_idx = ref_wf_adc_idx-1; % <== SET DESIRED CHANNEL TO TEST
    debug_wf_adc_idx = min(max(1,param.collate_equal.debug_wf_adc_idx),Nc);
    plot_bins = zero_surf_bin;
    
    clf(h_fig(1));
    pos = get(h_fig(1),'Position');
    set(h_fig(1),'Position',[pos(1:2) 700 800]);
    h_axes = subplot(3,1,1,'parent',h_fig(1));
    plot(h_axes(1), lp(wf_data(plot_bins,:,debug_wf_adc_idx) ./ wf_data(plot_bins,:,ref_wf_adc_idx)).','.')
    grid(h_axes(1),'on');
    title(h_axes(1),sprintf('Compare wf-adc pair %d to %d',debug_wf_adc_idx,ref_wf_adc_idx));
    ylabel(h_axes(1),'Relative power (dB)');
    h_axes(2) = subplot(3,1,2,'parent',h_fig(1));
    plot(h_axes(2), 180/pi*angle(wf_data(plot_bins,:,debug_wf_adc_idx) .* conj(wf_data(plot_bins,:,ref_wf_adc_idx)) ).','.')
    grid(h_axes(2),'on');
    ylabel(h_axes(2),'Relative angle (deg)');
    h_axes(3) = subplot(3,1,3,'parent',h_fig(1));
    plot(h_axes(3),180/pi*roll.');
    grid(h_axes(3),'on');
    ylabel(h_axes(3),'Roll angle (deg)');
    xlabel(h_axes(3),'Range line');
    
    clf(h_fig(2));
    h_axes(4) = subplot(3,1,1:2,'parent',h_fig(2));
    imagesc(lp(wf_data(:,:,ref_wf_adc_idx)),'parent',h_axes(4));
    ylabel(h_axes(4), 'Relative range bin');
    h_axes(5) = subplot(3,1,3,'parent',h_fig(2));
    plot(h_axes(5), time_rng(:,:).'*3e8/2);
    xlabel(h_axes(5), 'Range line');
    ylabel(h_axes(5), 'Surface range (m)')
    grid(h_axes(5), 'on');
    
    clf(h_fig(3));
    h_axes(6) = subplot(3,1,1:2,'parent',h_fig(3));
    h_plot = zeros(1,Nc);
    legend_str = cell(1,Nc);
    plot_mode = [0 0 0; hsv(7)];
    plot_bins = zero_surf_bin + (-1:1); % <== SET DESIRED RANGE BIN MULTILOOKING
    Nfir_dec = 11; % <== SET DESIRED ALONG TRACK MULTILOOKING
    total_coherence = zeros(1,Nx);
    for wf_adc = 1:Nc
      interferogram = mean(fir_dec(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx)) ./ abs(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx))),ones(1,Nfir_dec)/Nfir_dec,1));
      total_coherence = total_coherence + abs(interferogram);
    end
    total_coherence = ct_smooth(total_coherence,0.01);
    [~,max_coherence_rline] = max(total_coherence);
    ref_rline = min(max_coherence_rline,Nx); % <== SET DESIRED REFERENCE RANGE LINE
    for wf_adc = 1:Nc
      interferogram = mean(fir_dec(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx)) ./ abs(wf_data(plot_bins,:,wf_adc) ...
        .* conj(wf_data(plot_bins,:,ref_wf_adc_idx))),ones(1,Nfir_dec)/Nfir_dec,1));
      unwrapped_angle = angle(interferogram);
      unwrapped_angle = unwrap(unwrapped_angle);
      unwrapped_angle = unwrapped_angle - unwrapped_angle(ref_rline);
      h_plot(wf_adc) = ...
        plot(h_axes(6), 180/pi* unwrapped_angle, ...
        'Color', plot_mode(mod(wf_adc-1,length(plot_mode))+1,:), ...
        'LineStyle','none','Marker', '.');
      hold(h_axes(6), 'on');
      legend_str{wf_adc} = sprintf('%d-%d', wf_adc_list(wf_adc,1), wf_adc_list(wf_adc,2));
    end
    grid(h_axes(6), 'on');
    legend(h_plot,legend_str);
    ylabel(h_axes(6), 'Relative angle (deg)');
    title(h_axes(6), sprintf('Relative phase between channels, ref range line: %d', ref_rline))
    
    h_axes(7) = subplot(3,1,3,'parent',h_fig(3));
    plot(h_axes(7), 180/pi*roll.');
    grid(h_axes(7), 'on');
    ylabel(h_axes(7), 'Roll angle (deg)');
    xlabel(h_axes(7), 'Range line');
    
    linkaxes(h_axes,'x');
    if ~isempty(param.collate_equal.rlines)
      xlim(h_axes(1), [min(param.collate_equal.rlines) max(param.collate_equal.rlines)]);
    else
      xlim(h_axes(1), [1 size(wf_data,2)]);
    end
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_after_single_img_%02d',param.collate_equal.out_dir,img)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    fig_fn_dir = fileparts(fig_fn);
    if ~exist(fig_fn_dir,'dir')
      mkdir(fig_fn_dir);
    end
    saveas(h_fig(1),fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_after_single_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(1),fig_fn);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_after_surface_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(2),fig_fn);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_after_all_img_%02d',param.collate_equal.out_dir,img)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(3),fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_after_all_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(3),fig_fn);
    
    set(h_fig(1),'Position',pos);
    
    if enable_visible_plot
      keyboard
    end
  end
  
  %% Estimate time delay and amplitude and phase
  % =========================================================================
  peak_offset = NaN*zeros(Nc,Nx);
  peak_val = NaN*zeros(Nc,Nx);
  for rline_idx = 1:length(rlines)
    rline = rlines(rline_idx);
    if ~mod(rline_idx-1,10^floor(log10(length(rlines))-1))
      fprintf('  Estimating rline index %d of %d (rline %d of %d) (%s)\n', rline_idx, length(rlines), rline, Nx, datestr(now));
    end
    for wf_adc = 1:Nc
      if delay_method == 2
        % Time delay: cross correlation method with complex data
        in = wf_data(zero_surf_bin+search_bins,rline,wf_adc);
        ref_in = wf_data(zero_surf_bin+ref_bins,rline,ref_wf_adc_idx);
        [corr_out,lags] = xcorr(in, ref_in .* Hcorr_wind);
        corr_int = interpft(corr_out,Mt*length(corr_out));
        [peak_val(wf_adc,rline) peak_offset(wf_adc,rline)] = max(corr_int);
        peak_offset(wf_adc,rline) = (peak_offset(wf_adc,rline)-1)/Mt+1 ...
          + ref_bins(1) + search_bins(1) - 1 - zero_padding_offset;
        peak_val(wf_adc,rline) = abs( ...
          max(wf_data(zero_surf_bin+search_bins,rline,wf_adc)) ...
          ./ max(wf_data(zero_surf_bin+search_bins,rline,ref_wf_adc_idx))) ...
          .* exp(1i*angle(peak_val(wf_adc,rline)));
      elseif delay_method == 3
        % Time delay: cross correlation method with magnitude data
        error('Not finished');
        in = interpft(wf_data(zero_surf_bin+search_bins,rline,wf_adc), Mt*length(search_bins));
        ref_in = interpft(wf_data(zero_surf_bin+ref_bins,rline,ref_wf_adc_idx), Mt*length(search_bins));
        [corr_int,lags] = xcorr(abs(in), abs(ref_in) .* Hcorr_wind);
      elseif delay_method == 1
        % Time delay: threshold method
        threshold = lp(mean(abs(wf_data(param.collate_equal.noise_bin,:,ref_wf_adc_idx)).^2)) + 30;
        in = interpft(wf_data(:,rline,wf_adc), Mt*Nt);
        ref_in = interpft(wf_data(:,rline,ref_wf_adc_idx), Mt*Nt);
        
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
          
          [~,peak_idx] = max(wf_data(zero_surf_bin+search_bins,rline,ref_wf_adc_idx));
          peak_val(wf_adc,rline) = wf_data(zero_surf_bin+search_bins(peak_idx),rline,wf_adc) ...
            ./ wf_data(zero_surf_bin+search_bins(peak_idx),rline,ref_wf_adc_idx);
        end
      end
      
    end
  end
  
  if delay_method == 2
    peak_offset = bsxfun(@minus,peak_offset,peak_offset(ref_wf_adc_idx,:));
    % peak_val for the reference channel should be 1, so this line of code
    % should have no affect.
    peak_val = bsxfun(@times,peak_val,1./peak_val(ref_wf_adc_idx,:));
  end
  
  %% Plot final results
  % =========================================================================
  if any(strcmp('final',param.collate_equal.debug_plots))
    
    clf(h_fig(1));
    set(h_fig(1),'Name','Final Power');
    h_axes = axes('parent',h_fig(1));
    h_plot = zeros(1,Nc);
    legend_str = cell(1,Nc);
    plot_mode = [0 0 0; hsv(7)];
    for wf_adc = 1:Nc
      h_plot(wf_adc) = plot(h_axes,lp(ct_smooth(peak_val(wf_adc,:),0.01)), ...
        'Color', plot_mode(mod(wf_adc-1,length(plot_mode))+1,:), ...
        'LineStyle','none','Marker', '.');
      hold(h_axes, 'on');
      legend_str{wf_adc} = sprintf('%d-%d', wf_adc_list(wf_adc,1), wf_adc_list(wf_adc,2));
    end
    title(h_axes, 'Comparison of amplitudes relative to ref');
    grid(h_axes, 'on');
    legend(h_plot,legend_str,'location','NorthEastOutside');
    ylabel(h_axes, 'Relative power (dB)');
    xlabel(h_axes, 'Range line');
    
    clf(h_fig(2));
    set(h_fig(2),'Name','Final Phase');
    h_axes(2) = axes('parent',h_fig(2));
    h_plot = zeros(1,Nc);
    legend_str = cell(1,Nc);
    plot_mode = [0 0 0; hsv(7)];
    for wf_adc = 1:Nc
      h_plot(wf_adc) = plot(h_axes(2),angle(ct_smooth(peak_val(wf_adc,:),0.01)), ...
        'Color', plot_mode(mod(wf_adc-1,length(plot_mode))+1,:), ...
        'LineStyle','none','Marker', '.');
      hold(h_axes(2), 'on');
      legend_str{wf_adc} = sprintf('%d-%d', wf_adc_list(wf_adc,1), wf_adc_list(wf_adc,2));
    end
    title(h_axes(2), 'Comparison of phases relative to ref');
    grid(h_axes(2), 'on');
    legend(h_plot,legend_str,'location','NorthEastOutside');
    ylabel(h_axes(2), 'Relative angle (deg)');
    xlabel(h_axes(2), 'Range line');
    
    clf(h_fig(3));
    set(h_fig(3),'Name','Final Time Offset');
    h_axes(3) = axes('parent',h_fig(3));
    h_plot = zeros(1,Nc);
    legend_str = cell(1,Nc);
    plot_mode = [0 0 0; hsv(7)];
    for wf_adc = 1:Nc
      h_plot(wf_adc) = plot(h_axes(3),ct_smooth(peak_offset(wf_adc,:)*dt/Mt*1e9,0.01), ...
        'Color', plot_mode(mod(wf_adc-1,length(plot_mode))+1,:), ...
        'LineStyle','none','Marker', '.');
      hold(h_axes(3), 'on');
      legend_str{wf_adc} = sprintf('%d-%d', wf_adc_list(wf_adc,1), wf_adc_list(wf_adc,2));
    end
    title(h_axes(3), 'Comparison of time offsets relative to ref');
    grid(h_axes(3), 'on');
    legend(h_plot,legend_str,'location','NorthEastOutside');
    ylabel(h_axes(3), 'Relative time (ns)');
    xlabel(h_axes(3), 'Range line');
    
    linkaxes(h_axes,'x');
    if ~isempty(param.collate_equal.rlines)
      xlim(h_axes(1), [min(param.collate_equal.rlines) max(param.collate_equal.rlines)]);
    else
      xlim(h_axes(1), [1 size(wf_data,2)]);
    end
    
    pos = get(h_fig,'Position');
    set(h_fig(1),'Position',[pos{1}(1:2) 700 pos{1}(4)]);
    set(h_fig(2),'Position',[pos{2}(1:2) 700 pos{2}(4)]);
    set(h_fig(3),'Position',[pos{3}(1:2) 700 pos{3}(4)]);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_amp_img_%02d',param.collate_equal.out_dir,img)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    fig_fn_dir = fileparts(fig_fn);
    if ~exist(fig_fn_dir,'dir')
      mkdir(fig_fn_dir);
    end
    saveas(h_fig(1),fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_amp_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(1),fig_fn);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_phase_img_%02d',param.collate_equal.out_dir,img)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(2),fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_phase_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(2),fig_fn);
    
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_time_img_%02d',param.collate_equal.out_dir,img)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(3),fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_time_img_%02d',param.collate_equal.out_dir,img)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    saveas(h_fig(3),fig_fn);
    
    set(h_fig(1),'Position',pos{1});
    set(h_fig(2),'Position',pos{2});
    set(h_fig(3),'Position',pos{3});
    
    if enable_visible_plot
      keyboard
    end
  end
  
  %% Create outputs
  % =========================================================================
  
  % Sort results according to rx_paths since this is how they must be
  % entered into the parameter spreadsheet/structure.
  [rx_paths,rx_path_idxs] = sort(rx_paths);
  peak_offset = peak_offset(rx_path_idxs,:);
  peak_val = peak_val(rx_path_idxs,:);
  wf_adc_list = wf_adc_list(rx_path_idxs,:);
  Tsys = Tsys(rx_path_idxs);
  chan_equal_deg = chan_equal_deg(rx_path_idxs);
  chan_equal_dB = chan_equal_dB(rx_path_idxs);
  ref_wf_adc_idx = find(ref_wf_adc_idx == rx_path_idxs);
  
  rx_paths_all = 1:max(rx_paths);
  rx_paths_missing = setdiff(rx_paths_all,rx_paths);
  
  peak_offset(rx_paths,:) = peak_offset;
  peak_val(rx_paths,:) = peak_val;
  wf_adc_list(rx_paths,:) = wf_adc_list;
  Tsys(rx_paths) = Tsys;
  chan_equal_deg(rx_paths) = chan_equal_deg;
  chan_equal_dB(rx_paths) = chan_equal_dB;
  
  peak_offset(rx_paths_missing,:) = NaN;
  peak_val(rx_paths_missing,:) = NaN;
  wf_adc_list(rx_paths_missing,:) = NaN;
  Tsys(rx_paths_missing) = NaN;
  chan_equal_deg(rx_paths_missing) = NaN;
  chan_equal_dB(rx_paths_missing) = NaN;
  
  equal.peak_offset = peak_offset;
  equal.peak_val = peak_val;
  
  equal.Tsys_offset = nanmean(peak_offset(:,rlines),2)*dt;
  equal.chan_equal_deg_offset = angle(nanmean(peak_val(:,rlines),2)) * 180/pi;
  equal.chan_equal_dB_offset = lp(nanmean(abs(peak_val(:,rlines)).^2,2),1);
  
  equal.Tsys_offset_std = nanstd(peak_offset(:,rlines),[],2)*dt;
  equal.chan_equal_deg_offset_std = angle(nanstd(peak_val(:,rlines),[],2)) * 180/pi;
  equal.chan_equal_dB_offset_std = lp(nanstd(abs(peak_val(:,rlines)).^2,[],2),1);
  
  equal.Tsys_offset = reshape(equal.Tsys_offset,[1 Nc]);
  equal.chan_equal_deg_offset = reshape(equal.chan_equal_deg_offset,[1 Nc]);
  equal.chan_equal_dB_offset = reshape(equal.chan_equal_dB_offset,[1 Nc]);
  
  equal.Tsys_offset_std = reshape(equal.Tsys_offset_std,[1 Nc]);
  equal.chan_equal_deg_offset_std = reshape(equal.chan_equal_deg_offset_std,[1 Nc]);
  equal.chan_equal_dB_offset_std = reshape(equal.chan_equal_dB_offset_std,[1 Nc]);
  
  equal.Tsys = Tsys + equal.Tsys_offset;
  equal.chan_equal_deg = chan_equal_deg + equal.chan_equal_deg_offset;
  equal.chan_equal_dB = chan_equal_dB + equal.chan_equal_dB_offset;
  
  equal.chan_equal_deg = equal.chan_equal_deg - equal.chan_equal_deg(ref_wf_adc_idx);
  equal.chan_equal_deg = angle(exp(1i*equal.chan_equal_deg/180*pi))*180/pi;
  equal.chan_equal_dB = equal.chan_equal_dB - equal.chan_equal_dB(ref_wf_adc_idx);
  
  equal.old_Tsys_str = [mat2str(round(Tsys*1e9*100)/100), '/1e9'];
  equal.Tsys_str = [mat2str(round(equal.Tsys*1e9*100)/100), '/1e9'];
  
  equal.chan_equal_dB_str = mat2str(round(equal.chan_equal_dB*10)/10);
  
  equal.chan_equal_deg_str = mat2str(round(angle(exp(1i*equal.chan_equal_deg*pi/180))*180/pi*10)/10);
  
  % If rounded Tsys are inserted, this accounts for the phase shift expected
  % by doing this.
  equal.chan_equal_deg_with_Tsys = chan_equal_deg + equal.chan_equal_deg_offset + 360*fc*round((equal.Tsys-Tsys)*1e10)/1e10;
  
  equal.chan_equal_deg_with_Tsys_str = mat2str(round(angle(exp(1i*equal.chan_equal_deg_with_Tsys*pi/180))*180/pi*10)/10);
  
  %% Print Results
  % =========================================================================
  if any(strcmp('final',param.collate_equal.debug_plots))
    sw_version = current_software_version;
    
    diary_fn = [ct_filename_ct_tmp(param,'','collate_equal',sprintf('%s_table_img_%02d',param.collate_equal.out_dir,img)) '.txt'];
    fid = fopen(diary_fn,'wb');
    for fid = [1 fid]
      if fid == 1; fid_error = 2; else fid_error = fid; end;
      fprintf(fid,'%s rec:%d-%d ref_wf_adc_idx:%d git-hash:%s (%s)\n', param.day_seg, ...
        rlines(1), rlines(end), ref_wf_adc_idx, sw_version.rev, sw_version.cur_date_time);
      
      for wf_adc = rx_paths_all
        wf = wf_adc_list(wf_adc,1);
        adc = wf_adc_list(wf_adc,2);
        fprintf(fid,'rx %d\t', rx_paths_all(wf_adc));
      end
      fprintf(fid,'\n');
      for wf_adc = rx_paths_all
        wf = wf_adc_list(wf_adc,1);
        adc = wf_adc_list(wf_adc,2);
        fprintf(fid,'wf-adc\t', wf, adc);
      end
      fprintf(fid,'\n');
      for wf_adc = rx_paths_all
        wf = wf_adc_list(wf_adc,1);
        adc = wf_adc_list(wf_adc,2);
        fprintf(fid,'%2.0f-%3.0f\t', wf, adc);
      end
      fprintf(fid,'\n');
      fprintf(fid,'Offsets from Old Coefficients (rows: equal_dB, equal_deg, Tsys_ns)\n');
      fprintf(fid,'%.1f\t', equal.chan_equal_dB_offset); fprintf(fid,'\n');
      fprintf(fid,'%.1f\t', equal.chan_equal_deg_offset); fprintf(fid,'\n');
      fprintf(fid,'%.2f\t', 1e9*equal.Tsys_offset); fprintf(fid,'\n');
      fprintf(fid,'New Coefficients\n');
      fprintf(fid,'%.1f\t', equal.chan_equal_dB); fprintf(fid,'\n');
      fprintf(fid,'%.1f\t', equal.chan_equal_deg); fprintf(fid,'\n');
      fprintf(fid,'%.2f\t', 1e9*equal.Tsys); fprintf(fid,'\n');
      fprintf(fid,'New Coefficients if using the old Tsys (in spreadsheet format)\n');
      fprintf(fid,'%s\n',equal.chan_equal_dB_str);
      fprintf(fid,'%s\n',equal.chan_equal_deg_str);
      fprintf(fid,'%s\n',equal.old_Tsys_str);
      fprintf(fid,'New Coefficients if updating with this Tsys (in spreadsheet format)\n');
      fprintf(fid,'%s\n',equal.chan_equal_dB_str);
      fprintf(fid,'%s\n',equal.chan_equal_deg_with_Tsys_str);
      fprintf(fid,'%s\n',equal.Tsys_str);
    end
    fclose(fid);
    fprintf('Metric table: %s\n', diary_fn);
  end
  
  %% Save Results
  % =========================================================================
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
  
end

fprintf('  Done (%s)\n', datestr(now));
