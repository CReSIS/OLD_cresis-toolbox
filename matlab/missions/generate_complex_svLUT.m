% script generate_complex_svLUT.m
%
% Generates complex steering vector lookup table using roll data collected
% with an array system. First run collate_equal (coh_noise_tracker) on the
% dataset and then run this function. Transmit equalization mode should be
% used (one transmitter at a time).
%
% This script assumes that transmit equalization has been applied.
%
% Example: Run from run_generate_complex_svLUT.m
%
% Author: John Paden

physical_constants;

% Only supports a single image right now
img = 1;

%% Load waveform data
% =====================================================================
img = param.collate_equal.imgs{imgs_idx}(1);
wf_adcs = [];
for sub_img_idx = 1:length(param.collate_equal.imgs{imgs_idx})
  sub_img = param.collate_equal.imgs{imgs_idx}(sub_img_idx);
  
  if isempty(param.collate_equal.wf_adc_idxs)
    wf_adc_idxs = 1:size(param.analysis.imgs{sub_img},1);
  else
    wf_adc_idxs = param.collate_equal.wf_adc_idxs{imgs_idx}{sub_img_idx};
  end
  for wf_adc_idx = wf_adc_idxs
    wf = param.analysis.imgs{sub_img}(wf_adc_idx,1);
    adc = param.analysis.imgs{sub_img}(wf_adc_idx,2);
    wf_adcs = [wf_adcs; [wf adc]];
    
    % Load the waveform file
    fn_dir = fileparts(ct_filename_out(param,param.collate_equal.in_dir));
    fn = fullfile(fn_dir,sprintf('waveform_%s_wf_%d_adc_%d.mat', param.day_seg, wf, adc));
    fprintf('Loading %s (%s)\n', fn, datestr(now));
    waveform = load(fn);
    if wf_adc_idx == wf_adc_idxs(1)
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
param.gps_source = waveform.param_records.gps_source;

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

% Dimensions
Nt = size(wf_data,1);
Nx = size(wf_data,2);
Nc = size(wf_data,3);

% Determine which range lines will be used
rlines = param.collate_equal.rlines;
all_rlines = 1:Nx;
if isempty(rlines)
  rlines = all_rlines;
else
  rlines = intersect(rlines,all_rlines);
end

% Undo any transmit windowing that was applied
for wf_adc_idx = 1:Nc
  wf_data(:,:,wf_adc_idx) = wf_data(:,:,wf_adc_idx) ./ Hchan(wf_adc_idx);
end

if 0
  % DEBUG
  surf_bins = zero_surf_bin;
  figure(1); clf;
  subplot(3,1,1);
  plot(lp(wf_data(surf_bins,:,1)),'.')
  a1 = gca;
  subplot(3,1,2);
  plot(angle(wf_data(surf_bins,:,1) .* conj(wf_data(surf_bins,:,5)) ),'.')
  a2 = gca;
  subplot(3,1,3);
  plot(roll.'*180/pi);
  a3 = gca;
  linkaxes([a1 a2 a3],'x');
  figure(2); clf;
  imagesc(lp(wf_data(:,:,1)));
  return;
end

%% Apply motion compensation and channel equalization
% =========================================================================

%% 1. Create the wf-adc to receiver path mapping
rx_paths = zeros(1,size(wf_adcs,1));
for wf_adc_idx = 1:Nc
  wf = wf_adcs(wf_adc_idx,1);
  adc = wf_adcs(wf_adc_idx,2);
  rx_paths(wf_adc_idx) = param.radar.wfs(wf).rx_paths(adc);
end

%% 2. Determine time delay and phase correction for position and channel equalization
% drange: increased range is positive
% dtime: increased range is positive
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

if 0
  fig_offset = 10;
  chan_mapping = [1:Nc];
  rbin = zero_surf_bin; % Surface bin (look at imagesc(lp(wf_data(:,:,1))))
  figure(1); clf;
  imagesc(squeeze(angle(exp(-1i*2*pi*bsxfun(@times,fc,dtime)))).');
  h_axes = gca;
  title('Correction');
  
  figure(fig_offset+2); clf;
  myref = wf_data(rbin,:,chan_mapping(ref_ant));
  dd=(fir_dec(squeeze(bsxfun(@times,wf_data(rbin,:,chan_mapping),conj(myref))).',hanning(5).',1) );
  dd = bsxfun(@times,dd, exp(-j*angle(mean(dd,2))));
  var(angle(dd).')
  imagesc(angle(dd));
  h_axes(end+1) = gca;
  title('Angle Estimate');
  linkaxes(h_axes,'xy');
  
  figure(4); clf;
  plot(roll(1,:).'*180/pi);
  hold on;
  plot(pitch(1,:).'*180/pi);
  legend('Roll','Pitch');

  % Re-arrange correction terms
  %dtime = dtime(:,:,[4 5 6 1 2 3]);
end

%% 3. Apply time delay, phase, and amplitude correction
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

%% Retrack surface and circshift data so that surface lies in range bin 6
% =========================================================================
ref_bin = 6;
surf_bin = NaN*zeros(1,Nx);
ml_data = lp(fir_dec(abs(mean(wf_data,3)).^2,ones(1,5)/5,1));
for rline = 1:Nx
  cur_threshold = max([ml_data(1,rline)+7; ml_data(:,rline)-13]);
  tmp = find(ml_data(:,rline) > cur_threshold,1);
  if ~isempty(tmp)
    [~,max_offset] = max(ml_data(tmp+(0:2),rline));
    tmp = tmp-1 + max_offset;
    surf_bin(rline) = tmp;
  end
end

if debug_level == 2
  figure(2); clf;
  imagesc(ml_data);
  hold on;
  plot(surf_bin);
  xlabel('Range line');
  ylabel('Range bin');
  title('Multilooked data with surface tracker');
  h_axes = gca;
  figure(3); clf;
  plot(roll(1,:)*180/pi)
  linkaxes([h_axes gca],'x');
end

for rline = 1:Nx
  if ~isnan(surf_bin(rline))
    wf_data(:,rline,:) = circshift(wf_data(:,rline,:),[ref_bin-surf_bin(rline) 0 0]);
  end
end

if debug_level >= 1
  ml_data = lp(fir_dec(abs(mean(wf_data,3)).^2,ones(1,5)/5,1));
  number_h_fig = 1; figure(number_h_fig); clf;
  imagesc(ml_data);
  xlabel('Range line');
  ylabel('Range bin');
  title(sprintf('Multilooked data with surface aligned to range bin %d', ref_bin));
end

if debug_level == 2
  return
end

%% Create Lookup Table
% =========================================================================
% Bin all the roll measurements into 1 degree large bins
% Determine the roll angle bins that we will put each measurement into
[roll_binned,~,roll_idxs] = unique(round((roll(1,rlines))*180/pi));

figure(2); clf;
hist(roll_binned(roll_idxs), -90:90);
xlabel('Angle (deg)');
ylabel('Number of measurements');

% sv_table:
%  First row: roll angle
%  Rows 2 to Nc+1: complex values of surface return for channels 1 to Nc
%  Each column represents a different roll angle from roll_binned
sv_table = zeros(Nc, length(roll_binned));
power_table = zeros(size(sv_table));

for ant = 1:Nc
  powers = max(abs(double(wf_data(ref_bin+(-1:3),rlines,ant))).^2,[],1);
  complex_vals = mean(wf_data(ref_bin+[0:0],rlines,ant) ...
    .* exp(-1i*angle(wf_data(ref_bin+[0:0],rlines,ref_ant))),1);
  
  % Average all the data falling within each angle/roll bin
  for roll_idx = 1:length(roll_binned)
    sv_table(ant,roll_idx) = mean(complex_vals(roll_idxs == roll_idx));
    power_table(ant,roll_idx) = mean(powers(roll_idxs == roll_idx));
  end
end

if 0
  %% DEBUG
  close all
  figure(1); clf;
  imagesc(roll_binned,[],lp(power_table));
  colorbar;
  figure(2); clf;
  imagesc(roll_binned,[],angle(sv_table));
  colorbar;
  figure(30); clf;
  imagesc(lp(wf_data(:,:,1) ))
  return;
end
sv_table = sqrt(power_table) .* exp(1i*angle(sv_table));

%% Cull results (remove bins with too few samples or too large roll angle)
% =========================================================================
N = hist(roll_binned(roll_idxs), roll_binned);
good_mask = N > good_mask_min_samples;
good_mask = good_mask & roll_binned > good_mask_min_angle & roll_binned < good_mask_max_angle;
if 0
  %% DEBUG
  plot(roll_binned, good_mask)
end
roll_binned = roll_binned(good_mask);
sv_table = sv_table(:,good_mask);
power_table = power_table(:,good_mask);

%% Receiver equalization (force sv to be all ones at nadir)
% =========================================================================
nadir_idx = find(roll_binned==0);
rx_equalization = sv_table(:,nadir_idx);
fprintf('Equalization (deg):\n');
fprintf('%5.1f ', angle(rx_equalization)*180/pi);
fprintf('\n');
fprintf('Equalization (dB):\n');
fprintf('%5.1f ', 10*log10(abs(rx_equalization)));
fprintf('\n');
sv_table = sv_table ./ repmat(rx_equalization,[1 size(sv_table,2)]);
sv_deviation = sv_table;

%% Plot raw results
if debug_level == 1
  figure(3); clf;
  imagesc(roll_binned,[],angle(sv_deviation)*180/pi)
  h = colorbar;
  set(get(h,'YLabel'),'String','Phase (deg)')
  xlabel('Angle (deg)');
  ylabel('Antenna');
  
  figure(4); clf;
  plot(roll_binned,angle(sv_deviation)*180/pi)
  ylabel('Phase deviation (deg)')
  xlabel('Angle (deg)');
  legend_str = {};
  for ant = 1:Nc
    legend_str{ant} = sprintf('Ant %d', ant);
  end
  legend(legend_str,'location','best');
  
  figure(5); clf;
  plot(roll_binned,10*log10(abs(sv_deviation).^2))
  ylabel('Power (dB)');
  xlabel('Angle (deg)');
  for ant = 1:Nc
    legend_str{ant} = sprintf('Ant %d', ant);
  end
  legend(legend_str,'location','best');
end

%% Fit results to spatial filter
% =========================================================================
ky_relative = sind(roll_binned);
sv_deviation_filter = [];
sv_deviation_golay = [];
for ant = 1:size(sv_deviation,1)
  [B,A] = invfreqz(sv_deviation(ant,:),ky_relative,'complex',degrees_of_freedom,0);
  sv_deviation_filter(ant,:) = freqz(B,A,ky_relative);
  sv_deviation_golay(ant,:) = sgolayfilt(sv_deviation(ant,:),2,17,hanning(17));
  
  if debug_level == 1
    figure(100+ant); clf;
    plot(roll_binned,10*log10(abs(sv_deviation_filter(ant,:)).^2))
    hold on;
    plot(roll_binned,10*log10(abs(sv_deviation_golay(ant,:)).^2),'r-')
    plot(roll_binned,10*log10(abs(sv_deviation(ant,:)).^2),'.')
    grid on;
    ylabel('Power (dB)');
    xlabel('Angle (deg)');
    legend('Filter Fit','Savitzky-Golay Fit','Measured','Location','best');
    set(100+ant,'WindowStyle','docked')
    title(sprintf('Ant %d', ant));
    
    figure(200+ant); clf;
    plot(roll_binned,180/pi*angle(sv_deviation_filter(ant,:)))
    hold on;
    plot(roll_binned,180/pi*angle(sv_deviation_golay(ant,:)),'r-')
    plot(roll_binned,180/pi*angle(sv_deviation(ant,:)),'.')
    grid on;
    ylabel('Phase deviation (deg)')
    xlabel('Angle (deg)');
    legend('Filter Fit','Savitzky-Golay Fit','Measured','Location','best');
    set(200+ant,'WindowStyle','docked')
    title(sprintf('Ant %d', ant));
  end
end
if strcmp(fit_method,'filter')
  % Spatial filter fit
  fit_method_string = 'Spatial filter';
  sv_deviation_fit = sv_deviation_filter;
elseif strcmp(fit_method,'sgolay')
  % Savitsky-Golay filter fit
  fit_method_string = 'Savitsky-Golay';
  sv_deviation_fit = sv_deviation_golay;
else
  error('Unsupported filter method %s', fit_method);
end

%% Load the ideal steering vectors
physical_constants;
clear phase_centers;
for ant = 1:Nc
  wf = wf_adcs(ant,1);
  adc = wf_adcs(ant,2);
  tx_weights = param.radar.wfs(wf).tx_weights;
  rx_chan = param.radar.wfs(wf).rx_paths(adc);
  phase_centers(:,ant) = lever_arm(param, tx_weights, rx_chan);
end
[theta,sv_ideal] = array_proc_sv({'theta' roll_binned/180*pi}, fc, phase_centers(2,:).', phase_centers(3,:).');
sv_ideal = bsxfun(@(x,y) x./y, sv_ideal, sv_ideal(ref_ant,:));
rx_equalization = sv_ideal(:,nadir_idx);
sv_ideal = bsxfun(@(x,y) x./y, sv_ideal, rx_equalization);

%% Normalize angles and take sqrt of each result
% Normalize angle relative to zero roll
% Take the sqrt since the measurement includes both transmission and receiption
% through the same antenna and we want to know just one way.

sv_deviation_fit_angle = unwrap(angle(sv_deviation_fit).').';
sv_deviation_fit_angle = sv_deviation_fit_angle - repmat(sv_deviation_fit_angle(:,nadir_idx),[1 size(sv_deviation_fit_angle,2)]);
sv_deviation_fit = sqrt(abs(sv_deviation_fit)) .* exp(1i*0.5*sv_deviation_fit_angle);

sv_deviation_angle = unwrap(angle(sv_deviation).').';
sv_deviation_angle = sv_deviation_angle - repmat(sv_deviation_angle(:,nadir_idx),[1 size(sv_deviation_angle,2)]);
sv_deviation = sqrt(abs(sv_deviation)) .* exp(1i*0.5*sv_deviation_angle);

sv_deviation = bsxfun(@(x,y) x./y, sv_deviation, sv_deviation_fit(:,nadir_idx));
sv_deviation_fit = bsxfun(@(x,y) x./y, sv_deviation_fit, sv_deviation_fit(:,nadir_idx));

sv_ideal_angle = unwrap(angle(sv_ideal).').';
sv_ideal_angle = sv_ideal_angle - repmat(sv_ideal_angle(:,nadir_idx),[1 size(sv_ideal_angle,2)]);
sv_ideal = sqrt(abs(sv_ideal)) .* exp(1i*0.5*sv_ideal_angle);

%% Plot final results
if debug_level == 1
  number_h_fig = 6; figure(number_h_fig); clf;
  colors = hsv(Nc);
  for ant=1:Nc
    plot(roll_binned,10*log10(abs(sv_deviation_fit(ant,:) .* sv_ideal(ant,:)).^2),'Color',colors(ant,:));
    hold on;
  end
  for ant=1:Nc
    plot(roll_binned,10*log10(abs( sv_deviation(ant,:) .* sv_ideal(ant,:) ).^2),'.','Color',colors(ant,:));
  end
  grid on;
  ylabel('Power (dB)');
  xlabel('Angle (deg)');
  for ant = 1:Nc
    legend_str{ant} = sprintf('Ant %d', ant);
  end
  legend(legend_str,'location','best');
  xlim([plot_min_angle plot_max_angle]);
  title(sprintf('%s fitting method',fit_method_string));
  fprintf('Figure %d: Steering vector amplitude\n  Black: ideal, Colored: fit to measurements, Dots: measurements\n', number_h_fig);
  
  number_h_fig = 7; figure(number_h_fig); clf;
  nadir_idx = find(roll_binned==0);
  final_fit = unwrap(angle(sv_deviation_fit .* sv_ideal).').';
  final_fit = final_fit - repmat(final_fit(:,nadir_idx),[1 size(final_fit,2)]);
  final_meas = unwrap(angle(sv_deviation .* sv_ideal).').';
  final_meas = final_meas - repmat(final_meas(:,nadir_idx),[1 size(final_meas,2)]);
  ideal = unwrap(angle(sv_ideal).').';
  ideal = ideal - repmat(ideal(:,nadir_idx),[1 size(ideal,2)]);
  colors = hsv(Nc);
  for ant = 1:Nc
    plot(roll_binned,final_fit(ant,:)*180/pi,'Color',colors(ant,:));
    hold on;
  end
  for ant = 1:Nc
    plot(roll_binned,final_meas(ant,:)*180/pi,'.','Color',colors(ant,:));
    plot(roll_binned,ideal(ant,:)*180/pi,'k-');
  end
  grid on;
  ylabel('Phase (deg)')
  xlabel('Angle (deg)');
  for ant = 1:Nc
    legend_str{ant} = sprintf('Ant %d', ant);
  end
  legend(legend_str,'location','best');
  xlim([plot_min_angle plot_max_angle]);
  title(sprintf('%s fitting method',fit_method_string));
  fprintf('Figure %d: Steering vector phases\n  Black: ideal, Colored: fit to measurements, Dots: measurements\n', number_h_fig);
end

%% Save Results
% =========================================================================
fprintf('Saving %s\n', output_fn);
output_fn_dir = fileparts(output_fn);
if ~exist(output_fn_dir,'dir')
  mkdir(output_fn_dir);
end
save(output_fn,'surf_bin','tx_ant','roll_binned','sv_deviation','sv_deviation_fit','sv_ideal','param');

return;
