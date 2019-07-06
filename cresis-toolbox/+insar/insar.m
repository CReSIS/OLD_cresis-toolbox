
if 1
  wf = 2;
  insar_mode = 2; % 1 to find equalization coefficients, 2 to process data
  
  equalization_rlines = [];

  if wf == 1
    rbins = [];
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
  elseif wf == 2
    rbins = 280:420;
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
  elseif wf == 3
    equalization = 10.^(zeros(1,15)/20) .* exp(1i*([122.5 121.6 126.7 105.5 130.1 120.1 128.7 -0.0 -134.2 121.1 36.4 125.8 -171.5 -1.0 128.4]-[42.8 45.5 48.2 50.9 53.6 56.3 58.9 0.0 2.7 5.4 8.1 93.6 96.3 99.0 101.7]/2)/180*pi);
    rbins = 420:500;
  end
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_005_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_005_wf%d.mat',wf));
%   end
%   if ispc
%     fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_067_wf%d.mat',wf));
%   else
%     fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_20140429_01_067_wf%d.mat',wf));
%   end
  if ispc
    fn = fullfile('X:/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_combine_wf%d.mat',wf));
    equalization = [equalization equalization equalization equalization];
  else
    fn = fullfile('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/',sprintf('rds_thule_combine_wf%d.mat',wf));
    equalization = [equalization equalization equalization equalization];
  end
  master_idx = 8;
  % For combined file

elseif 0
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/north.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/north.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 0
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/middle.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/middle.mat';
  end
  master_idx = 10;
  rbins = 20:300;
elseif 0
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/south.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/south.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 0
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_south.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_south.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 0
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_north.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_TOdtu/CSARP_insar/iceland_north.mat';
  end
  master_idx = 1;
  rbins = 20:300;
elseif 1
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/good_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/good_line.mat';
  end
  master_idx = 1;
  rbins = 20:100;
elseif 0
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/medium_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/medium_line.mat';
  end
  master_idx = 2;
  rbins = 20:100;
else
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/bad_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/bad_line.mat';
  end
  master_idx = 2;
  rbins = 20:100;
end

proj = geotiffinfo(ct_filename_gis([],fullfile('greenland','Landsat-7','Greenland_natural_90m.tif')));

physical_constants;

load(fn);

%% Plot Results
% =========================================================================
h_fig_map = figure(100); clf;
h_plot_map = [];
h_legend_map = {};
hold on;
axis('equal');
h_fig_elev = figure(101); clf;
h_plot_elev = [];
h_legend_elev = {};
hold on;
xlabel('Range line');
ylabel('WGS-84 elevation (m)');
grid on;

h_data_axes = [];
for pass_idx = 1:length(pass)
  figure(pass_idx); clf;
  set(pass_idx,'WindowStyle','docked')
  if isempty(rbins)
    rbins = 1:size(pass(pass_idx).data,1);
  end
  if 1
    imagesc(lp(pass(pass_idx).data(rbins,:)))
    colormap(1-gray(256));
  else
    % Form interferogram (couple options)
    complex_data = pass(pass_idx).data(rbins,:);
    % Plot interferogram
    imagesc(hsv_plot(complex_data,-90));
    colormap(hsv(256))
    h_colorbar = colorbar;
    caxis([-pi pi])
    set(get(h_colorbar,'ylabel'),'string','angle (rad)')
    
  end
  h_data_axes(end+1) = gca;
  
  % Apply GPS time offset
  if 0
    pass_idx = 5;
    time_offset = -5;
    pass(pass_idx).ecef = interp1(pass(pass_idx).gps_time,pass(pass_idx).ecef.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
    pass(pass_idx).x = interp1(pass(pass_idx).gps_time,pass(pass_idx).x.',pass(pass_idx).gps_time+time_offset,'linear','extrap').';
  end
  
  pass(pass_idx).ecef = pass(pass_idx).origin;
  for rline = 1:size(pass(pass_idx).origin,2)
    pass(pass_idx).ecef(:,rline) = pass(pass_idx).ecef(:,rline) ...
      + [pass(pass_idx).x(:,rline) pass(pass_idx).y(:,rline) pass(pass_idx).z(:,rline)]*pass(pass_idx).pos(:,rline);
  end
  [pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev] = ecef2geodetic(referenceEllipsoid('wgs84'), ...
    pass(pass_idx).ecef(1,:), pass(pass_idx).ecef(2,:), pass(pass_idx).ecef(3,:));
  
  figure(h_fig_map);
  if 0
    h_plot = plot(pass(pass_idx).lon, pass(pass_idx).lat,'.');
    color = get(h_plot,'Color');
    h_text = text(pass(pass_idx).lon(1), pass(pass_idx).lat(1), sprintf('%d', pass_idx), 'Color', color);
  else
    [pass(pass_idx).proj_x,pass(pass_idx).proj_y] = projfwd(proj,pass(pass_idx).lat,pass(pass_idx).lon);
    h_plot_map(end+1) = plot(pass(pass_idx).proj_x/1e3, pass(pass_idx).proj_y/1e3,'.');
    h_legend_map{end+1} = sprintf('%d',pass_idx);
    color = get(h_plot_map(end),'Color');
    h_text = text(pass(pass_idx).proj_x(1)/1e3, pass(pass_idx).proj_y(1)/1e3, sprintf('%d', pass_idx), 'Color', color);
    xlabel('X (km)');
    ylabel('Y (km)');
    grid on;
  end
  
  figure(h_fig_elev);
  if master_idx == pass_idx
    h_plot_elev(end+1) = plot(pass(pass_idx).elev,'LineWidth',2,'UserData',pass_idx);
  else
    h_plot_elev(end+1) = plot(pass(pass_idx).elev,'UserData',pass_idx);
  end
  h_legend_elev{end+1} = sprintf('%d',pass_idx');
end
linkaxes(h_data_axes,'xy');
legend(h_plot_map,h_legend_map);

%% Co-Register Results
% =========================================================================


if 1
  % Option 1: Use a single pass as the reference.
  ref = pass(master_idx);

else
  % Option 2: Take a single pass as master and construct along track vectors relative
  % to it. Fit a polynomial to all the data using the along track vector
  % from step 1. Use this as the reference. This may be useful to do if no
  % single pass follows the center of the tube of passes.
  
  % With reference, run SAR coord system in a special mode where the mean is
  % not taken across Lsar, but instead each point is directly passed to the
  % output.
end

along_track = geodetic_to_along_track(ref.lat,ref.lon,ref.elev);
ref.surface_bin = interp1(ref.wfs(ref.wf).time, 1:length(ref.wfs(ref.wf).time), ref.surface);

if 0
  h_fig_ref_idx = figure(102); clf;
  hold on;
end

data = [];
for pass_idx = 1:length(pass)
    
  pass(pass_idx).ref_idx = zeros(1,size(pass(pass_idx).origin,2));
  last_idx = 0;
  for rline = 1:size(pass(pass_idx).ecef,2)
    % offset: offset from position rline of pass pass_idx from every point
    %         on the ref line
    offset = bsxfun(@minus, pass(pass_idx).ecef(:,rline), ref.ecef);
    % dist: converts offset into along-track distance of the slave line
    dist = offset.'*pass(pass_idx).x(:,rline);
    % min_idx: finds the point on the reference line that is closest to
    %          this point
    [min_dist,min_idx] = min(abs(dist));
    %       if min_idx == last_idx
    %         keyboard
    %       end
    last_idx = min_idx;
    pass(pass_idx).ref_idx(rline) = min_idx;
    
    % x_offset: along-track offset on master line of the closest point
    x_offset = offset(:,min_idx).'*ref.x(:,min_idx);
    % Compute FCS of slave point in master line coordinate system
    %   FCS: flight (aka SAR) coordinate system
    pass(pass_idx).along_track(rline) = along_track(min_idx) + x_offset;
    pass(pass_idx).ref_y(rline) = offset(:,min_idx).'*ref.y(:,min_idx);
    pass(pass_idx).ref_z(rline) = offset(:,min_idx).'*ref.z(:,min_idx);

    if 0
      % Compute the location of all pixels from this range line in ECEF
      pass(pass_idx).wfs(pass(pass_idx).wf).time;
      time = pass(pass_idx).time(2)-pass(pass_idx).time(1);
      
      range = time * c/2;
      range(time>pass(pass_idx).surface(rline)) = pass(pass_idx).surface(rline)*c/2 ...
        + (time(time>pass(pass_idx).surface(rline)) - pass(pass_idx).surface(rline))*c/2/sqrt(er_ice);
      
      pixels = pass(pass_idx).ecef(:,rline) + pass(pass_idx).z(:,rline)*range;
      
      % Compute the location of all pixels from this range line in the master
      % line coordinate system FCS.
      Tfcs = [ref.x, ref.y, ref.z];
      pass(pass_idx).pixels(:,:,rline) = (pixels - ref.origin) / Tfcs;
    end
  end
  
  pass(pass_idx).along_track_slave = geodetic_to_along_track(pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev);;
  
  if 0
    %% Debug plot showing indexes for alignment of passes
    figure(h_fig_ref_idx);
    plot(pass(pass_idx).ref_idx)
    drawnow;
  end
  
  % Resample images and position vectors onto a common along-track axes
  % 1. Oversample slave data by 10x in along track
  Mx = 10;
  Nx = size(pass(pass_idx).data,2);
  data_oversample = interpft(pass(pass_idx).data.',Mx*Nx);
  % 2. Interpolate to find the oversampled slave axes
  along_track_oversample = interp1(0:Nx-1, ...
    pass(pass_idx).along_track, (0:Nx*Mx-1)/Mx,'linear','extrap');
  % 3. Interpolate oversampled slave data onto master along track axes
  pass(pass_idx).ref_data = interp1(along_track_oversample, ...
    data_oversample, along_track,'linear','extrap').';
  
  pass(pass_idx).ref_y = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_y, along_track,'linear','extrap').';
  pass(pass_idx).ref_z = interp1(pass(pass_idx).along_track, ...
    pass(pass_idx).ref_z, along_track,'linear','extrap').';

  if insar_mode == 1
    % Motion compensation of FCS z-motion
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*pass(pass_idx).wfs(pass(pass_idx).wf).freq*dt) );
    end
  elseif insar_mode == 2
    % Co-register images using GPS and nadir squint angle assumption
    %
    % Motion compensation of FCS z-motion without center frequency so there
    % is no phase shift.
    freq = pass(pass_idx).wfs(pass(pass_idx).wf).freq;
    freq = freq - freq(1); % Remove center frequency offset
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*freq*dt) );
    end
  elseif 0
    % Co-register images using cross-correlation
    keyboard
  end

  % Match time axis to master_idx
  pass(pass_idx).ref_data = interp1(pass(pass_idx).wfs(pass(pass_idx).wf).time, pass(pass_idx).ref_data, pass(master_idx).wfs(pass(master_idx).wf).time, 'linear', 0);
  
  if 0
    % Normalize surface phase
    Nt = size(pass(pass_idx).ref_data,1);
    Nx = size(pass(pass_idx).ref_data,2);
    H = pass(master_idx).ref_data(round(ref.surface_bin)+(0:Nx-1)*Nt) .* conj(pass(pass_idx).ref_data(round(ref.surface_bin)+(0:Nx-1)*Nt));
    H = exp(1i*angle(H));
    pass(pass_idx).ref_data = bsxfun(@times,pass(pass_idx).ref_data,H);
  end
  
  % Concatenate data into a single matrix
  data = cat(3,data,pass(pass_idx).ref_data);
end

% Apply equalization
% -----------------------
if insar_mode == 2
  equalization = reshape(equalization,[1 1 numel(equalization)]);
  data = bsxfun(@times,data,1./equalization);
end

h_data_axes = [];
new_equalization = [];
for pass_idx = 1:length(pass)
  % Plot interferograms
  % -----------------------
  figure(pass_idx); clf;
  set(pass_idx,'WindowStyle','docked')
  if 0
    imagesc(lp(data(rbins,:,pass_idx)))
    colormap(1-gray(256));
    ylabel('Range bin');
    xlabel('Range line');
    title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    caxis([-90 8]);
  else
    % Form interferogram (couple options)
    complex_data = fir_dec(data(rbins,:,pass_idx) .* conj(data(rbins,:,master_idx)),ones(1,11)/11,1);
    if isempty(equalization_rlines)
      new_equalization(pass_idx) = mean(complex_data(:)); % equalization only valid when motion compensation with phase is used
    else
      new_equalization(pass_idx) = mean(mean(complex_data(:,equalization_rlines))); % equalization only valid when motion compensation with phase is used
    end
    if insar_mode == 1
      complex_data = complex_data ./ new_equalization(pass_idx);
    end
    % Plot interferogram
    imagesc(hsv_plot(complex_data,-90));
    colormap(hsv(256))
    h_colorbar = colorbar;
    caxis([-pi pi])
    set(get(h_colorbar,'ylabel'),'string','angle (radians)');
    ylabel('Range bin');
    xlabel('Range line');
    title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    
    if 0
      coherence = fir_dec(data(rbins,:,pass_idx) .* conj(data(rbins,:,master_idx)) ./ abs(data(rbins,:,pass_idx)) ./ abs(data(rbins,:,master_idx))  ,ones(1,11)/11,1);
      imagesc(abs(coherence));
      colormap(1-gray(256))
      h_colorbar = colorbar;
      caxis([0 1])
      set(get(h_colorbar,'ylabel'),'string','Coherence');
      ylabel('Range bin');
      xlabel('Range line');
      title(sprintf('%s_%03d %d',pass(pass_idx).param_sar.day_seg,pass(pass_idx).param_sar.load.frm,pass(pass_idx).direction),'interpreter','none')
    end
    
  end
  h_data_axes(end+1) = gca;
  
end
fprintf('%.1f ', lp(new_equalization)-mean(lp(new_equalization)));
fprintf('\n');
fprintf('%.1f ', angle(new_equalization)*180/pi)
fprintf('\n');
linkaxes(h_data_axes,'xy');
if insar_mode == 1
  return
end

%% Array Processing

% Package data to call array_proc.m
% 1. Data
% 2. Trajectory and attitude
% 3. Array processing parameters
data = {permute(data,[1 2 4 5 3])};

param.array = [];
param.array.method = 1;
param.array.Nsv = 128;
% param.array = rmfield(param.array,'Nsv');
param.array.theta = linspace(-3,3,128);
% param.array.Nsv = {'theta',linspace(-10,10,256)/180*pi};
param.array.Nsrc = 2;
param.array.bin_rng = [-2:2];
param.array.line_rng = [-20:20];
param.array.dbin = 1;
param.array.dline = 11;
param.array.freq_rng = 1;
h_fig_baseline = figure(200); clf;
h_plot_baseline = [];
h_legend_baseline = {};
for pass_idx = 1:length(pass)
  param.array.fcs{1}{pass_idx}.pos = along_track;
  param.array.fcs{1}{pass_idx}.pos(2,:) = pass(pass_idx).ref_y;
  param.array.fcs{1}{pass_idx}.pos(3,:) = pass(pass_idx).ref_z;
  param.array.fcs{1}{pass_idx}.base_line ...        
    = sqrt( (pass(pass_idx).ref_z - pass(master_idx).ref_z).^2 ...
      + (pass(pass_idx).ref_y - pass(master_idx).ref_y).^2 );
    
  h_plot_baseline(end+1) = plot(param.array.fcs{1}{pass_idx}.base_line);
  h_legend_baseline{end+1} = sprintf('%d',pass_idx);
  hold on;

  param.array.fcs{1}{pass_idx}.surface = ref.surface;
end
xlabel('Range line');
ylabel('Baseline (m)');
grid on;
legend(h_plot_baseline,h_legend_baseline);

param.array.wfs.time = ref.wfs(ref.wf).time;
dt = param.array.wfs.time(2)-param.array.wfs.time(1);
param.array_proc.bin0 = param.array.wfs.time/dt;
param.array.sv_fh = @array_proc_sv;
param.array.wfs.fc = ref.wfs(ref.wf).fc;
param.array.imgs = {[ones(length(pass),1), (1:length(pass)).']};
param.array.tomo_en = true;

%%
array_proc_methods;
param = array_proc(param);
param.array.method = STANDARD_METHOD;
[param_array0,result0] = array_proc(param,data);
% param.array.method = MVDR_METHOD;
% [param_array1,result1] = array_proc(param,data);
param.array.method = MUSIC_METHOD;
[param_array2,result2] = array_proc(param,data);

figure(101); clf;
imagesc(lp(result0.img))
title('Periodogram')
colormap(1-gray(256))
h_axes = gca;

% figure(102); clf;
% imagesc(lp(result1.img))
% title('MVDR')
% colormap(1-gray(256))
% h_axes(end+1) = gca;

figure(103); clf;
imagesc(lp(result2.img))
title('MUSIC')
colormap(1-gray(256))
h_axes(end+1) = gca;

figure(104); clf;
imagesc(lp(fir_dec(abs(data{1}(:,:,master_idx)).^2, ones(size(param.array.line_rng)), param.array.dline, ...
  1-param.array.line_rng(1), size(data{1}(:,:,master_idx),2)-length(param.array.line_rng)+1)))
title('Single Channel');
colormap(1-gray(256))
h_axes(end+1) = gca;

linkaxes(h_axes,'xy');

%% Save Results
[fn_dir,fn_name] = fileparts(fn);

Tomo = result0.tomo;
Data = result0.img;
GPS_time = ref.gps_time(param_array0.array_proc.lines);
Latitude = ref.lat(param_array0.array_proc.lines);
Longitude = ref.lon(param_array0.array_proc.lines);
Elevation = ref.elev(param_array0.array_proc.lines);
Roll = ref.roll(param_array0.array_proc.lines);
Pitch = ref.pitch(param_array0.array_proc.lines);
Heading = ref.heading(param_array0.array_proc.lines);
Surface = ref.surface(param_array0.array_proc.lines);
Bottom = nan(size(Surface));
param_sar = pass(master_idx).param_sar;
param_records = pass(master_idx).param_records;
param_array = param_array0;
Time = param.array.wfs.time(param_array0.array_proc.bins);
file_version = '1';
fn_mat = fullfile(fn_dir,[fn_name '_standard.mat']);
save('-v7.3',fn_mat,'Tomo','Data','Latitude','Longitude','Elevation','GPS_time', ...
  'Surface','Bottom','Time','param_array','param_records', ...
  'param_sar', 'Roll', 'Pitch', 'Heading', 'file_version');

Tomo = result2.tomo;
Data = result2.img;
GPS_time = ref.gps_time(param_array2.array_proc.lines);
Latitude = ref.lat(param_array2.array_proc.lines);
Longitude = ref.lon(param_array2.array_proc.lines);
Elevation = ref.elev(param_array2.array_proc.lines);
Roll = ref.roll(param_array2.array_proc.lines);
Pitch = ref.pitch(param_array2.array_proc.lines);
Heading = ref.heading(param_array2.array_proc.lines);
Surface = ref.surface(param_array2.array_proc.lines);
Bottom = nan(size(Surface));
param_sar = pass(master_idx).param_sar;
param_records = pass(master_idx).param_records;
param_array = param_array2;
Time = param.array.wfs.time(param_array2.array_proc.bins);
file_version = '1';
fn_mat = fullfile(fn_dir,[fn_name '_music.mat']);
save('-v7.3',fn_mat,'Tomo','Data','Latitude','Longitude','Elevation','GPS_time', ...
  'Surface','Bottom','Time','param_array','param_records', ...
  'param_sar', 'Roll', 'Pitch', 'Heading', 'file_version');

fn_map = fullfile(fn_dir,[fn_name '_map.fig']);
saveas(h_fig_map,fn_map);
fn_elev = fullfile(fn_dir,[fn_name '_elev.fig']);
saveas(h_fig_elev,fn_elev);
fn_baseline = fullfile(fn_dir,[fn_name '_baseline.fig']);
saveas(h_fig_baseline,fn_baseline);
fn_standard = fullfile(fn_dir,[fn_name '_standard.fig']);
saveas(101,fn_standard);
% fn_mvdr = fullfile(fn_dir,[fn_name '_mvdr.fig']);
% saveas(102,fn_mvdr);
fn_music = fullfile(fn_dir,[fn_name '_music.fig']);
saveas(103,fn_music);
fn_single = fullfile(fn_dir,[fn_name '_single.fig']);
saveas(104,fn_single);

return



%% Other plots and plot setup


figure(500); clf;
for idx = 1800:-1:600%1:50:size(result1.img,3)
  imagesc(lp(result0.img(20:120,:,idx)))
  colormap(jet(256));
  caxis([-80 -20])
  idx
  pause
end

%%

figure; plot(array_param2.theta*180/pi,lp(result1.img(20+54-1,:,1291)))
grid on
xlim([-90 90])
xlabel('Direction of arrival (deg)');
ylabel('Relative power (dB)');
title('Bad line, MVDR passes [2 4 6], rline 1291');

figure(h_fig_ref_idx); h_axes = gca;
%figure(h_fig_map); h_axes(end+1) = gca;
figure(h_fig_elev); h_axes(end+1) = gca;
linkaxes(h_axes,'x');

h_fig_combined = figure(103); clf;
imagesc(lp(mean(data(rbins,:,[1 6]),3)));
colormap(1-gray(256));
set(h_fig_combined,'WindowStyle','Docked');
hold on;
plot(ref.surface_bin-rbins(1)+1);

surf_data = zeros(length(pass),size(data,2));
for pass_idx = 1:length(pass)
  for rline = 1:size(data,2)
    surf_data(pass_idx,rline) = data(round(ref.surface_bin(rline)),rline,pass_idx);
  end
end

h_fig_angle = figure(104); clf;
complex_data = surf_data(6,:) .* conj(surf_data(1,:));
complex_data(lp(complex_data) < -96) = NaN;
plot(180/pi*angle(complex_data),'.')
xlabel('Range line');
ylabel('Angle (deg)');
grid on;
