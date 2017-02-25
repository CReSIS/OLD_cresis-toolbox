
if 1
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/good_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/good_line.mat';
  end
  master_idx = 1;
elseif 1
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/medium_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/medium_line.mat';
  end
  master_idx = 2;
else
  if ispc
    fn = 'X:/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/bad_line.mat';
  else
    fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/bad_line.mat';
  end
  master_idx = 2;
end

proj = geotiffinfo(ct_filename_gis([],fullfile('greenland','Landsat-7','Greenland_natural_90m.tif')));

physical_constants;

load(fn);
rbins = 20:120;

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
    h_plot_elev(end+1) = plot(pass(pass_idx).elev,'LineWidth',2);
  else
    h_plot_elev(end+1) = plot(pass(pass_idx).elev);
  end
  h_legend_elev{end+1} = sprintf('%d',pass_idx');
  
  pass(pass_idx).ecef = [];
  [pass(pass_idx).ecef(1,:),pass(pass_idx).ecef(2,:),pass(pass_idx).ecef(3,:)] ...
    = geodetic2ecef(pass(pass_idx).lat/180*pi, pass(pass_idx).lon/180*pi, pass(pass_idx).elev, WGS84.ellipsoid);
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
ref.surface_bin = interp1(ref.wfs.time, 1:length(ref.wfs.time), ref.surface);

if 0
  h_fig_ref_idx = figure(102); clf;
  hold on;
end

h_data_axes = [];
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

  if 1
    % Motion compensation of FCS z-motion
    for rline = 1:size(pass(pass_idx).ref_data,2)
      % Convert z-offset into time-offset assuming nadir DOA
      dt = pass(pass_idx).ref_z(rline)/(c/2);
      pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
        .*exp(1i*2*pi*pass(end).wfs.freq*dt) );
    end
  end
  
  % Concatenate data into a single matrix
  data = cat(3,data,pass(pass_idx).ref_data);
  
  % Plot interferograms
  % -----------------------
  figure(pass_idx); clf;
  set(pass_idx,'WindowStyle','docked')
  if 0
    imagesc(lp(pass(pass_idx).ref_data(rbins,:)))
    colormap(1-gray(256));
  else
    % Form interferogram (couple options)
    complex_data = pass(pass_idx).ref_data(rbins,:) .* conj(ref.data(rbins,:));
    % Plot interferogram
    imagesc(hsv_plot(complex_data,-50));
    colormap(hsv(256))
    h_colorbar = colorbar;
    caxis([-pi pi])
    set(get(h_colorbar,'ylabel'),'string','angle (radians)');
    ylabel('Range bin');
    xlabel('Range line');
    
  end
  h_data_axes(end+1) = gca;
  
end
linkaxes(h_data_axes,'xy');

%% Array Processing

% Package data to call array_proc.m
% 1. Data
% 2. Trajectory and attitude
% 3. Array processing parameters
data = {permute(data,[1 2 4 5 3])};

array_param.method = 1;
array_param.Nsv = 64;
array_param.Nsig = 2;
array_param.bin_rng = [0];
array_param.rline_rng = [-21:21];
array_param.dbin = 1;
array_param.dline = 1;
array_param.freq_rng = 1;
h_fig_baseline = figure(200); clf;
h_plot_baseline = [];
h_legend_baseline = {};
for pass_idx = 1:length(pass)
  array_param.fcs{1}{pass_idx}.pos = along_track;
  array_param.fcs{1}{pass_idx}.pos(2,:) = pass(pass_idx).ref_y;
  array_param.fcs{1}{pass_idx}.pos(3,:) = pass(pass_idx).ref_z;
  array_param.fcs{1}{pass_idx}.base_line ...        
    = sqrt( (pass(pass_idx).ref_z - pass(master_idx).ref_z).^2 ...
      + (pass(pass_idx).ref_y - pass(master_idx).ref_y).^2 );
    
  h_plot_baseline(end+1) = plot(array_param.fcs{1}{pass_idx}.base_line);
  h_legend_baseline{end+1} = sprintf('%d',pass_idx);
  hold on;

  array_param.fcs{1}{pass_idx}.surface = ref.surface;
end
xlabel('Range line');
ylabel('Baseline (m)');
grid on;
legend(h_plot_baseline,h_legend_baseline);

array_param.wfs.time = ref.wfs.time;
array_param.sv_fh = @array_proc_sv;
array_param.wfs.fc = pass(1).wfs.fc;
array_param.imgs = {[ones(length(pass),1), (1:length(pass)).']};
array_param.three_dim.en = true;


%%
array_param.method = 0;
[array_param0,result0] = array_proc(array_param,data);
array_param.method = 1;
[array_param1,result1] = array_proc(array_param,data);
array_param.method = 2;
[array_param2,result2] = array_proc(array_param,data);

figure(11); clf;
imagesc(lp(result0.val))
colormap(1-gray(256))
h_axes = gca;

figure(12); clf;
imagesc(lp(result1.val))
colormap(1-gray(256))
h_axes(end+1) = gca;

figure(13); clf;
imagesc(lp(result2.val))
colormap(1-gray(256))
h_axes(end+1) = gca;

figure(14); clf;
imagesc(lp(fir_dec(abs(data{1}(:,:,master_idx)).^2, ones(size(array_param.rline_rng)), 1, ...
  1-array_param.rline_rng(1), size(data{1}(:,:,master_idx),2)-length(array_param.rline_rng)+1)))
colormap(1-gray(256))
h_axes(end+1) = gca;

linkaxes(h_axes,'xy');

%% Save Results
[fn_dir,fn_name] = fileparts(fn);
fn_map = fullfile(fn_dir,[fn_name '_map.fig']);
saveas(h_fig_map,fn_map);
fn_elev = fullfile(fn_dir,[fn_name '_elev.fig']);
saveas(h_fig_elev,fn_elev);
fn_baseline = fullfile(fn_dir,[fn_name '_baseline.fig']);
saveas(h_fig_baseline,fn_baseline);
fn_standard = fullfile(fn_dir,[fn_name '_standard.fig']);
saveas(1,fn_standard);
fn_mvdr = fullfile(fn_dir,[fn_name '_mvdr.fig']);
saveas(2,fn_mvdr);
fn_music = fullfile(fn_dir,[fn_name '_music.fig']);
saveas(3,fn_music);
fn_single = fullfile(fn_dir,[fn_name '_single.fig']);
saveas(4,fn_single);

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
