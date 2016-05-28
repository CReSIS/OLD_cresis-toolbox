
fn = '/cresis/snfs1/dataproducts/ct_data/rds/2016_Greenland_G1XB/CSARP_insar/good_line.mat';

load(fn);
master_idx = 1;
rbins = 20:120;

h_fig_map = figure(100); clf;
hold on;
h_fig_elev = figure(101); clf;
hold on;

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
  h_plot = plot(pass(pass_idx).lon, pass(pass_idx).lat,'.');
  color = get(h_plot,'Color');
  h_text = text(pass(pass_idx).lon(1), pass(pass_idx).lat(1), sprintf('%d', pass_idx), 'Color', color);
%   plot(pass(pass_idx).lon(1), pass(pass_idx).lat(1), 'o');
  
  figure(h_fig_elev);
  if master_idx == pass_idx
    plot(pass(pass_idx).elev,'LineWidth',2);
  else
    plot(pass(pass_idx).elev);
  end
  
  pass(pass_idx).ecef = [];
  [pass(pass_idx).ecef(1,:),pass(pass_idx).ecef(2,:),pass(pass_idx).ecef(3,:)] ...
    = geodetic2ecef(pass(pass_idx).lat/180*pi, pass(pass_idx).lon/180*pi, pass(pass_idx).elev, WGS84.ellipsoid);
end
linkaxes(h_data_axes,'xy');

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

%%
h_fig_ref_idx = figure(102); clf;
hold on;

h_data_axes = [];
data = [];
for pass_idx = 1:length(pass)
    
  pass(pass_idx).ref_idx = zeros(1,size(pass(pass_idx).origin,2));
  last_idx = 0;
  for rline = 1:size(pass(pass_idx).ecef,2)
    offset = bsxfun(@minus, pass(pass_idx).ecef(:,rline), ref.ecef);
    dist = offset.'*pass(pass_idx).x(:,rline);
    [min_dist,min_idx] = min(abs(dist));
    %       if min_idx == last_idx
    %         keyboard
    %       end
    last_idx = min_idx;
    pass(pass_idx).ref_idx(rline) = min_idx;
    
    x_offset = offset(:,min_idx).'*ref.x(:,min_idx);
    pass(pass_idx).along_track(rline) = along_track(min_idx) + x_offset;
    pass(pass_idx).ref_y(rline) = offset(:,min_idx).'*ref.y(:,min_idx);
    pass(pass_idx).ref_z(rline) = offset(:,min_idx).'*ref.z(:,min_idx);

  end
  
  pass(pass_idx).along_track_slave = geodetic_to_along_track(pass(pass_idx).lat,pass(pass_idx).lon,pass(pass_idx).elev);;
  
  figure(h_fig_ref_idx);
  plot(pass(pass_idx).ref_idx)
  drawnow;

  % Resample images onto a common along-track axes
  Mx = 10;
  Nx = size(pass(pass_idx).data,2);
  data_oversample = interpft(pass(pass_idx).data.',Mx*Nx);
  along_track_oversample = interp1(0:Nx-1, ...
    pass(pass_idx).along_track, (0:Nx*Mx-1)/Mx,'linear','extrap');
  pass(pass_idx).ref_data = interp1(along_track_oversample, ...
    data_oversample, along_track,'linear','extrap').';
  
  % Motion compensation
  for rline = 1:size(pass(pass_idx).ecef,2)
    dt = pass(pass_idx).ref_z(rline)/(c/2);
    pass(pass_idx).ref_data(:,rline) = ifft(fft(pass(pass_idx).ref_data(:,rline)) ...
      .*exp(1i*2*pi*pass(end).wfs.freq*dt) );
  end
  data = cat(3,data,pass(pass_idx).ref_data);
  
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

%%
h_fig_angle = figure(104); clf;
complex_data = surf_data(6,:) .* conj(surf_data(1,:));
complex_data(lp(complex_data) < -96) = NaN;
plot(180/pi*angle(complex_data),'.')
xlabel('Range line');
ylabel('Angle (deg)');
grid on;
