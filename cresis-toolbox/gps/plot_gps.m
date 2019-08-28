function gps = plot_gps(fn, plot_color, plot_est_heading)
% gps = plot_gps(fn OR gps, plot_color, plot_est_heading)
%
% Function for testing the quality of the GPS data.
%
% fn = Matlab file with standard GPS structure in it.
%   OR gps struct from read_gps_reveal, read_gps_applanix, etc.
% plot_color = If empty or left blank, all plots are cleared.
%   If specified (e.g. 'r', 'g', etc), some of the plots will be
%   added to rather than cleared.
% plot_est_heading = Plot estimated heading, takes a long time
%   don't pass or leave empty for default which is false
%
% Example:
%
%   % Compare two GPS sources:
%   fn = '/cresis/scratch1/mdce/csarp_support/gps/2009_Antarctica_DC8_GPS/20091018_ALL_pos.mat';
%   plot_gps(fn);
%   fn = '/cresis/scratch1/mdce/csarp_support/gps/2009_Antarctica_DC8_DGPSwINS/20091018_ALL_pos.mat';
%   plot_gps(fn,'r');
%
%   % Use GPS struct
%   fn = '/cresis/data2/MCoRDS/2009_Chile/GPS/IWG1_110209.log';
%   gps = read_gps_reveal(fn);
%   plot_gps(gps);
%
%   fn = '/cresis/scratch1/mdce/csarp_support/gps/2009_Antarctica_DC8_DGPSwINS/Javadsbet_02Nov09_PPP_Revised.out'
%   gps = read_gps_applanix(fn,struct('year',2009,'month',11,'day',2));
%   plot_gps(gps);
%
% Author: John Paden
%
% See also read_gps_applanix.m, read_gps_reveal.m,
%   make_gps_2009_antarctica_DC8_DGPSwINS.m,
%   make_gps_2009_antarctica_DC8_GPS.m

if ~exist('plot_est_heading','var') || isempty(plot_est_heading)
  plot_est_heading = false;
end

if ~exist('plot_color','var') || isempty(plot_color)
  plot_color = 'b';
  clear_figures = true;
else
  clear_figures = false;
end

if ischar(fn)
  gps = load(fn);
else
  gps = fn;
end

fprintf('Start: %s GPS\n', datestr(epoch_to_datenum(gps.gps_time(1))));
fprintf('Stop: %s GPS\n', datestr(epoch_to_datenum(gps.gps_time(end))));

leap_sec = utc_leap_seconds(gps.gps_time(1));
fprintf('Start: %s UTC\n', datestr(epoch_to_datenum(gps.gps_time(1) - leap_sec)));
fprintf('Stop: %s UTC\n', datestr(epoch_to_datenum(gps.gps_time(end) - leap_sec)));

% Get the GPS seconds of day to sync to radar
gps_time_datenum = epoch_to_datenum(gps.gps_time);
[year month day hour minute sec] = datevec(gps_time_datenum);

GPS_sod = (day-day(1))*86400+hour*3600+minute*60+sec;  % GPS seconds of day

fig = 7;
if isfield(gps,'sync_gps_time')
  fig = fig+1;
end
if isfield(gps,'radar_time')
  fig = fig+1;
end
h_fig = get_figures(fig,true);
fig = 0;
h_axes = [];

%% Longitude vs Latitude Plot
fig = fig+1;
figure(h_fig(fig));
if clear_figures
  clf(h_fig(fig));
  h_axes(fig) = axes;
else
  tmp = get(h_fig(fig),'Children');
  if isa(tmp,'matlab.graphics.axis.Axes')
    h_axes(fig) = tmp(1);
  else
    h_axes(fig) = axes;
  end
end
plot(h_axes(fig), gps.lon, gps.lat, [plot_color '.']);
hold(h_axes(fig),'on');
grid(h_axes(fig),'on');
xlabel(h_axes(fig), 'Longitude, E (deg)');
ylabel(h_axes(fig), 'Latitude, N (deg)');
set(h_fig(fig),'Name','Trajectory');

med_sample_rate = median(diff(gps.gps_time));

skip_idxs = find(diff(gps.gps_time) > med_sample_rate*1.5);
fprintf('Approx. number of skips in gps time: %d\n', length(skip_idxs))
plot(h_axes(fig),gps.lon(skip_idxs), gps.lat(skip_idxs), 'ro');

repeat_idxs = find(diff(gps.gps_time) <= 0);
fprintf('Approx. number of repeats (non-monotonically increasing points in gps time):\n  %d\n', length(repeat_idxs))
plot(h_axes(fig),gps.lon(repeat_idxs), gps.lat(repeat_idxs), 'go');

%% Elevation Plot
fig = fig+1;
figure(h_fig(fig));
if clear_figures
  clf(h_fig(fig));
  h_axes(fig) = axes;
else
  tmp = get(h_fig(fig),'Children');
  if isa(tmp,'matlab.graphics.axis.Axes')
    h_axes(fig) = tmp(1);
  else
    h_axes(fig) = axes;
  end
end
plot(h_axes(fig),GPS_sod,gps.elev,plot_color);
hold(h_axes(fig),'on');
grid(h_axes(fig),'on');
xlabel(h_axes(fig),'GPS\_sod');
ylabel(h_axes(fig),'Elevation (m)');
set(h_fig(fig),'Name','Elevation');

%% GPS Time Plot
fig = fig+1;
figure(h_fig(fig));
if clear_figures
  clf(h_fig(fig));
  h_axes(fig) = axes;
else
  tmp = get(h_fig(fig),'Children');
  if isa(tmp,'matlab.graphics.axis.Axes')
    h_axes(fig) = tmp(1);
  else
    h_axes(fig) = axes;
  end
end
plot(h_axes(fig),diff(gps.gps_time));
hold(h_axes(fig),'on');
grid(h_axes(fig),'on');
xlabel(h_axes(fig),'Index');
ylabel(h_axes(fig),'Time step (sec)');
set(h_fig(fig),'Name','Diff GPS_time');

%% Roll Plot
fig = fig+1;
figure(h_fig(fig));
if clear_figures
  clf(h_fig(fig));
  h_axes(fig) = axes;
else
  tmp = get(h_fig(fig),'Children');
  if isa(tmp,'matlab.graphics.axis.Axes')
    h_axes(fig) = tmp(1);
  else
    h_axes(fig) = axes;
  end
end
plot(h_axes(fig),GPS_sod,gps.roll*180/pi,plot_color);
hold(h_axes(fig),'on');
grid(h_axes(fig),'on');
xlabel(h_axes(fig),'GPS sec of day (sec)');
ylabel(h_axes(fig),'Roll (deg)');
set(h_fig(fig),'Name','Roll');

%% Pitch Plot
fig = fig+1;
figure(h_fig(fig));
if clear_figures
  clf(h_fig(fig));
  h_axes(fig) = axes;
else
  tmp = get(h_fig(fig),'Children');
  if isa(tmp,'matlab.graphics.axis.Axes')
    h_axes(fig) = tmp(1);
  else
    h_axes(fig) = axes;
  end
end
plot(h_axes(fig),GPS_sod,gps.pitch*180/pi,plot_color);
hold(h_axes(fig),'on');
grid(h_axes(fig),'on');
xlabel(h_axes(fig),'GPS sec of day (sec)');
ylabel(h_axes(fig),'Pitch (deg)');
set(h_fig(fig),'Name','Pitch');

%% Heading Plot
fig = fig+1;
figure(h_fig(fig));
if clear_figures
  clf(h_fig(fig));
  h_axes(fig) = axes;
else
  tmp = get(h_fig(fig),'Children');
  if isa(tmp,'matlab.graphics.axis.Axes')
    h_axes(fig) = tmp(1);
  else
    h_axes(fig) = axes;
  end
end
plot(h_axes(fig),GPS_sod,gps.heading*180/pi,[plot_color '.']);
hold(h_axes(fig),'on');
grid(h_axes(fig),'on');
xlabel(h_axes(fig),'GPS sec of day (sec)');
ylabel(h_axes(fig),'True heading (deg)');
set(h_fig(fig),'Name','Heading');

if plot_est_heading
  % Estimated heading
  est_heading = zeros(size(gps.heading));
  % Only do every Mx records, because this is slow
  Mx = 100;
  for idx = 2:Mx:length(gps.heading)
    if ~mod(idx-2,1000)
      fprintf('Estimated heading idx %d of %d (%s sec)\n', idx, length(gps.heading), datestr(now));
    end
    [north,east] = geodetic_to_utm(gps.lat(idx-1:idx),gps.lon(idx-1:idx));
    est_heading(idx) = atan2(diff(north), diff(east));
  end
  plot(h_axes(fig),GPS_sod(2:Mx:end), est_heading(2:Mx:end)*180/pi, [plot_color 'x']);
end

%% Speed Plot
fig = fig+1;
figure(h_fig(fig));
if clear_figures
  clf(h_fig(fig));
  h_axes(fig) = axes;
else
  tmp = get(h_fig(fig),'Children');
  if isa(tmp,'matlab.graphics.axis.Axes')
    h_axes(fig) = tmp(1);
  else
    h_axes(fig) = axes;
  end
end
along_track = geodetic_to_along_track(gps.lat,gps.lon,gps.elev);
plot(h_axes(fig),GPS_sod(2:end),diff(along_track) ./ diff(gps.gps_time),[plot_color '.']);
hold(h_axes(fig),'on');
grid(h_axes(fig),'on');
xlabel(h_axes(fig),'GPS sec of day (sec)');
ylabel(h_axes(fig),'Speed (m/s)');
title(h_axes(fig),'Should not be zero in flight, should be a very smooth function, need to filter if not')
set(h_fig(fig),'Name','Speed');

%% Sync GPS Time Plot
if isfield(gps,'sync_gps_time')
  fig = fig+1;
  figure(h_fig(fig));
  if clear_figures
    clf(h_fig(fig));
    h_axes(fig) = axes;
  else
    tmp = get(h_fig(fig),'Children');
    if isa(tmp,'matlab.graphics.axis.Axes')
      h_axes(fig) = tmp(1);
    else
      h_axes(fig) = axes;
    end
  end
  plot(h_axes(fig),diff(gps.sync_gps_time));
  hold(h_axes(fig),'on');
  grid(h_axes(fig),'on');
  xlabel(h_axes(fig),'Index');
  ylabel(h_axes(fig),'Time step (sec)');
  set(h_fig(fig),'Name','Diff Sync GPS_time');
  
  repeat_idxs = find(diff(gps.sync_gps_time) <= 0);
  fprintf('Approx. number of repeats (non-monotonically increasing points in sync gps time):\n  %d\n', length(repeat_idxs))
end

%% Sync GPS Time Plot
if isfield(gps,'radar_time')
  fig = fig+1;
  figure(h_fig(fig));
  if clear_figures
    clf(h_fig(fig));
    h_axes(fig) = axes;
  else
    tmp = get(h_fig(fig),'Children');
    if isa(tmp,'matlab.graphics.axis.Axes')
      h_axes(fig) = tmp(1);
    else
      h_axes(fig) = axes;
    end
  end
  plot(h_axes(fig),diff(gps.radar_time));
  hold(h_axes(fig),'on');
  grid(h_axes(fig),'on');
  xlabel(h_axes(fig),'Index');
  ylabel(h_axes(fig),'Time step (sec)');
  set(h_fig(fig),'Name','Diff Radar Time');
  
  repeat_idxs = find(diff(gps.radar_time) <= 0);
  fprintf('Approx. number of repeats (non-monotonically increasing points in radar time):\n  %d\n', length(repeat_idxs))
end
