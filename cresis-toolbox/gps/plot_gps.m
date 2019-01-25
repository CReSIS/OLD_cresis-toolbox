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

figure(1);
if clear_figures
  clf;
else
  hold on;
end
plot(gps.lon, gps.lat, [plot_color '.']);
xlabel('Longitude, E (deg)');
ylabel('Latitude, N (deg)');
grid on;
set(1,'Name','Trajectory');

figure(2);
if clear_figures
  clf;
else
  hold on;
end
plot(GPS_sod,gps.elev,plot_color);
xlabel('GPS\_sod');
ylabel('Elevation (m)');
grid on;
set(2,'Name','Elevation');

figure(3); clf;
plot(diff(gps.gps_time));
xlabel('Index');
ylabel('Time step (sec)');
grid on;
set(3,'Name','Diff GPS_time');

med_sample_rate = median(diff(gps.gps_time));

skip_idxs = find(diff(gps.gps_time) > med_sample_rate*1.5);
fprintf('Approx. number of skips: %d\n', length(skip_idxs))
figure(1);
hold on;
plot(gps.lon(skip_idxs), gps.lat(skip_idxs), 'ro');
hold off;

repeat_idxs = find(diff(gps.gps_time) <= 0);
fprintf('Approx. number of repeats: %d\n', length(repeat_idxs))
figure(1);
hold on;
plot(gps.lon(repeat_idxs), gps.lat(repeat_idxs), 'go');
hold off;

figure(4);
if clear_figures
  clf;
else
  hold on;
end
plot(GPS_sod,gps.roll*180/pi,plot_color);
xlabel('GPS sec of day (sec)');
ylabel('Roll (deg)');
grid on;
set(4,'Name','Roll');

figure(5);
if clear_figures
  clf;
else
  hold on;
end
plot(GPS_sod,gps.pitch*180/pi,plot_color);
xlabel('GPS sec of day (sec)');
ylabel('Pitch (deg)');
grid on;
set(5,'Name','Pitch');

figure(6);
if clear_figures
  clf;
else
  hold on;
end
plot(GPS_sod,gps.heading*180/pi,[plot_color '.']);
xlabel('GPS sec of day (sec)');
ylabel('True heading (deg)');
grid on;
set(6,'Name','Heading');

figure(7);
if clear_figures
  clf;
else
  hold on;
end
along_track = geodetic_to_along_track(gps.lat,gps.lon,gps.elev);
plot(GPS_sod(2:end),diff(along_track) ./ diff(gps.gps_time),[plot_color '.']);
xlabel('GPS sec of day (sec)');
ylabel('Speed (m/s)');
grid on;
title('Should not be zero in flight, should be a very smooth function, need to filter if not')
set(7,'Name','Speed');

if plot_est_heading
  tic;
  % Estimated heading
  est_heading = zeros(size(gps.heading));
  % Only do every Mx records, because this is slow
  Mx = 100;
  for idx = 2:Mx:length(gps.heading)
    if ~mod(idx-2,1000)
      fprintf('Estimated heading idx %d of %d (%.1f sec)\n', idx, length(gps.heading), toc);
    end
    [north,east] = geodetic_to_utm(gps.lat(idx-1:idx),gps.lon(idx-1:idx));
    est_heading(idx) = atan2(diff(north), diff(east));
  end
  figure(6); hold on;
  plot(GPS_sod(2:Mx:end), est_heading(2:Mx:end)*180/pi, [plot_color 'x']);
end

return;




