function plot_lidar(fn, record_fn, plot_color)
% plot_lidar(fn OR lidar, record_fn, plot_color)
%
% Function for testing the quality of the GPS data.
% 
% fn = Matlab file with standard LIDAR structure in it.
%   OR lidar struct from read_lidar_atm.m
% record_fn = optional records file to load (default is [])
% plot_color = If empty or left blank, all plots are cleared.
%   If specified (e.g. 'r', 'g', etc), some of the plots will be
%   added to rather than cleared.
%
% Example:
%
% % Single file example
% atm_fns = '/cresis/data1/NASA/2009_Antarctica_DC8_ATM/20091102/091102_203702_smooth_nadir5seg_50pt';
% lidar = read_lidar_atm(atm_fns,struct('year',2009,'month',11,'day',02));
% plot_lidar(lidar);
% 
% % Multiple file example
% atm_fns = get_filenames('/cresis/data1/NASA/2009_Antarctica_DC8_ATM/20091112/','091112','','');
% records_fn = '/cresis/scratch1/mdce/csarp_support/records/mcords/2009_Antarctica_DC8/records_20091112_seg1.mat';
% lidar = read_lidar_atm(atm_fns,struct('year',2009,'month',11,'day',12));
% plot_lidar(lidar,records_fn);
%
% Author: John Paden
% 
% See also read_gps_applanix.m, read_gps_reveal.m,
%   make_gps_2009_antarctica_DC8_DGPSwINS.m,
%   make_gps_2009_antarctica_DC8_GPS.m

if ~exist('record_fn','var')
  record_fn = [];
end

if ~exist('plot_color','var') || isempty(plot_color)
  plot_color = 'b';
  clear_figures = true;
else
  clear_figures = false;
end

if ischar(fn)
  lidar = load(fn);
else
  lidar = fn;
end

if ~isempty(record_fn)
  records = load(record_fn);
  if isfield(records,'records')
    records = records.records;
  end
  min_time = min(min(lidar.gps_time),min(records.gps_time));
else
  min_time = min(lidar.gps_time);
end

figure(1);
if clear_figures
  clf;
else
  hold on;
end
plot(lidar.lon, lidar.lat, [plot_color '.']);
if ~isempty(record_fn)
  hold on;
  plot(records.lon,records.lat, ['g','o']);
  hold off;
end
xlabel('Longitude, E (deg)');
ylabel('Latitude, N (deg)');

figure(2);
if clear_figures
  clf;
else
  hold on;
end
plot(lidar.gps_time-min_time,lidar.surface,[plot_color,'o']);
xlabel('Index');
ylabel('Elevation (m)');

if ~isempty(record_fn) && isfield(records,'surface')
  physical_constants;
  hold on;
  plot(records.gps_time-min_time, records.elev - records.surface*c/2, ['r','x']);
  hold off;
end

figure(3); clf;
plot(diff(lidar.gps_time));
xlabel('Index');
ylabel('Time step (sec)');

med_sample_rate = median(diff(lidar.gps_time));

skip_idxs = find(diff(lidar.gps_time) > med_sample_rate*1.5);
fprintf('Approx. number of skips: %d\n', length(skip_idxs))
figure(1);
hold on;
plot(lidar.lon(skip_idxs), lidar.lat(skip_idxs), 'ro');
hold off;

repeat_idxs = find(diff(lidar.gps_time) <= 0);
fprintf('Approx. number of repeats: %d\n', length(repeat_idxs))
figure(1);
hold on;
plot(lidar.lon(repeat_idxs), lidar.lat(repeat_idxs), 'go');
hold off;

fprintf('Start: %s GPS\n', datestr(epoch_to_datenum(lidar.gps_time(1))));
fprintf('Stop: %s GPS\n', datestr(epoch_to_datenum(lidar.gps_time(end))));

leap_sec = utc_leap_seconds(lidar.gps_time(1));
fprintf('Start: %s UTC\n', datestr(epoch_to_datenum(lidar.gps_time(1) - leap_sec)));
fprintf('Stop: %s UTC\n', datestr(epoch_to_datenum(lidar.gps_time(end) - leap_sec)));

% Get the GPS seconds of day to sync to radar
gps_time_datenum = epoch_to_datenum(lidar.gps_time);
[year month day hour minute sec] = datevec(gps_time_datenum);

GPS_sod = (day-day(1))*86400+hour*3600+minute*60+sec;  % GPS seconds of day

return;




