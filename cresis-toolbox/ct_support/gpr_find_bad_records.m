% script gpr_find_bad_records
%
% GPR profiles usually contain stops, 90+ deg sharp turns, etc that may be
% undesirable for SAR processing. Use this script to help find the 
% bad data records associated with these maneuvers and remove them.

%% User Settings
% records_fn = '/cresis/snfs1/dataproducts/csarp_support/records/accum/2013_Antarctica_Ground/records_20140103_06.mat';
% geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/antarctica/Landsat-7/Antarctica_LIMA_480m.tif';
records_fn = '/cresis/snfs1/dataproducts/csarp_support/records/accum/2015_Greenland_Ground/records_20150501_02.mat';
geotiff_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/mzl7geo_90m_lzw.tif';

bad_vel_threshold = 0.5;
bad_heading_diff_threshold = 1;

%% Automated Section
records = load(records_fn);

proj = geotiffinfo(geotiff_fn);

[x,y] = projfwd(proj,records.lat,records.lon);

%% Fabricating a heading now (pulled from make_gps.m)
gps = records;
along_track = geodetic_to_along_track(gps.lat,gps.lon,gps.elev);
rlines = get_equal_alongtrack_spacing_idxs(along_track,0.75);
physical_constants;
est_heading = size(gps.heading);
clear origin heading east north;
for rline_idx = 1:length(rlines)
  rline = rlines(rline_idx);
  if rline_idx < length(rlines)
    rline_end = rlines(rline_idx+1);
  else
    rline_end = length(along_track);
  end
  [origin(1),origin(2),origin(3)] = geodetic2ecef(gps.lat(rline)/180*pi,gps.lon(rline)/180*pi,gps.elev(rline),WGS84.ellipsoid);
  [heading(1),heading(2),heading(3)] = geodetic2ecef(gps.lat(rline_end)/180*pi,gps.lon(rline_end)/180*pi,gps.elev(rline_end),WGS84.ellipsoid);
  heading = heading - origin;
  % Determine east vector
  [east(1) east(2) east(3)] = lv2ecef(1,0,0,gps.lat(rline)/180*pi,gps.lon(rline)/180*pi,gps.elev(rline),WGS84.ellipsoid);
  east = east - origin;
  % Determine north vector
  [north(1) north(2) north(3)] = lv2ecef(0,1,0,gps.lat(rline)/180*pi,gps.lon(rline)/180*pi,gps.elev(rline),WGS84.ellipsoid);
  north = north - origin;
  % Determine heading
  est_heading(rline:rline_end) = atan2(dot(east,heading),dot(north,heading));
end

% Slow velocity mask
vel = diff(along_track) ./ diff(records.gps_time);
vel_mask = vel < bad_vel_threshold;
vel_mask(end+1) = 0;

% Fast heading change mask
heading_diff = diff(est_heading);
heading_mask = abs(heading_diff) > bad_heading_diff_threshold;
heading_mask(end+1) = 0;

% Plot records
figure; clf;
h_plots = [];
h_plots(end+1) = plot(x,y);
hold on;
fprintf('\nRed is slow velocity records\n');
h_plots(end+1) = plot(x(vel_mask),y(vel_mask),'r.');
fprintf('\nGreen is fast heading change records\n');
h_plots(end+1) = plot(x(heading_mask),y(heading_mask),'g.');

% Manual tool for removing records
clip_vectors(h_plots);

fprintf('\nRemove bad records (press F1 in plot for help). Once done, type dbcont.\n');
keyboard

% Find the bad records
ydata = get(h_plots(1),'YData');
good_mask = ~isnan(ydata);

figure; clf;
plot(records.lon(good_mask),records.lat(good_mask));

fprintf('Check to make sure you are happy with the results. Once done, type dbcont.\n');
keyboard

gpr_remove_bad_records(records_fn,good_mask);

return
