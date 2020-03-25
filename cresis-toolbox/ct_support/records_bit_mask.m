function records_bit_mask(param,param_override)
%
% GPR profiles usually contain stops, 90+ deg sharp turns, etc that may be
% undesirable for SAR processing. Use this script to help find the
% bad data records associated with these maneuvers and remove them.
%
% param = struct with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Author: Nick Holschuh, John Paden
%
% See also: run_records_bit_mask.m, records_bit_mask.m

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Input Checks: records_bit_mask
% =====================================================================

% Here we set threshold values for the velocity and heading
if ~isfield(param.records_bit_mask,'bad_vel_threshold') || isempty(param.records_bit_mask.bad_vel_threshold)
  param.records_bit_mask.bad_vel_threshold = 0.25;
end
bad_vel_threshold = param.records_bit_mask.bad_vel_threshold;

if ~isfield(param.records_bit_mask,'bad_heading_diff_threshold') || isempty(param.records_bit_mask.bad_heading_diff_threshold)
  param.records_bit_mask.bad_heading_diff_threshold = 1;
end
bad_heading_diff_threshold = param.records_bit_mask.bad_heading_diff_threshold;

if ~isfield(param.records_bit_mask,'manual_masking') || isempty(param.records_bit_mask.manual_masking)
  param.records_bit_mask.manual_masking = 1;
end

if ~isfield(param.records_bit_mask,'debug_plot') || isempty(param.records_bit_mask.debug_plot)
  param.records_bit_mask.debug_plot = 1;
end

%% Load records and background geotiff for GUI

geotiff_fn = ct_filename_gis(param,param.records.frames.geotiff_fn);

records = records_load(param);

proj = geotiffinfo(geotiff_fn);

[x,y] = projfwd(proj,records.lat,records.lon);

%% Fabricating a heading now (pulled from gps_make.m)
gps = records;
along_track = geodetic_to_along_track(gps.lat,gps.lon);
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

% This either initiates the manual removal of records or does it
% automatically based on thresholds.
if param.records_bit_mask.manual_masking == 1
  % Plot records
  figure; clf;
  h_plots = [];
  h_plots(end+1) = plot(x,y);
  hold on;
  fprintf('\nRed is slow velocity records\n');
  h_plots(end+1) = plot(x(vel_mask),y(vel_mask),'r.');
  fprintf('\nGreen is fast heading change records\n');
  if ~any(heading_mask)
    h_plots(end+1) = plot(NaN,NaN,'g.');
  else
    h_plots(end+1) = plot(x(heading_mask),y(heading_mask),'g.');
  end
  % Manual tool for removing records
  clip_vectors(h_plots);
  
  fprintf('\nRemove bad records (press F1 in plot for help). Once done, type dbcont.\n');
  keyboard
  
  % Find the bad records
  ydata = get(h_plots(1),'YData');
  good_mask = ~isnan(ydata);
  
else
  good_mask = [~heading_mask & ~vel_mask];
end


if param.records_bit_mask.debug_plot == 1
  figure; clf;
  plot(records.lon(good_mask),records.lat(good_mask));
  title([num2str(sum(good_mask)),' of ',num2str(length(good_mask)),' remaining']);
  fprintf('Check to make sure you are happy with the results. Once done, type dbcont.\n');
  keyboard
end

%%%%%%%% Here we use the good_mask info to produce records.bit_mask.
%%%%%%%% This preserves data already masked for bad headers, but adds masks
%%%%%%%% to stationary or incorrectly moving records.

if ~isfield(records,'bit_mask')
  records.bit_mask = zeros(size(records.offset));
end

for board_idx = 1:size(records.offset,1)
  % For records that are good and are considered bad by this analysis, set
  % them as bad.
  mask = [records.bit_mask(board_idx,:) == 0 & ~good_mask];
  records.bit_mask(board_idx,mask) = 2;
  
  % For records that are bad and are now considered good by this analysis,
  % set them as good.
  mask = [records.bit_mask(board_idx,:) == 2 & good_mask];
  records.bit_mask(board_idx,mask) = 0;
end

records_fn = ct_filename_support(param,'','records');
fprintf('  Saving records %s\n', records_fn);
ct_save(records_fn,'-struct','records');


