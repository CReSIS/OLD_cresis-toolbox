function new_atm_lidar = prep_lidar_atm(atm_lidar,maxgap)
%
% atm_lidar = prep_atm_lidar(lidar)
%
% Prepares ATM lidar for interpolation onto other data (RDS)
%   (1) Removes NaNs in ATM data
%   (2) Finds large along track gaps in ATM data (Add's NANs)
%
%   atm_lidar: data structure from read_lidar_atm.m
%   maxgap = maximum gap (in meters) to prep for interp. (500m works) 
%   new_atm_lidar: lidar data structure with "prepped" data.
%
% Author: Kyle W. Purdon
%
% see also read_lidar_atm

% Remove NANs from data
good_lidar_idxs = ~isnan(atm_lidar.gps_time);
atm_lidar.gps_time = atm_lidar.gps_time(good_lidar_idxs);
atm_lidar.surface = atm_lidar.surface(good_lidar_idxs);
atm_lidar.lat = atm_lidar.lat(good_lidar_idxs);
atm_lidar.lon = atm_lidar.lon(good_lidar_idxs);
if ~isempty(isnan(atm_lidar.surface))
  good_lidar_idxs = ~isnan(atm_lidar.surface);
  atm_lidar.gps_time = atm_lidar.gps_time(good_lidar_idxs);
  atm_lidar.surface = atm_lidar.surface(good_lidar_idxs);
  atm_lidar.lat = atm_lidar.lat(good_lidar_idxs);
  atm_lidar.lon = atm_lidar.lon(good_lidar_idxs);
end

% Find Gaps Based on Along_Track Distance
atm_lidar.along_track = geodetic_to_along_track(atm_lidar.lat,atm_lidar.lon,atm_lidar.surface);
atm_lidar.lidar_gaps = (diff(atm_lidar.along_track) > maxgap);

% Prepare Data for Interp1 if gaps exist
if sum(atm_lidar.lidar_gaps) ~= 0
  % Find Start Indexes of Gaps
  atm_lidar.lidar_gap_start_idxs = strfind([0 atm_lidar.lidar_gaps 0],[0 1]);
  atm_lidar.lidar_gap_start_idxs = atm_lidar.lidar_gap_start_idxs-1;
  if atm_lidar.lidar_gap_start_idxs(1) == 0
    atm_lidar.lidar_gap_start_idxs(1) = 1;
  end
  
  % Find End Indexes of Gaps
  atm_lidar.lidar_gap_end_idxs = strfind([0 atm_lidar.lidar_gaps 0],[1 0]);
  if atm_lidar.lidar_gap_end_idxs(end) > length(atm_lidar.lidar_gaps)
    atm_lidar.lidar_gap_end_idxs(end) = length(atm_lidar.lidar_gaps);
  end
  
  % Add NaN's to lidar at Gap_IDXS (Start & End)
  atm_lidar.surface([atm_lidar.lidar_gap_start_idxs atm_lidar.lidar_gap_end_idxs]) = NaN;
end

% Populate New ATM Output
new_atm_lidar.gps_time = atm_lidar.gps_time;
new_atm_lidar.lat = atm_lidar.lat;
new_atm_lidar.lon = atm_lidar.lon;
new_atm_lidar.surface = atm_lidar.surface;
new_atm_lidar.track = atm_lidar.track;


end