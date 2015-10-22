function ops_atm_param = atm_to_ops(atm_fn,records_fn)
% ops_atm_param = atm_to_ops(atm_fn,records_fn)
% Convert L2 ATM  to the database insert param structure.
%
% Input:
%   atm_fn: Absolute path/s to a atm_L2_fn (smooth_nadir3seg_50pt) file/s.
%   records_fn: Absolute path to a records file (that matches the atm file)
%
% Output:
%   ops_atm_param: structure array with fields
%     geometry.coordinates = double array of format ([lon lat elevation])
%     properties.gps_time = double array
%     properties.twtt = double array (converted from ATM surface elevation)
%     properties.type = cell of strings ('auto')
%     properties.quality = integer array (always 1 for ATM data)
%     properties.lyr_name = string ('atm')
%
% Author: Kyle W. Purdon
%
% see also read_lidar_atm.m, physical_constants.m

% LOAD ATM
fprintf('\tReading ATM files ...\n');
lidar = read_lidar_atm(atm_fn);
fprintf('\tSyncing Records and ATM data ...\n');

% LOAD RECORDS
records = load(records_fn);
if isfield(records,'records')
  % Adding this check to support old record (pre-cr1) file format
  records = records.records;
end

% INTERPOLATE LIDAR SURF/LAT/LON TO RECORDS GPS
records.lidar.lat = interp1(lidar.gps_time,lidar.lat,records.gps_time);
records.lidar.lon = interp1(lidar.gps_time,lidar.lon,records.gps_time);
records.lidar.surface = interp1(lidar.gps_time,lidar.surface,records.gps_time);

% CONVERT ELEVATION/SURFACE TO TWTT
physical_constants;
records.lidar.twtt = (records.elev-records.lidar.surface)/(c/2);

% GET NON-NAN INDEXES ON LIDAR SURFACE
good_idxs = ~isnan(records.lidar.surface);

% RE-SAMPLE BACK TO 40M ATM ALONG-TRACK
decim_idxs = get_equal_alongtrack_spacing_idxs(geodetic_to_along_track(records.lidar.lat(good_idxs),records.lidar.lon(good_idxs),records.elev(good_idxs)),40);

% FORMAT output
ops_atm_param.properties.gps_time = records.gps_time(decim_idxs);
ops_atm_param.geometry.coordinates = [records.lidar.lon(decim_idxs); records.lidar.lat(decim_idxs); records.elev(decim_idxs)]';
ops_atm_param.properties.twtt = records.lidar.twtt(decim_idxs);
ops_atm_param.properties.type = ones(1,length(records.gps_time(decim_idxs)))*2;
ops_atm_param.properties.lyr_name = 'atm';
ops_atm_param.properties.quality = ones(1,length(records.gps_time(decim_idxs)));

end