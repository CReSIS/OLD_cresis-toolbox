function opsAtmData = atmToOps(atmFn,settings)
% opsAtmData = atmToOps(atmFn,recordsFn)
%
% 1. Converts L2 ATM to the OPS layer format.
% 2. Interpolates L2 ATM onto a fixed scale based on point paths in the OPS database.
% 3. Removes any large gaps in the L2 ATM (see also data_gaps_check.m)
% 4. Removes any duplicate points in the L2 ATM
%
% Input:
%   atmFn: Absolute path/s to a atm_L2_fn (smooth_nadir3seg_50pt) file/s.
%   atmParam = struct that gives absolute reference of what day the data wastaken
%     .year
%     .month
%     .day
%     .time_reference = 'gps' or 'utc'
%     .fixed_field_width = logical (default is false), uses fixed column widths to load the file
%   settings: structure with the following fields
%     .location = string ('arctic','antarctic')
%     .seasonName = string
%     .sysName = string ('rds','accum',...)
%
% Output:
%   opsAtmData = structure with fields:
%       properties.point_path_id = integer array
%       properties.user = string
%       properties.twtt = double array (converted from ATM surface elevation)
%       properties.type = integer arry 2:auto
%       properties.quality = integer array 1:good
%       properties.lyr_name = string ('atm')
%
% Author: Kyle W. Purdon
%
% see also read_lidar_atm.m

% LOAD THE ATM FILE
lidar = read_lidar_atm(atmFn);

% FIND DUPLICATES IN ATM DATA
[~,notDupIdxs] = unique(lidar.gps_time);
newGpsTime = nan(size(lidar.gps_time));
newGpsTime(notDupIdxs) = lidar.gps_time(notDupIdxs);
lidar.gps_time = newGpsTime;
clear newGpsTime;

% REMOVE NAN FROM ATM DATA
keepIdxs = (sum(combine(~isnan(lidar.gps_time),~isnan(lidar.surface)))>=1);
lidar.gps_time = lidar.gps_time(keepIdxs);
lidar.lat = lidar.lat(keepIdxs);
lidar.lon = lidar.lon(keepIdxs);
lidar.surface  = lidar.surface(keepIdxs);

% GET OPS PATH INFORMATION
opsCmd;
pathParam.properties.location = settings.location;
pathParam.properties.season = settings.seasonName;
pathParam.properties.start_gps_time = min(lidar.gps_time);
pathParam.properties.stop_gps_time = max(lidar.gps_time);
pathParam.properties.nativeGeom = true;
[~,pathData] = opsGetPath(settings.sysName,pathParam);

% FIND GAPS IN DATA
lidarAlongTrack = geodetic_to_along_track(lidar.lat,lidar.lon,lidar.surface);
pathAlongTrack = geodetic_to_along_track(pathData.properties.Y,pathData.properties.X,pathData.properties.elev);
dataGapIdxs = data_gaps_check_mex(pathAlongTrack,lidarAlongTrack,50,20);

% INTERPOLATE LIDAR ONTO OPS PATH AND STORE OUPTUTS
opsAtmData.properties.point_path_id = pathData.properties.id;
opsAtmData.properties.twtt = (pathData.properties.elev-interp1(lidar.gps_time,lidar.surface,pathData.properties.gps_time))/(299792458.0/2);
opsAtmData.properties.type = ones(size(pathData.properties.id))*2;
opsAtmData.properties.quality = ones(size(pathData.properties.id));
opsAtmData.properties.lyr_name = 'atm';

% REMOVE GAPS FROM OUTPUTS
opsAtmData.properties.point_path_id = opsAtmData.properties.point_path_id(~dataGapIdxs);
opsAtmData.properties.twtt = opsAtmData.properties.twtt(~dataGapIdxs);
opsAtmData.properties.type = opsAtmData.properties.type(~dataGapIdxs);
opsAtmData.properties.quality = opsAtmData.properties.quality(~dataGapIdxs);

clear lidar pathData lidarAlongTrack pathAlongTrack % CLEAN UP

end