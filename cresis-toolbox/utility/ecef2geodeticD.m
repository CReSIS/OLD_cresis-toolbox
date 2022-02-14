function [lat, lon, h] = ecef2geodeticD(x, y, z, ellipsoid)

% ecef2geodetic but output (lat and lon) are in degrees
% Author: Hara Madhav Talasila

[lat, lon, h] = ecef2geodetic(x, y, z, ellipsoid);
lat = lat*180/pi;
lon = lon*180/pi;

end