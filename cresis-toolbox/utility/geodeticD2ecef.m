function [x, y, z] = geodeticD2ecef(lat, lon, h, ellipsoid)

% geodetic2ecef but input (lat and lon) are in degrees
% Author: Hara Madhav Talasila

[x, y, z] = geodetic2ecef(lat*pi/180, lon*pi/180, h, ellipsoid);

end