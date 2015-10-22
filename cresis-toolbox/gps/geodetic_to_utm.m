function [east,north,mstruct] = geodetic_to_utm(lat,lon)
% [east,north,mstruct] = geodetic_to_utm(lat,lon)
%
% Converts from geodetic latitude/longitude to UTM. Returns the
% mstruct required to convert back and forth (mstruct.zone
% gives the UTM zone).
%
% lat = size X array (deg)
% lon = size X array (deg)
%
% east = size X array (m)
% north = size X array (m)
%
% Authors: John Paden
%
% See also: minvtran.m, mfwdtran.m

mstruct = defaultm('utm');
mstruct.zone = utmzone(lat,lon);
mstruct = defaultm(mstruct);
% geoid = utmgeoid(mstruct.zone);
% mstruct.geoid = geoid(1,:);
physical_constants;
mstruct.geoid = WGS84.ellipsoid;

[east,north] = mfwdtran(mstruct,lat,lon);

return;
