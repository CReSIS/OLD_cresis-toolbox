function [x,y,mstruct] = geodetic_to_stereographic(lat,lon,param)
% [x,y,mstruct] = geodetic_to_stereographic(lat,lon,param)
%
% Converts from geodetic latitude/longitude to stereographic. Returns the
% mstruct required to convert back and forth.
%
% lat = size X array (deg)
% lon = size X array (deg)
%
% x = size X array (m)
% y = size X array (m)
%
% Authors: John Paden
%
% See also: minvtran.m, mfwdtran.m

mstruct = defaultm('ups');
if mean(lat(:)) > 0
  mstruct.zone = 'north';
else
  mstruct.zone = 'south';
end
mstruct = defaultm(mstruct);
mstruct.falsenorthing = 0;
mstruct.falseeasting = 0;
physical_constants;
mstruct.geoid = WGS84.ellipsoid;

[x,y] = mfwdtran(mstruct,lat,lon);

return;
