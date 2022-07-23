function d_final = distance_geodetic(varargin)

% function d_final = distance_geodetic(varargin)
% input (lat and lon) are in degrees, default elll = WGS84.ellipsoid
% d_final = distance_geodetic(lat1, lon1, lat2, lon2)
% d_final = distance_geodetic(lat1, lon1, lat2, lon2, elll)
% d_final = distance_geodetic(lat1, lon1, elev1, lat2, lon2, elev2)
% d_final = distance_geodetic(lat1, lon1, elev1, lat2, lon2, elev2, elll)
% Author: Hara Madhav Talasila

WGS84 = physical_constants('WGS84');

switch nargin
  case 4
    lat1 = varargin{1};
    lon1 = varargin{2};
    lat2 = varargin{3};
    lon2 = varargin{4};
    elll = WGS84.ellipsoid;
    elev1 = 0;
    elev2 = 0;
  case 5
    lat1 = varargin{1};
    lon1 = varargin{2};
    lat2 = varargin{3};
    lon2 = varargin{4};
    elll = varargin{5};
    elev1 = 0;
    elev2 = 0;
  case 6
    lat1 = varargin{1};
    lon1 = varargin{2};
    elev1 = varargin{3};
    lat2 = varargin{4};
    lon2 = varargin{5};
    elev2 = varargin{6};
    elll = WGS84.ellipsoid;
  case 7
    lat1 = varargin{1};
    lon1 = varargin{2};
    elev1 = varargin{3};
    lat2 = varargin{4};
    lon2 = varargin{5};
    elev2 = varargin{6};
    elll = varargin{7};
end

[A(1,:), A(2,:), A(3,:)] = geodeticD2ecef(lat1, lon1, elev1, elll);
[B(1,:), B(2,:), B(3,:)] = geodeticD2ecef(lat2, lon2, elev2, elll);

d_final = vecnorm(A-B);

end % EOF