function altra = fcs_local(varargin)

% function altra = fcs_local(varargin)
% create local flight coordinate system from records
% Author: Hara Madhav Talasila

switch nargin
  case 1
    records = varargin{1};
  case 2
    records = varargin{1};
    target  = varargin{2};
  otherwise
    warning('find_multiple: argin not supported !!!');
    return;
end

debug_plots_en = 0;
[c, WGS84] = physical_constants('c', 'WGS84');

%% Create local FCS (Flight Coordinate System)

altra = [];
% 1. Compute along-track vector (geodetic_to_along_track)
[altra.along_track, altra.lat, altra.lon, altra.elev] = ...
  geodetic_to_along_track(records.lat, records.lon, records.elev);

% 2. Find ECEF of trajectory
[altra.x, altra.y, altra.z] = geodeticD2ecef(altra.lat, altra.lon, altra.elev, WGS84.ellipsoid);

% 3. For each point in the trajectory, find the XYZ unit vectors in the
% flight coordinate system.
% X, Y, Z are 3x(rec_len) matrices of (rec_len) unit vectors
X = [];
U = [];
N = [];
Z = [];
Y = [];

% X is along track
X = [diff(altra.x); diff(altra.y); diff(altra.z)]; % 3x(rec_len-1)
X = X./vecnorm(X);
% extend last unit vector to the last point % 3x(rec_len)
X = [X, X(:,end)];

% U is up vector
[U(1,:), U(2,:), U(3,:)] = enu2ecef( 0, 0, 1, altra.lat, altra.lon, altra.elev, WGS84.ellipsoid);
U = U - [altra.x; altra.y; altra.z];
U = U./vecnorm(U);

% N is normal to the plane with X, U and Zenith
N = cross(X, U, 1);
N = N./vecnorm(N);

% Z is Zenith
Z = cross(N, X);
Z = Z./vecnorm(Z);

% Y is left ==> X Y Z are right-handed coord system
Y = cross(Z, X);

if debug_plots_en
  figure;
  plot(vecnorm(Y),'x'); hold on;
end

Y = Y./vecnorm(Y);

if debug_plots_en
  plot(vecnorm(Y),'o');
end

altra.X = X;
altra.Y = Y;
altra.Z = Z;

altra.U = U;
altra.N = N;

clear X U N Z Y

