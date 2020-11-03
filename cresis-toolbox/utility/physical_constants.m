function [varargout] = physical_constants(varargin)
% function [varargout] = physical_constants(varargin)
% Adds desired or all physical constants to caller's workspace
% Outputs desired constants to assigned arguments
%
% Usage:  
%   physical_constants;
%       Adds all constants to workspace. General use.
%   physical_constants('c');
%   physical_constants('c','er_ice');
%   physical_constants('WGS84');
%   physical_constants('wgs84');
%   physical_constants('gearthmassprod');
%       Adds only selected constant(s) to workspace
%       case-insensitive and uses default case used in this script
%   [light_vel] = physical_constants('c');
%   [light_vel, grav_const] = physical_constants('c','G');
%       Assigns according to the input and ouput arguments
% 
% Authors: John Paden, Hara Madhav Talasila
%
% See Also: proj_load_standard.m

%% Physical Constants
% =========================================================================

% Permittivity of free-space (F*m^-1), http://en.wikipedia.org/wiki/Permittivity
e0 = 8.8541878176e-12;
% Permeability of free-space (N*A^-2), http://en.wikipedia.org/wiki/Permeability_%28electromagnetism%29
u0 = 4e-7*pi;
% Speed of light in free-space
c = 1/sqrt(e0*u0);
%c = 2.997924580003452e+08;
% Boltzmann's constant (J * K^-1 * Hz^-1), http://en.wikipedia.org/wiki/Boltzmann_constant
BoltzmannConst = 1.380650524e-23;
% Earth Radius (m) from WGS-84 ellipsoid (larger: equatorial radius,
% smaller: polar radius)
earthRadius = 0.5*(6378137+6356752.31424518);
% Earth mass (kg), http://en.wikipedia.org/wiki/Mass_of_the_Earth
earthMass = 5.9742e24;
% Gravitational constant (N * m^2 * kg^-2), http://en.wikipedia.org/wiki/Gravity_constant
G = 6.6742e-11;
% Note: G*earthMass = 398600.5 km^3/s^2 according to satellite orbit book,
% which does match the values above.
GearthMassProd = 398600.5e9;

% WGS84 ellipsoid parameters [semimajor sqrt(e2)] for geodetic2ecef
WGS84.semimajor = 6378137;
WGS84.semiminor = 6356752.314245;
WGS84.flattening = 298.257223563;
WGS84.eccentricity = sqrt(0.00669437999013);
WGS84.ellipsoid = [WGS84.semimajor WGS84.eccentricity];

% GRS80 ellipsoid parameters [semimajor sqrt(e2)] for geodetic2ecef
GRS80.semimajor = 6378137;
GRS80.flattening = 298.25722210;
GRS80.semiminor = GRS80.semimajor - 1/GRS80.flattening * GRS80.semimajor;
GRS80.eccentricity = sqrt(1 - GRS80.semiminor^2/GRS80.semimajor^2);
GRS80.ellipsoid = [GRS80.semimajor GRS80.eccentricity];

% Basic ice properties
er_ice = 3.15;
% er_ice = 1;

%% Automated section
% =========================================================================

local_var_list = who;
if nargin == 0 % default mode as before
  for idx = 1:length(local_var_list)
    if ~any(strcmpi('varargin',local_var_list{idx}))
      % Adds all physical constants to caller's workspace
      assignin('caller',local_var_list{idx},eval(local_var_list{idx}));
    end
  end
else % Get only required constants
  for idx = 1:length(varargin)
    loc = find(strcmpi(varargin{idx},local_var_list));
    if loc && nargout==nargin % assigns physical constant to out variable
      varargout{idx} = eval(local_var_list{loc});
    else % adds physical constant to caller's workspace
      assignin('caller',local_var_list{loc},eval(local_var_list{loc}));
    end
  end
end

end
