function soil(mv)
%
% Moisture (volume fraction)
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m, salinity_to_volume_fraction.m

format compact; format long;
if (ispc)
  path(path,'P:\prism\radar\radarSimulator\');
else
  path(path,'/projects/prism/radar/radarSimulator/');
end

physicalConstants;

er = 3 - j*0.05;
er = 2.7 - j*0.022*2.7*1.5;
er = 2.6 - j*0.015*2.6*1.5;
er = 2.8 - j*0.035*2.8*1.5
f = 210e6;
w = 2*pi*f;

gamma = j*w*sqrt(u0*e0*er);
alpha = real(gamma);

penetrationDepth = 1./alpha

return;


