function density = densityProfile(A,densityInit,T,z)
% density = densityProfile(A,densityInit,T,z)
%
% A = mean annual accumulation rate (m water equiv per year)
% densityInit = initial density (mg/m^3)
% T = mean annual temperature (K)
% z = depth (m)
%
% Herron and Langway, Journal of Glaciology, 1980.

densityIce = 0.917;
densityCritical = 0.55;

b = 0.5;

R = 8.314;
k0 = 11*exp(-10160./(R*T));
k1 = 575*exp(-21400./(R*T));

% Depth where density = densityCritical Mg/m^3
zCritical = 1./(densityIce.*k0) .* (log(densityCritical./(densityIce-densityCritical)) ...
  - log(densityInit./(densityIce-densityInit)));

density = zeros(size(z));

ind = find(z <= zCritical);
Z0 = exp(densityIce.*k0.*z(ind) + log(densityInit./(densityIce-densityInit)));
density(ind) = densityIce.*Z0 ./ (1+Z0);

ind = find(z > zCritical);
Z1 = exp(densityIce.*k1.*(z(ind)-zCritical)./A.^b ...
  + log(densityCritical./(densityIce-densityCritical)));
density(ind) = densityIce.*Z1 ./ (1+Z1);

return;

% ------------------------------------------------------------------------------
% Examples
% ------------------------------------------------------------------------------
z = 0:70;
density = densityProfile(0.28,0.3/100,-30+273.15,z);
plot(density,z);
set(get(gcf,'CurrentAxes'),'ydir','reverse');


z = 0:70;
density = densityProfile(0.5,0.35,-21.5+273.15,z);
plot(density,z);
set(get(gcf,'CurrentAxes'),'ydir','reverse');



