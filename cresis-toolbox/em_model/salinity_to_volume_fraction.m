function CF = salinity_to_volume_fraction(temp,method)
% CF = salinity_to_volume_fraction(temp,method)
%
% temp = temperature in Celcius
% method = 1, 2, or 3 (default is 1)
%   See code for descriptions.
%
% Returns constant for relating sea ice salinity to volume fraction:
%   vf = CF * salinity
% where salinity is in parts per thousand.
%
% Equation relating vf to salinity of sea-ice, brine, and density of
% sea-ice and brine is from
%   Ulaby, Moore, Fung, Vol 3, Appendix E, pg 2047, Brine volume fraction.
%
% Equation relating salinity of brine to temperature is from
%   Ulaby, Moore, Fung, Vol 3, Appendix E, pg 2045, dielectric constant of
%   brine
%
% Equation relating salinity of sea ice and temperature to volume fraction
% of brine is from
%   Ulaby, Moore, Fung, Vol 3, Appendix E, pg 2048, Brine volume fraction.
%
% Equation relating rho_brine to brine salinity and temperature is from
%   http://en.wikibooks.org/wiki/Methods_Manual_for_Salt_Lake_Studies/Salinity/measuring_brine_density
%   If you wish to convert temperature corrected SG to total dissolved solids (TDS g/L), a curve has been derived from Baseggio's paper on the composition of seawater derived hyperaline brines (Baseggio, 1974)
%   TDS (g/L) = -91897 SG4 + 403869 SG3 - 663919 SG2 + 485355 SG - 133408
%
% Pure water is 1000 kg per m^3 at 4 C
% One Liter is 1 dm^3 or 1 kg so salinity 35 g/kg is 35 g/L when the solution is water at 4 C
%   and 35 g/kg is 0.035 kg/kg or 35 parts per thousand by mass
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m, salinity_to_volume_fraction.m

%% Input Checking
if ~exist('method','var')
  method = 1;
end

%% Compute conversion factor
if method == 1
  %% More accurate, temperature ranges
  % Really only valid for -0.5 C > temp > -22.5 C
  CF = zeros(size(temp));
  CF(temp > -2.06) = 1e-3*(-52.56./temp(temp > -2.06) - 2.28);
  CF(temp <= -2.06 & temp > -8.2) = 1e-3*(-45.917./temp(temp <= -2.06 & temp > -8.2) - 0.930);
  CF(temp <= -8.2) = 1e-3*(-43.795./temp(temp <= -8.2) - 1.189);
elseif method == 2
  %% Combined, less accurate equation
  % Really only valid for -0.5 C > temp > -22.5 C
  CF = 1e-3*(-49.185./temp - 0.532);
elseif method == 3
  %% From components
  % Really only valid for -2 C > temp > -43.2 C
  % Very slow to compute due to function numerical minimization
  
  rho_ice = 0.916; % grams/cm^3
  
  % Compute Sb in parts per thousand (g/L)
  Sb = zeros(size(temp));
  n = find(temp>-8.2);
  Sb(n) = 1.725 - 18.756*temp(n) - 0.3964*temp(n).^2;
  n = find(temp <= -8.2 & temp > -22.9);
  Sb(n) = 57.041 - 9.929*temp(n) - 0.16204*temp(n).^2 - 0.002396*temp(n).^3;
  n = find(temp <= -22.9 & temp > -36.8);
  Sb(n) = 242.94 + 1.5299*temp(n) + 0.0429*temp(n).^2;
  n = find(temp <= -36.8);
  Sb(n) = 508.18 + 14.535*temp(n) + 0.2018*temp(n).^2;
  
  % Density of sea water ranges from 1020 to 1029 kg/m^3
  rho_b = zeros(size(Sb));
  for n = 1:numel(Sb)
    %Sb_func = @(rho_b) -91897*rho_b.^4 + 403869*rho_b.^3 - 663919*rho_b.^2 + 485355*rho_b.^1 - 133408 - Sb;
    Sb_func = inline(sprintf('abs(-91897*rho_b.^4 + 403869*rho_b.^3 - 663919*rho_b.^2 + 485355*rho_b.^1 - 133408 - %f)', Sb(n)));
    rho_b(n) = fminsearch(Sb_func,1.03);
  end
  rho_b = rho_b - 0.001*temp/3;
  
  CF = 1./Sb .* rho_ice ./ rho_b;
end

return;

temp = linspace(-22.5,-0.5,101);
CF = salinity_to_volume_fraction(temp,1);
plot(-temp,CF);
hold on;
CF = salinity_to_volume_fraction(temp,2);
plot(-temp,CF,'r');
CF = salinity_to_volume_fraction(temp,3);
plot(-temp,CF,'g');
hold off;
xlabel('temperature (-C)');
ylabel('ratio brine volume fraction/salinity');


