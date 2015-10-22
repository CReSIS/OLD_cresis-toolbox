function normality = salinity_to_normality(salinity)
% normality = salinity_to_normality(salinity)
%
% Returns normality given a salinity according to eqn E.20 in Ulaby,
% Moore, Fung, Vol 3, Appendix E.
%
% salinity = salinity of the salt water solution
%
% Salinity is the total mass of solid salt in grams dissolved in one
% kilogram of solution. Usually measured in parts per thousand on a weight
% basis.
%
% Oceans have an average salinity of 32.54 parts per thousand.  So
% salinity = 0.03254.
%   
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m

A = 0.9141;
normality = A*salinity.*(1.707e-2 + 1.205e-5*salinity + 4.058e-9*salinity.^2);

return;
