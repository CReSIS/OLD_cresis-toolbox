function er = rock(T,freq,type)
% er = rock(T,freq)
%
% T = temperature in Kelvin
% freq = frequency in Hz
% type = type of rock (choose from table)
%
% e0*er = e0*(e' - j*e'')
% 
% Type  Rock            Ref
%    1  Granite         1
%    2  Granite         1
%    3  Granite         1
%    4  Granite         1
%    5  Sandstone       1
%    6  Sandstone       1
%    7  Limestone       1
%    8  Gneiss          1
%    9  Granite-Gneiss  1
% 
% References
% 1. Parkhomenko, E. I., "Electrical Properties of Rocks," Translated and Edited by Keller,
%    Plenum Press, New York, 1967, pp. 38-41.
%
% ------------------------------------------------------------------------------
% REFERENCES FROM CHRIS ALEN
% Granite
% Campbell, M.J. and J. Ulrichs, "Electrical properties of rocks and their significance for lunar radar
% observations," Journal of Geophysical Research, 74, pp. 5867-5881, 1969.
% 
% Granite & Limestone
%  Olhoeft, G.R., "Electrical properties of rocks," in Physical Properties of  Rocks and Minerals,
%  Touloukian, Y. S., Judd, W. R., and Roy, R. F., eds.: New York, McGraw-Hill, p. 257-330, 1981.
% 
% Clay
% Wakita, Y. and Y. Yamaguchi, "Estimation of the soil permittivity and conductivity by a GPR antenna,"
% Proceedings of the 6th International Conference on Ground Penetrating Radar (GPR '96), Sendai, Japan, pp.
% 123-126, 1996.
% 
% Ulaby, F.T., R.K. Moore, and A.K. Fung, Microwave Remote Sensing: Active and Passive, Volume III: From
% theory to applications, Artech House,
% Dedham, Massachusetts, 1986.
% 
% Dobson, M.C., F.T. Ulaby, M.T. Hallikainen, and M.A. El-Rayes, "Microwave dielectric behavior of wet
% soil-part II: dielectric mixing models," IEEE
% Transactions on Geoscience and Remote Sensing, 23(1), pp. 35-46, 1985.
% 
% Peplinski, N.R., F.T. Ulaby, M.C. Dobson, "Dielectric properties of soils in the 0.3-1.3-GHz range," IEEE
% Transactions on Geoscience and Remote
% Sensing, 33(3), pp. 803-807, 1995.
% ------------------------------------------------------------------------------
%
% Look at Carl Leuschen's PhD
%
% Look at Matthew Peters, 2005, JGR, Table 1 (very good table/refs for bedrock properties)
%
% Author: John Paden
%
% See also brine.m, sea_ice.m, rock.m, salt_water.m, soil.m, water.m
%   and salinity_to_normality.m, salinity_to_volume_fraction.m

% Eventually:
% type,1DCELL{name1,name2,etc},2DMATRIX{frequency,dielectric,conductivity,lossTan,moisture content,temperature}
% type,name,dielectric,frequency,moisture content
erRock = {1  'Granite'          5.42  5e5   0
          2  'Granite'          4.74  5e5   0
          3  'Granite'          5.06  5e5   0
          4  'Granite'          4.5   5e5   0
          5  'Sandstone'        4.66  5e5   0
          6  'Sandstone'        3.96  5e5   0
          7  'Limestone'        7.3   NaN   0
          8  'Gneiss'          11.5   NaN  NaN
          9  'Granite-Gneiss'   8.5   5e7  NaN};

er = erRock{type,3};

return;

