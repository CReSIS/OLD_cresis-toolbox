function [WGS84elev_axis, surface_calc, AGL] = time2elev(t_axis, flight_elev, surf_idx)

% [dd.elev_axis, surface_calc(xx), AGL(xx)] = time2elev(dd.Time, dd.Elevation, surf_idx)
% create the elevation_axis
c = physical_constants('c');
er = physical_constants('er_ice');

Nt = length(t_axis);
air_idxs = 1:surf_idx;
AGL = t_axis(surf_idx) *c/2;
sub_idxs = min(surf_idx+1,Nt):Nt;
range_axis = [ t_axis(air_idxs) .*c/2 ; ...
    AGL+(t_axis(sub_idxs)-t_axis(surf_idx)) .*(c/2/sqrt(er))];
WGS84elev_axis = flight_elev - range_axis;
surface_calc = flight_elev - AGL;