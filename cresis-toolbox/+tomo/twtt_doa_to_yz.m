function [y,z] = twtt_doa_to_yz(time,doa,theta,surface,er_ice)
% [y,z] = twtt_doa_to_yz(time,doa,theta,surface,er_ice)
%
% Converts twtt,doa coordinates (the standard synthetic aperture radar
% processing coordinate system) to Cartesian coordinates.
%
% theta and doa: +90 is y/left, 0 is -z/down, -90 is -y/right
% Nt: number of range (fast-time) bins
% Nx: number of range lines (along-track)/records
% Nsig: number of targets/sources/signals per fast-time/range bin
% Nsv: number of steering vectors/directions of arrival
%
% time: Nt by 1 [seconds]
%   Fast time vector to targets (used for each of the Nx positions since
%   the time vector to each target is not expected to change in that
%   dimension). The time vector is the two way travel time to the target
% doa: Nsig by Nx array [radians]
%   Direction of arrival to targets
% theta: Nsv by 1 vector [radians]
%   Direction of arrival vector corresponding to the surface matrix
% surface: Nsv by Nx array or 1 by Nx [seconds]
%   The two way travel time to the ice surface.
% er_ice: scalar containing relative dielectric of ice
%
% y,z: Nt by Nsig by Nx matrix [meters]
%   Target positions in Cartesian coordinates (y is left, z is up)
%
% mdata = load(DATAFILE);
% theta = reshape(mdata.param_combine.array_param.theta,[Nsig 1]);
% Nsig = length(theta); Nx = length(mdata.GPS_time);
% [y,z] = tomo.twtt_doa_to_yz(mdata.Time,repmat(theta,[1 Nx]),theta,mdata.twtt,3.15);
%
% Author: John Paden

% Load in standard constants like speed of light, c
physical_constants;

% For each pixel in the tomographic image in the twtt,doa coordinate
% system, we will convert the position to y,z accounting for refraction at
% the ice surface.
Nx = size(doa,2);
Nsig = size(doa,1);
Nt = size(time,1);
Nsv = length(theta);
y = NaN*zeros(Nt,Nsig,Nx);
z = NaN*zeros(Nt,Nsig,Nx);
warning('off','MATLAB:interp1:NaNinY');
for rline = 1:Nx
  rline
  
  % Convert surface for this along-track slope to y,z coordinates
  surf_y = sin(theta) .* surface(:,rline) * c/2;
  surf_z = -cos(theta) .* surface(:,rline) * c/2;
   
  for doa_idx = 1:Nsig
    upper_bound = find(theta > doa(doa_idx,rline),1);
    if ~isempty(upper_bound) & doa(doa_idx,rline) > theta(1)
      % Only update y,z if it falls within the theta support (i.e. we only
      % interpolate and do not extrapolate)
      
      % Determine the time to the surface and the position of incidence
      twtt_to_surf = interp1(theta, surface(:,rline), doa(doa_idx,rline));
      inc_y = sin(doa(doa_idx,rline)) * twtt_to_surf * c/2;
      inc_z = -cos(doa(doa_idx,rline)) * twtt_to_surf * c/2;
      
      % For each time position up to the surface, finding y,z is a cylindrical to
      % cartesian coordinate conversion.
      air_mask = time < twtt_to_surf;
      y(air_mask,doa_idx,rline) = sin(doa(doa_idx,rline)) * time(air_mask) * c/2;
      z(air_mask,doa_idx,rline) = -cos(doa(doa_idx,rline)) * time(air_mask) * c/2;
      
      % Interpolation indices
      interp_idxs = max(1,upper_bound-3) : min(Nsv,upper_bound+2);
      interp_idxs = interp_idxs(~isnan(surf_y(interp_idxs)) & ~isnan(surf_z(interp_idxs)));
      
      % Ensure there are still enough points to polyfit
      if length(interp_idxs) < 2
        continue;
      end
      
      % Determine the surface slope vector (ignoring along-track slope)
      p = polyfit(surf_y(interp_idxs),surf_z(interp_idxs),1);
      theta_slope = atan(p(1));
      theta_inc = doa(doa_idx,rline) + theta_slope;
      if abs(theta_inc) >= pi/2
        continue;
      end
      theta_tx = asin(sin(theta_inc)/sqrt(er_ice));
      
      ice_mask = time > twtt_to_surf;
      y(ice_mask,doa_idx,rline) = inc_y + sin(doa(doa_idx,rline)) ...
        * (time(ice_mask)-twtt_to_surf) * c/2;
      z(ice_mask,doa_idx,rline) = inc_z - cos(doa(doa_idx,rline)) ...
        * (time(ice_mask)-twtt_to_surf) * c/2;
      
    end
    
  end
  
end
warning('on','MATLAB:interp1:NaNinY');

return;

