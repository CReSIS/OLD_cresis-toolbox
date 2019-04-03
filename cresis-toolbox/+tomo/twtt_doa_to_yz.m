function [y,z] = twtt_doa_to_yz(doa,theta,surface,er_ice,twtt,doa_method_flag,doa_limits)
% [y,z] = tomo.twtt_doa_to_yz(doa,theta,surface,er_ice,twtt)
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
% doa: Nsig by Nx array [radians]
%   Direction of arrival to targets
% theta: Nsv by 1 vector [radians]
%   Direction of arrival vector corresponding to the surface matrix
% surface: Nsv by Nx array or 1 by Nx [seconds]
%   The two way travel time to the ice surface.
% er_ice: scalar containing relative dielectric of ice
% twtt: Nsig by Nx array [seconds]
%   The two way travel to the target
%
% y,z: Nt by Nsig by Nx matrix [meters]
%   Target positions in Cartesian coordinates (y is left, z is up)
% y,z: Nsv by Nx matrix [meters]
%   Target positions in Cartesian coordinates (y is left, z is up)
%
% mdata = load(DATAFILE);
% Nx = length(mdata.GPS_time);
% theta = mdata.param_combine.array_param.theta; theta = reshape(theta,[Nsig 1]);
% [y,z] = tomo.twtt_doa_to_yz(mdata.Time,repmat(theta,[1 Nx]),theta,mdata.twtt,3.15);
%
% Author: John Paden

% Load in standard constants like speed of light, c
physical_constants;

% For each pixel in the tomographic image in the twtt,doa coordinate
% system, we will convert the position to y,z accounting for refraction at
% the ice surface.
if ~exist('doa_method_flag','var')
  doa_method_flag = false;
end
Nx = size(twtt,2);
if ~doa_method_flag
  Nsig = size(twtt,1);
  Nsv = length(theta);
end
y = NaN*zeros(size(twtt));
z = NaN*zeros(size(twtt));
warning('off','MATLAB:interp1:NaNinY');
if ~doa_method_flag
  %% Beamforming method
  for rline = 1:Nx
    % Convert surface for this along-track slope to y,z coordinates
    surf_y = sin(theta) .* surface(:,rline) * c/2;
    surf_z = -cos(theta) .* surface(:,rline) * c/2;
    p = polyfit(surf_y,surf_z,1);
    
    % Determine the time to the surface and the position of incidence
    twtt_to_surf = interp1(theta, surface(:,rline), doa(:,rline));
    for doa_idx = 1:Nsig
      inc_y = sin(doa(doa_idx,rline)) * twtt_to_surf(doa_idx) * c/2;
      inc_z = -cos(doa(doa_idx,rline)) * twtt_to_surf(doa_idx) * c/2;
      
      % For each time position up to the surface, finding y,z is a cylindrical to
      % cartesian coordinate conversion.
      if twtt(doa_idx,rline) <= twtt_to_surf(doa_idx)
        y(doa_idx,rline) = sin(doa(doa_idx,rline)) * twtt(doa_idx,rline) * c/2;
        z(doa_idx,rline) = -cos(doa(doa_idx,rline)) * twtt(doa_idx,rline) * c/2;
        
      elseif twtt(doa_idx,rline) > twtt_to_surf(doa_idx)
        
        % Interpolation indices for surface
        upper_bound = find(theta > doa(doa_idx,rline),1);
        if isempty(upper_bound)
          upper_bound = Nsv;
        end
        
        % IDEAL: Determine the local surface slope vector... surface needs to be smooth or the gradient/slope will be noisy
        %  (ignoring along-track slope since SAR focussing should accomodate for)
        % NOT USING IDEAL: surface slope approximated by single polynomial
        % for each range line
        % theta_slope = atan(p(1));
        % NOT USING IDEAL: surface slope approximated by zero
        
        theta_slope = 0;
        theta_inc = doa(doa_idx,rline) + theta_slope;
        if abs(theta_inc) >= pi/2
          continue;
        end
        theta_tx = asin(sin(theta_inc)/sqrt(er_ice));
        
        y(doa_idx,rline) = inc_y + sin(theta_tx) ...
          * (twtt(doa_idx,rline)-twtt_to_surf(doa_idx)) * c/2/sqrt(er_ice);
        z(doa_idx,rline) = inc_z - cos(theta_tx) ...
          * (twtt(doa_idx,rline)-twtt_to_surf(doa_idx)) * c/2/sqrt(er_ice);
      else
        % Surface does not exist for this angle of arrival
      end
      
    end
    
  end
else
  %% DOA method
  for rline = 1:Nx
    % Prepare doa (or theta), surface, and twtt for this specific range-line
    tmp_doa = doa(:,:,rline);
    tmp_doa = tmp_doa(~isnan(tmp_doa));
    
    % Setting tmp_doa entries that are outside the pre-specified DOA limits
    % to NaN upsets interp1. So, process them now, but remove in the end.
    doa_trim_idx = find(tmp_doa<doa_limits(1) | tmp_doa>doa_limits(2));
    if ~isempty(tmp_doa)
      [tmp_doa tmp_doa_idx] = sort(tmp_doa,'ascend');
      tmp_theta = tmp_doa;
      
      tmp_surface = surface(:,rline);
      tmp_surface = tmp_surface(~isnan(tmp_surface));
      tmp_surface = tmp_surface(tmp_doa_idx);
      
      tmp_twtt = twtt(:,rline);
      tmp_twtt = tmp_twtt(~isnan(tmp_twtt));
      tmp_twtt = tmp_twtt(tmp_doa_idx);
      
      Nsig = length(tmp_twtt);
      Nsv = Nsig;
      % Convert surface for this along-track slope to y,z coordinates
      surf_y = sin(tmp_theta) .* tmp_surface * c/2;
      surf_z = -cos(tmp_theta) .* tmp_surface * c/2;
      p = polyfit(surf_y,surf_z,1);
      
      % Determine the time to the surface and the position of incidence
      twtt_to_surf = interp1(tmp_theta, tmp_surface, tmp_doa);
      for doa_idx = 1:Nsig
        inc_y = sin(tmp_doa(doa_idx)) * twtt_to_surf(doa_idx) * c/2;
        inc_z = -cos(tmp_doa(doa_idx)) * twtt_to_surf(doa_idx) * c/2;
        
        % For each time position up to the surface, finding y,z is a cylindrical to
        % cartesian coordinate conversion.
        if tmp_twtt(doa_idx) <= twtt_to_surf(doa_idx)
          y(doa_idx,rline) = sin(tmp_doa(doa_idx)) * tmp_twtt(doa_idx) * c/2;
          z(doa_idx,rline) = -cos(tmp_doa(doa_idx)) * tmp_twtt(doa_idx) * c/2;
          
        elseif tmp_twtt(doa_idx) > twtt_to_surf(doa_idx)
          
          % Interpolation indices for surface
          upper_bound = find(tmp_theta > tmp_doa(doa_idx),1);
          if isempty(upper_bound)
            upper_bound = Nsv;
          end
          
          % IDEAL: Determine the local surface slope vector... surface needs to be smooth or the gradient/slope will be noisy
          %  (ignoring along-track slope since SAR focussing should accomodate for)
          % NOT USING IDEAL: surface slope approximated by single polynomial
          % for each range line
          % theta_slope = atan(p(1));
          % NOT USING IDEAL: surface slope approximated by zero
          
          theta_slope = 0;
          theta_inc = tmp_doa(doa_idx) + theta_slope;
          if abs(theta_inc) >= pi/2
            continue;
          end
          theta_tx = asin(sin(theta_inc)/sqrt(er_ice));
          
          y(doa_idx,rline) = inc_y + sin(theta_tx) ...
            * (tmp_twtt(doa_idx)-twtt_to_surf(doa_idx)) * c/2/sqrt(er_ice);
          z(doa_idx,rline) = inc_z - cos(theta_tx) ...
            * (tmp_twtt(doa_idx)-twtt_to_surf(doa_idx)) * c/2/sqrt(er_ice);
        else
          % Surface does not exist for this angle of arrival
        end
        
      end
      % Set points outside the desired DOA limits to NaN.
      y(doa_trim_idx,rline) = NaN;
      z(doa_trim_idx,rline) = NaN;
    end
  end
end
warning('on','MATLAB:interp1:NaNinY');

return;

