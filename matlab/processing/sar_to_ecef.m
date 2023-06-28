function [x,y,z] = sar_to_ecef(origin,y_axis,z_axis,time,doa,surf,coord_system)
% [x,y,z] = sar_to_ecef(origin,y,z,time,doa,surf,coord_system)
%
% origin: position of measurement in earth centered earth fixed
%   coordinate system(ECEF), (3 by Nx matrix, meters)
% y_axis,z_axis: 3 by Nx matrices giving the y_axis and z_axis of the flight
%   coordinate system which the doa is in reference to. These axes are
%   written in ECEF coordinates. They should be unit length vectors.
%   A doa of zero points toward -z and increases toward positive y.
% surf: round trip time delay to surface (1 by Nx vector, seconds)
% time: time is round trip time delay to target (Nt by 1 vector, seconds)
% doa: corresponding direction of arrival to target (radians)
%   doa matrix must be Nt by Nx by ? (where ? can be any number of
%   additional dimensions)
% coord_system: 'fcs' or 'ecef'
%
% For coord_system == 'ecef':
%  x,y,z: x,y,z are matrices of the same size as doa which give the
%    targets' positions in ECEF.
%  x,y,z: x,y,z are matrices of the same size as doa which give the
%    targets' positions in FCS (x is just along-track).

physical_constants;
doa_size = size(doa);
time = repmat(time,[1 doa_size(2:end)]);
surf = repmat(surf,[doa_size(1) 1 doa_size(3:end)]);

%% Find DOA values which are outside real view (i.e. > or < pi/2)
badmask = abs(doa) > pi/2;
doa(badmask) = NaN;

% er_ice = 1.78;
% Compute the DOA in ice
doa_ice = asin(sin(doa)/sqrt(er_ice));

%H = height above surface (m)
H = surf * c/2;

% R0 = slant range to surface in the DOA (m)
R0 = H./cos(doa);

% R1 = slant range from surface to target in the DOA in ice (m)
R1 = (time - surf./cos(doa)) * c/2/sqrt(er_ice);

% Determine which targets are above the surface
above_mask = time < surf ./ cos(doa);

% Initialize arrays
y_fcs = zeros(size(doa));
z_fcs = zeros(size(doa));

%% Solve for FCS for targets above the surface
y_fcs(above_mask) = time(above_mask)*c/2 .* sin(doa(above_mask));
z_fcs(above_mask) = -time(above_mask)*c/2 .* cos(doa(above_mask));

%% Solve for FCS for targets below the surface

% Add in the offsets to the surface
y_fcs(~above_mask) = R0(~above_mask) .* sin(doa(~above_mask));
z_fcs(~above_mask) = -R0(~above_mask) .* cos(doa(~above_mask));

% Add in the offsets to the target in the ice
y_fcs(~above_mask) = y_fcs(~above_mask) + R1(~above_mask) .* sin(doa_ice(~above_mask));
z_fcs(~above_mask) = z_fcs(~above_mask) - R1(~above_mask) .* cos(doa_ice(~above_mask));

%% Create outputs

if strcmp(coord_system,'fcs')
  x_fcs = [0 sqrt(diff(origin(1,:)).^2 + diff(origin(2,:)).^2 + diff(origin(3,:)).^2)];
  x_fcs = repmat(x_fcs, [doa_size(1) 1 doa_size(3:end)]);
  
  x = x_fcs;
  y = y_fcs;
  z = z_fcs;
else
  
  %% Convert FCS to ECEF
  x = repmat(origin(1,:), [doa_size(1) 1 doa_size(3:end)]) ...
    + repmat(y_axis(1,:), [doa_size(1) 1 doa_size(3:end)]) .* y_fcs ...
    + repmat(z_axis(1,:), [doa_size(1) 1 doa_size(3:end)]) .* z_fcs;
  
  y = repmat(origin(2,:), [doa_size(1) 1 doa_size(3:end)]) ...
    + repmat(y_axis(2,:), [doa_size(1) 1 doa_size(3:end)]) .* y_fcs ...
    + repmat(z_axis(2,:), [doa_size(1) 1 doa_size(3:end)]) .* z_fcs;
  
  z = repmat(origin(3,:), [doa_size(1) 1 doa_size(3:end)]) ...
    + repmat(y_axis(3,:), [doa_size(1) 1 doa_size(3:end)]) .* y_fcs ...
    + repmat(z_axis(3,:), [doa_size(1) 1 doa_size(3:end)]) .* z_fcs;
  
end

return
