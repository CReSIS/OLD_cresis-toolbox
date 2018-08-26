function surf_model = surface_gen(param)
% surf_model = surface_gen(param)
%
% param
%  .z: structure describing the statistics of the layer height
%    .mean: mean surface height (either a scalar or an Ny length vector)
%    .rms_height: RMS height of surface
%    .corr_length_x: Correlation length of surface in x-dimension
%    .corr_length_y: Correlation length of surface in y-dimension
%    .pdf: string containing pdf to use (e.g. 'gaussian', 'exponential')
%  .dx: sample spacing in x-dimension
%  .x_range: 2-element vector specifying the start/stop of the x-axis of
%    the surface model
%  .dy: sample spacing in y-dimension
%  .y_range: 2-element vector specifying the start/stop of the y-axis of
%    the surface model
%  .rcs_in: structure describing RCS statistics
%     .mean
%     .var
%     .pdf
%
% surf_model:
%  .dem: digital elevation model for surface which is Ny by Nx matrix
%  .x_axis: x-axis of the surface model is 1 by Nx
%  .y_axis: y-axis of the surface model is Ny by 1
%  .rcs: radar cross section of each pixel in DEM
%
% Author: Sean Holloway, John Paden

surf_model = [];

surf_model.x_axis = param.x_range(1) : param.dx : param.x_range(end);
surf_model.y_axis = (param.y_range(1) : param.dy : param.y_range(end)).';

xsize = length(surf_model.x_axis);
ysize = length(surf_model.y_axis);

%% Create random surface height and radar cross section (RCS)
surf_model.dem = param.z.rms_height .* randn(ysize, xsize);

%% Filter surface height in x-dimension
if param.dx/param.z.corr_length_x < 1
  [bx, ax] = butter(2, param.dx/param.z.corr_length_x);
  for targ = 1:ysize
    surf_model.dem(targ,:) = sqrt(param.z.corr_length_x/param.dx*1.2)...
      *filtfilt(bx,ax,surf_model.dem(targ,:));
  end
end

%% Filter surface height along y-dimension
if param.dy/param.z.corr_length_y < 1
  [by, ay] = butter(2, param.dy/param.z.corr_length_y);
  for targ = 1:xsize
    surf_model.dem(:,targ) = sqrt(param.z.corr_length_y/param.dy*1.12)...
      *filtfilt(by,ay,surf_model.dem(:,targ));
  end
end

%% Add in mean surface height
surf_model.dem = param.z.mean + surf_model.dem;

%% Add in mean RCS
% surf_model.rcs = param.rcs_in.mean + surf_model.rcs;

surf_model.x = param.x;
surf_model.y = param.y;
% interp2 assignes NaN to all queries that lie outside the domain of the
% sample points in all interpolation methods, except for 'spline' and
% 'makima', where interp2 does extrapolation. This problem happened few
% times, so I changed the method to 'spline', or you can add
% interp2(...,'linear',extrapval), where extrapval is a constant value that
% gets assigned to all queries outside the domain of the sample points.
% surf_model.z = interp2(surf_model.x_axis, surf_model.y_axis, surf_model.dem, surf_model.x, surf_model.y,'spline');
% surf_model.rcs = interp2(surf_model.x_axis, surf_model.y_axis, surf_model.rcs, surf_model.x, surf_model.y,'spline');
surf_model.z = interp2(surf_model.x_axis, surf_model.y_axis, surf_model.dem, surf_model.x, surf_model.y,'linear',0);
surf_model.rcs = param.rcs_in.mean + sqrt(param.rcs_in.var/2) .* (randn(length(surf_model.y), length(surf_model.x)) + 1i*randn(length(surf_model.y), length(surf_model.x)));

% surf_model.rcs = interp2(surf_model.x_axis, surf_model.y_axis, surf_model.rcs, surf_model.x, surf_model.y,'linear',0);
surf_model.Nx = length(param.x);

return;
