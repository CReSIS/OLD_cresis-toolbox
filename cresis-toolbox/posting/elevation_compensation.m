function [mdata,depth_axis] = elevation_compensation(mdata,param)
% function [mdata,depth_axis] = elevation_compensation(mdata,param)
%
% mdata = L1B data structure from load_L1B.m. Requires
%   .Data is Nt x Nx matrix
%   .Time is Nt x 1 vector
%   .Elevation is 1 x Nx vector
%   .Surface is 1 x Nx vector
% param = structure controlling compensation
%  .update_surf = logical, default false, if true the function uses
%    threshold and sidelobe fields to track the surface again
%    This will cause "mdata.Surface" to be updated.
%  .filter_surf = logical, default false, applies an along-track median
%    filter to the surface.
%    This will cause "mdata.Surface" to be updated.
%  .threshold = used with update_surf (specified in power magnitude
%    assuming mdata.Data is in power magnitude units), this threshold
%    will be used for tracking the surface.  Default is estimated from
%    the data.
%  .sidelobe = special case which helps the tracker avoid tracking side
%    lobes at the relative power specified or below. Actual threshold
%    used is:
%     max of param.threshold and max_value_of_current_rangeline*param.sidelobe
%    Default is "10^(-13/10)
%  .elev_comp = elevation compensation mode (same as publish_echogram.m)
%  .er_ice = default is 3.15, vector of dielectrics
%  .er_depth = default is 0, depth vector to align with er_ice, first element
%     must always be zero
%  .depth = string to be evaluated in the calculation of depth_axis.
%     This specifies which WGS-84 elevation bins to include in case a
%     subset of the data is desired. A vector "Surface_Elev" is available
%     which contains the ice surface elevation at each range line.
%     DEFAULT: '[-inf inf]';
%     EXAMPLE: '[min(Surface_Elev)-600 max(Surface_Elev)+20]';
%     EXAMPLE: '[1200 1800]';

% mdata = updated data structure
%  .Data = updated to be gridded on constant WGS-84 axes
%  .Elevation, .Surface_Elev = updated to account for compensation
%    of .Data on constant grid (only for mode 3)
% depth_axis = Depends on elevation compensation mode
%   mode 2: depth axis
%   mode 3: indices specified by param.depth
%
% Examples: See plot_L1B.m
%
% Author: John Paden
%
% See also: load_L1B.m, plot_L1B.m

if ~exist('param','var')
  param = [];
end

if ~isfield(param,'update_surf')
  param.update_surf = false;
end

if ~isfield(param,'filter_surf')
  param.filter_surf = false;
end

if ~isfield(param,'threshold')
  param.threshold = median(max(mdata.Data)) / 10;
end

if ~isfield(param,'sidelobe')
  param.sidelobe = 10^(-13/10);
end

if ~isfield(param,'er_depth')
  param.er_depth = 0;
end
param.er_depth = reshape(param.er_depth,[length(param.er_depth) 1]);

if ~isfield(param,'er_ice')
  param.er_ice = 3.15;
end
param.er_ice = reshape(param.er_ice,[length(param.er_ice) 1]);

if ~isfield(param,'er_freq')
  param.er_freq = 150e6;
end

if ~isfield(param,'depth')
  param.depth = '[-inf inf]';
end

physical_constants;

if param.update_surf
  %% Update surface:
  for rline = 1:size(mdata.Data,2)
    new_surf = find(mdata.Data(:,rline) ...
      > max(param.threshold,max(mdata.Data(:,rline))*param.sidelobe),1);
    if isempty(new_surf)
      mdata.Surface_Bin(rline) = NaN;
    else
      mdata.Surface_Bin(rline) = new_surf;
    end
  end
  mdata.Surface = interp1(1:length(mdata.Time),mdata.Time,mdata.Surface_Bin);
  nan_mask = isnan(mdata.Surface);
  mdata.Surface = interp1(find(~nan_mask),mdata.Surface(~nan_mask),1:length(mdata.Surface));
end

if param.filter_surf
  %% Filter surface
  if length(mdata.Surface) >= 100
    [Bfilt,Afilt] = butter(2,0.02);
    surf_filt = filtfilt(Bfilt,Afilt,mdata.Surface);
    tmp = polyval(polyfit(1:51,reshape(mdata.Surface(1:51),[1 51]),2),1:51);
    surf_filt(1:50) = tmp(1:50) - tmp(51) + surf_filt(51);
    tmp = polyval(polyfit(1:51,reshape(mdata.Surface(end-50:end),[1 51]),2),1:51);
    surf_filt(end-49:end) = tmp(2:51) - tmp(1) + surf_filt(end-50);
    mdata.Surface = surf_filt;
  else
    surf_filt = medfilt1(mdata.Surface,round(length(mdata.Surface)/20)*2+1);
    mdata.Surface = surf_filt;
  end
else
  surf_filt = mdata.Surface;
end

%% Remove data before zero time
negative_bins = mdata.Time < 0;
mdata.Time = mdata.Time(~negative_bins);
mdata.Data = mdata.Data(~negative_bins,:);

if param.elev_comp == 1
  %% Relative elevation compensation
  error('Not supported');
elseif param.elev_comp == 2
  %% Compensate for surface variations (i.e. flatten to low pass filtered
  % version of surface)
  
  % Determine start/stop depths
  
  %% Limit plotted depths according to input param.depth
  % Example: param.depth = '[-2 20]';
  depth_range = eval(param.depth);
  
  % Determine sampling
  max_er = max(param.er_ice);
  dt = mdata.Time(2) - mdata.Time(1);
  dz = dt * c/2/sqrt(max_er);

  % Depth axis
  depth = (round(depth_range(1)/dz)*dz : dz : depth_range(end)).';
  
  % Depth time axis
  depth_time = size(depth);
  % Above surface
  depth_time(depth <= 0) = depth(depth <= 0) / (c/2);
  % Below surface and within dielectric profile defined depth
  if length(param.er_depth) > 1
    TWtime = genPropProfileFromPerm(param.er_depth,param.er_ice,param.er_freq);
    profile_idxs = depth > 0 & depth < param.er_depth(end);
    depth_time(profile_idxs) = interp1(param.er_depth, [0; TWtime], depth(profile_idxs));
  end
  % Below surface and below dielectric profile defined depth
  const_idxs = depth >= param.er_depth(end);
  depth_time(const_idxs) = TWtime(end) + (depth(const_idxs) - param.er_depth(end)) / (c/2/sqrt(param.er_ice(end)));
  
  % Re-interpolate data to constant depth axis
  newData = zeros(length(depth),size(mdata.Data,2));
  for rline = 1:size(mdata.Data,2)
    newData(:,rline) = interp1(mdata.Time, mdata.Data(:,rline), ...
      mdata.Surface(rline) + depth_time);
  end
  mdata.Data = newData;
  
  depth_axis = depth;
  
elseif param.elev_comp == 3
  %% Elevation compensate to WGS-84 y-axis
  
  if length(param.er_ice) > 1
    error('Vector form of er_ice not supported');
  end
  
  %% Create elevation axis to interpolate to
  max_elev = max(mdata.Elevation);
  min_elev = min(mdata.Elevation - surf_filt*c/2 - (mdata.Time(end)-surf_filt)*c/2/sqrt(param.er_ice));
  dt = mdata.Time(2)-mdata.Time(1);
  dr = dt * c/2 / sqrt(param.er_ice);
  elev_axis = max_elev:-dr:min_elev;
  mdata.Elevation_Fasttime = elev_axis;
  
  % Zero pad data to create space for interpolated data
  zero_pad_len = length(elev_axis) - length(mdata.Time);
  mdata.Data = cat(1,mdata.Data,zeros(zero_pad_len,size(mdata.Data,2)));
  
  % Determine the corrections to apply to elevation and layers
  dRange = max_elev - mdata.Elevation;
  dBins = round(dRange / (c/2) / dt);
  dtime = dRange/(c/2);
  
  warning off
  for rline = 1:size(mdata.Data,2)
    % Determine elevation bins before surface
    surf_elev = mdata.Elevation(rline) - surf_filt(rline) * c/2;
    dt_air = dr/(c/2);
    time0 = -(max_elev - mdata.Elevation(rline))/(c/2);
    last_air_idx = find(elev_axis > surf_elev,1,'last');
    new_time = (time0 + dt_air*(0:last_air_idx-1)).';
    if isempty(last_air_idx)
      % Radar is on the surface of the ice
      last_air_idx = 0;
    end
    if last_air_idx < length(elev_axis)
      % Determine elevation bins after surface
      dt_ice = dr/(c/2/sqrt(param.er_ice));
      first_ice_idx = last_air_idx + 1;
      time0 = surf_filt(rline) + (surf_elev - elev_axis(first_ice_idx))/(c/2/sqrt(param.er_ice));
      new_time = cat(1,new_time, (time0 + dt_ice*(0:length(elev_axis)-length(new_time)-1)).');
    end
    mdata.Data(:,rline) = interp1(mdata.Time, mdata.Data(1:length(mdata.Time),rline), new_time, 'linear',0);
    mdata.Elevation(rline) = mdata.Elevation(rline) + dRange(rline);
    mdata.Surface(rline) = mdata.Surface(rline) + dtime(rline);
    %mdata.Bottom(rline) = mdata.Bottom(rline) + dtime(rline);
  end
  warning on
  
  %% Limit plotted depths according to input param.depth
  DSurface = mdata.Elevation - mdata.Surface*c/2;
  Surface_Elev = DSurface;
  mdata.Surface_Elev = DSurface;
  % Example: param.depth = '[min(Surface_Elev) - 15 max(Surface_Elev)+3]';
  % Example: param.depth = '[100 120]';
  depth_range = eval(param.depth);
  depth_axis = find(elev_axis >= depth_range(1) & elev_axis <= depth_range(end));
else
  error('Not supported');
end

return;
