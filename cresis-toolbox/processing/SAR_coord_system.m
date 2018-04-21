function fcs = SAR_coord_system(param,gps,ref,along_track,output_along_track)
% fcs = SAR_coord_system(param,gps,ref,along_track,output_along_track)
%
% Develops flight coordinate system (FCS) for SAR data. Also produces inputs
% used for motion_comp.m.
%
% param
%  .squint = vector in flight coordinate system, pointing from phase center
%    to the center of the scene, motion compensation, drange, is corrected
%    along this vector
%    default is [0 0 -1] or straight down
%  .type = scalar integer
%    0: no motion compensation
%    1. motion compensation to reference
%    2. motion compensation to line fitted to reference
%  .Lsar = length of SAR aperture (m), filter length to obtain
%    heading estimate used to develop FCS
% gps = struct with fields lat,lon,elev,roll,pitch,heading
%    1 by Nx vectors
%    Not required if param.mocomp_mode is 0 or 1
%    lat,lon in degrees
%    elev in meters
%    roll, pitch, heading in radians
% ref = struct with fields lat,lon,elev,roll,pitch,heading
%    1 by Nx vectors
%    Not required if param.mocomp_mode is 0 or 1
%    lat,lon in degrees
%    elev in meters
%    roll, pitch, heading in radians
% along_track = 1 by Nx along-track vector
% output_along_track = 1 by No SAR output along-track vector positions
%
% fcs (flight-coordinate-system structure)
%   .x = x-axis of FCS for each output position, 3 by No array
%   .y = y-axis of FCS for each output position, 3 by No array
%   .z = z-axis of FCS for each output position, 3 by No array
%   .origin = origin of FCS for each output position, 3 by No array
%   .pos = phase center of measurement in FCS for each output position,
%     3 by No array
%   .roll = platform roll
%   .pitch = platform pitch
%   .heading = platform heading
%
% Flight coordinate system (FCS)
%   origin = first sample
%   x = along-track
%   z = pointing toward local up and orthogonal to x
%   y = completes right handed coordinate system
%   fcs = flight coordinate system
%
% Author: John Paden

fcs.Lsar = param.Lsar;

physical_constants;

[ref_ecef(1,:),ref_ecef(2,:),ref_ecef(3,:)] = ...
  geodetic2ecef(ref.lat/180*pi,ref.lon/180*pi,ref.elev,WGS84.ellipsoid);

[ecef(1,:),ecef(2,:),ecef(3,:)] ...
  = geodetic2ecef(gps.lat/180*pi,gps.lon/180*pi,gps.elev,WGS84.ellipsoid);

if param.type == 4 || param.type == 5
  % Savitsky Golay filtering of reference positions
  dx_in = median(diff(along_track));
  sgolayfilt_F = 2*round(fcs.Lsar / dx_in) + 1;
  sv_ecef = sgolayfilt(ref_ecef.', 1, sgolayfilt_F, hanning(sgolayfilt_F)).';
end

%% Create flight coordinate system (MFCS) at each output position
fcs.x = zeros(3,length(output_along_track));
fcs.y = zeros(3,length(output_along_track));
fcs.z = zeros(3,length(output_along_track));
fcs.origin = zeros(3,length(output_along_track));
fcs.pos = zeros(3,length(output_along_track));
fcs.roll = zeros(1,length(output_along_track));
fcs.pitch = zeros(1,length(output_along_track));
fcs.heading = zeros(1,length(output_along_track));
fcs.surface = zeros(1,length(output_along_track));
fcs.bottom= zeros(1,length(output_along_track));
good_rline = ones(1,length(output_along_track));

for out_rline = 1:length(output_along_track)
  if ~mod(out_rline-1,10)
    fprintf('%d of %d (%s)\n', out_rline, length(output_along_track), datestr(now));
  end
  % For this output range line determine the input range lines
  rlines_in = find(along_track >= output_along_track(out_rline)-param.Lsar/2 ...
    & along_track <= output_along_track(out_rline)+param.Lsar/2);
  
  if length(rlines_in) < 2
    % Sometimes there are gaps in the data, we will use neighboring points
    % to fill in these gaps
    good_rline(out_rline) = 0;
    continue;
  end
  
  fcs.origin(:,out_rline) = mean(ref_ecef(:,rlines_in),2);
  
  fit_ecef = [];
  
  %% Use the SAR aperture length to create an averaged FCS
  if param.type == 4 || param.type == 5
    %% Compensation to savitsky golay
    [~,rlines_fit] = min(abs(along_track-output_along_track(out_rline)));
    if rlines_fit < length(along_track)
      rlines_fit = rlines_fit + [0 1]; 
    else
      rlines_fit = rlines_fit + [-1 0]; 
    end
    fit_ecef = sv_ecef(:,rlines_fit);
  elseif param.type == 3
    %% Compensation to line fit to reference (param.type == 3)
    fit_ecef(1,:) = polyval(polyfit(1:size(ref_ecef,2),ref_ecef(1,:),1),rlines_in);
    fit_ecef(2,:) = polyval(polyfit(1:size(ref_ecef,2),ref_ecef(2,:),1),rlines_in);
    fit_ecef(3,:) = polyval(polyfit(1:size(ref_ecef,2),ref_ecef(3,:),1),rlines_in);
  elseif param.type == 2
    %% Compensation to piece-wise line fit to reference (param.type == 2)
    [p,~,mu] = polyfit(rlines_in,ref_ecef(1,rlines_in),1);
    fit_ecef(1,:) = polyval(p,rlines_in,[],mu);
    [p,~,mu] = polyfit(rlines_in,ref_ecef(2,rlines_in),1);
    fit_ecef(2,:) = polyval(p,rlines_in,[],mu);
    [p,~,mu] = polyfit(rlines_in,ref_ecef(3,rlines_in),1);
    fit_ecef(3,:) = polyval(p,rlines_in,[],mu);
  elseif param.type >= 0 && param.type <= 1
    %% Compensation to reference (param.type == 1)
    fit_ecef = ref_ecef(:,rlines_in);
  else
    error('Motion compensation type must be 0, 1, 2, or 3');
  end
  
  fcs_x = [fit_ecef(1,end)-fit_ecef(1,1); fit_ecef(2,end)-fit_ecef(2,1); fit_ecef(3,end)-fit_ecef(3,1)];
  fcs_x = fcs_x / sqrt(dot(fcs_x,fcs_x));
  
  [start_geo(1),start_geo(2),start_geo(3)] = ecef2geodetic(fit_ecef(1),fit_ecef(2),fit_ecef(3),WGS84.ellipsoid);
  [fit_ecef_up(1,:),fit_ecef_up(2,:),fit_ecef_up(3,:)] ...
    = geodetic2ecef(start_geo(1),start_geo(2),start_geo(3)+1,WGS84.ellipsoid);
  
  local_up = [fit_ecef_up(1)-fit_ecef(1); fit_ecef_up(2)-fit_ecef(2); fit_ecef_up(3)-fit_ecef(3)];
  local_up = local_up / sqrt(dot(local_up,local_up));
  fcs_z = local_up;
  fcs_z = fcs_z - fcs_x*dot(fcs_z,fcs_x);
  fcs_z = fcs_z / sqrt(dot(fcs_z,fcs_z));
  
  % y-axis completes RH coordinate system
  fcs.y(:,out_rline) = cross(fcs_z, fcs_x);
  fcs.x(:,out_rline) = fcs_x;
  fcs.z(:,out_rline) = fcs_z;
  
  % The following operation computes the linear transformation matrix:
  %   fcs_T = [fcs.x(:,out_rline) fcs.y(:,out_rline) fcs.z(:,out_rline)];
  % Uses the "A\B" (help slash) nomenclature to inverse it:
  %   fcs_Tinv = inv(fcs_T);
  fcs.pos(:,out_rline) = [fcs.x(:,out_rline) fcs.y(:,out_rline) fcs.z(:,out_rline)] \ (mean(ecef(:,rlines_in),2) - fcs.origin(:,out_rline));
  
  fcs.surface(:,out_rline) = mean(gps.surface(rlines_in));
  if isfield(gps,'bottom')
    fcs.bottom(:,out_rline) = mean(gps.bottom(rlines_in));
  else
    fcs.bottom(:,out_rline) = NaN;
  end
  
  fcs.roll(:,out_rline) = mean(gps.roll(rlines_in));
  fcs.pitch(:,out_rline) = mean(gps.pitch(rlines_in));
  fcs.heading(:,out_rline) = mean(gps.heading(rlines_in));
end
% GPS time is handles differently since it is used to synchronize
% to other instruments
[~,unique_idxs] = unique(along_track);
fcs.gps_time = interp1(along_track(unique_idxs),gps.gps_time(unique_idxs),output_along_track);

% Fill in all gaps in the output (gaps in output are caused by gaps in
% the input data). In other words, there was no input data for the particular
% output position (where input data is determined by Lsar). We will
% interpolate from neighboring samples to fill in the gaps in the FCS.
interp_rline = find(good_rline == 0,1);
while ~isempty(interp_rline)
  if interp_rline == 1
    % Gap is at first range line so we have to look to see if there are any
    % good data points after it and copy from that
    rline = find(good_rline(interp_rline:end)==1,1) + interp_rline - 1;
    if isempty(rline)
      error('Data gap extends across entire chunk. Consider breaking segment into two at this gap or increasing the chunk size.');
    end
    interp_rline = interp_rline:rline-1;
    fcs.x(:,interp_rline) = repmat(fcs.x(:,rline),[1 size(interp_rline,2)]);
    fcs.y(:,interp_rline) = repmat(fcs.y(:,rline),[1 size(interp_rline,2)]);
    fcs.z(:,interp_rline) = repmat(fcs.z(:,rline),[1 size(interp_rline,2)]);
    fcs.origin(:,interp_rline) = repmat(fcs.origin(:,rline),[1 size(interp_rline,2)]);
    fcs.pos(:,interp_rline) = repmat(fcs.pos(:,rline),[1 size(interp_rline,2)]);
    fcs.roll(:,interp_rline) = repmat(fcs.roll(:,rline),[1 size(interp_rline,2)]);
    fcs.pitch(:,interp_rline) = repmat(fcs.pitch(:,rline),[1 size(interp_rline,2)]);
    fcs.heading(:,interp_rline) = repmat(fcs.heading(:,rline),[1 size(interp_rline,2)]);
    fcs.surface(:,interp_rline) = repmat(fcs.surface(:,rline),[1 size(interp_rline,2)]);
    fcs.bottom(:,interp_rline) = repmat(fcs.bottom(:,rline),[1 size(interp_rline,2)]);
  else
    % Gap is not at first range line, so there is a preceding good point.
    % We look to see if there is a good point after this gap. If there is
    % not then we just use the preceding good point. If there is, then we
    % interpolate between the two good points.
    rline = find(good_rline(interp_rline:end)==1,1) + interp_rline - 1;
    % The last good data point is always the range line before this one
    rline0 = interp_rline-1;
    if isempty(rline)
      rline = rline0;
      interp_rline = interp_rline:length(output_along_track);
      fcs.x(:,interp_rline) = repmat(fcs.x(:,rline),[1 size(interp_rline,2)]);
      fcs.y(:,interp_rline) = repmat(fcs.y(:,rline),[1 size(interp_rline,2)]);
      fcs.z(:,interp_rline) = repmat(fcs.z(:,rline),[1 size(interp_rline,2)]);
      fcs.origin(:,interp_rline) = repmat(fcs.origin(:,rline),[1 size(interp_rline,2)]);
      fcs.pos(:,interp_rline) = repmat(fcs.pos(:,rline),[1 size(interp_rline,2)]);
      fcs.roll(:,interp_rline) = repmat(fcs.roll(:,rline),[1 size(interp_rline,2)]);
      fcs.pitch(:,interp_rline) = repmat(fcs.pitch(:,rline),[1 size(interp_rline,2)]);
      fcs.heading(:,interp_rline) = repmat(fcs.heading(:,rline),[1 size(interp_rline,2)]);
      fcs.surface(:,interp_rline) = repmat(fcs.surface(:,rline),[1 size(interp_rline,2)]);
      fcs.bottom(:,interp_rline) = repmat(fcs.bottom(:,rline),[1 size(interp_rline,2)]);
      break;
    end
    interp_rline = rline0+1:rline-1;
    fcs.x(:,interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.x(:,[rline0 rline]).',output_along_track(interp_rline)).';
    fcs.y(:,interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.y(:,[rline0 rline]).',output_along_track(interp_rline)).';
    fcs.z(:,interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.z(:,[rline0 rline]).',output_along_track(interp_rline)).';
    fcs.origin(:,interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.origin(:,[rline0 rline]).',output_along_track(interp_rline)).';
    fcs.pos(:,interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.pos(:,[rline0 rline]).',output_along_track(interp_rline)).';
    fcs.roll(interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.roll([rline0 rline]),output_along_track(interp_rline));
    fcs.pitch(interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.pitch([rline0 rline]),output_along_track(interp_rline));
    fcs.heading(interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.heading([rline0 rline]),output_along_track(interp_rline));
    fcs.surface(interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.surface([rline0 rline]),output_along_track(interp_rline));
    fcs.bottom(interp_rline) = interp1(output_along_track([rline0 rline]), ...
      fcs.bottom([rline0 rline]),output_along_track(interp_rline));
  end
  % Find next gap that needs to be interpolated over
  interp_rline = interp_rline(end)+1;
  interp_rline = find(good_rline(interp_rline:end) == 0,1) + interp_rline - 1;
end

return;


