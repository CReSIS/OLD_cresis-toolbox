function decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,spacing)
% decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,spacing)
%
% Given a gps structure OR an along track vector (usually meters traveled at each point
% returned from geodetic_to_along_track.m) and a desired spacing
% (in same units as along_track), this function returns the indices
% in the gps structure or along_track which create a minimum spacing between points
% given by the spacing variable.
%
% Mode 1:
% along_track = structure of gps coordinates
%   .lat = N length vector of latitude north in degrees
%   .lon = N length vector of longitude east in degrees
%   .elev = N length vector of elevation relative to WGS84 in meters
% spacing = scalar representing the desired spacing between positions
%
% Mode 2:
% along_track = N length vector of along-track positions
% spacing = scalar representing the desired spacing between positions
%
% decim_idxs = indices of along_track.[lat/lon/elev] or of along_track
% that create the desired spacing
%
% Author: John Paden
%
% See also: geodetic_to_along_track.m

if isstruct(along_track)
  % Mode 1
  gps = along_track;
  
  if isempty(gps.lat)
    decim_idxs = [];
    return;
  end

  physical_constants;
  [ecef(1,:) ecef(2,:) ecef(3,:)] = geodetic2ecef(gps.lat/180*pi, gps.lon/180*pi, gps.elev, WGS84.ellipsoid);
  
  good_mask = zeros(size(gps.lat));
  good_mask(1) = 1;
  cur_pos = ecef(:,1);
  for cur_idx = 2:length(gps.lat)
    if sqrt(sum(abs(ecef(:,cur_idx) - cur_pos).^2)) > spacing
      cur_pos = ecef(:,cur_idx);
      good_mask(cur_idx) = 1;
    end
  end
  decim_idxs = find(good_mask);
else
  % Mode 2
  if isempty(along_track)
    decim_idxs = [];
    return;
  end
  
  good_mask = zeros(size(along_track));
  good_mask(1) = 1;
  cur_pos = along_track(1);
  for cur_idx = 2:length(along_track)
    if along_track(cur_idx) >= cur_pos + spacing
      cur_pos = along_track(cur_idx);
      good_mask(cur_idx) = 1;
    end
  end
  decim_idxs = find(good_mask);
end

return;
