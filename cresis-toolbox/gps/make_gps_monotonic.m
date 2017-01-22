function [gps,error_flag] = make_gps_monotonic(gps)
% [gps,error_flag] = make_gps_monotonic(gps)
%
% Check gps structure for non-monotonically increasing GPS time and return
% a corrected gps structure.
%
% Author: John Paden

error_flag = false;

[~,sort_idxs] = sort(gps.gps_time);
if any(sort_idxs ~= 1:length(sort_idxs))
  warning('GPS time is not monotonically increasing. Manual inspection is suggested. May need to run make_gps_monotonic.m in make_gps_SEASON.m for this particular file.');
  gps.gps_time = gps.gps_time(sort_idxs);
  gps.lat = gps.lat(sort_idxs);
  gps.lon = gps.lon(sort_idxs);
  gps.elev = gps.elev(sort_idxs);
  gps.roll = gps.roll(sort_idxs);
  gps.pitch = gps.pitch(sort_idxs);
  gps.heading = gps.heading(sort_idxs);
  error_flag = true;
end

good_mask = [true, (diff(gps.gps_time) > 0)];

if ~all(good_mask)
  warning('GPS time has repeated records. Manual inspection is suggested. May need to run make_gps_monotonic.m in make_gps_SEASON.m for this particular file.');
  
  gps.gps_time = gps.gps_time(good_mask);
  gps.lat = gps.lat(good_mask);
  gps.lon = gps.lon(good_mask);
  gps.elev = gps.elev(good_mask);
  gps.roll = gps.roll(good_mask);
  gps.pitch = gps.pitch(good_mask);
  gps.heading = gps.heading(good_mask);
  error_flag = true;
end

end
