function gps = merge_sync_gps(gps,merge_type)
% gps = merge_sync_gps(gps,merge_type)
%
% Merges sync NMEA and regular GPS data when sync NMEA has a longer extent than
% the regular GPS data. This function is called from make_gps.m.
% 
% Author: Logan Smith

if merge_type == 1 % sync_gps vector is longer than gps vector
  warning('gps time vector is shorter than the sync gps time vector. Merging gps variables with sync gps variables...')
  Ng = length(gps.gps_time);
  Ns = length(gps.sync_gps_time);
  gps_start_idx = find(gps.sync_gps_time >= gps.gps_time(1),1);
  gps_stop_idx = find(gps.sync_gps_time <= gps.gps_time(end),1,'last');
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    tmp_ext = Ns-gps_stop_idx;
    tmp_len = Ng+tmp_ext;
    sync_start_idx = gps_stop_idx+1;
    sync_stop_idx = Ns;
  elseif gps_stop_idx == Ns % gps_time ends after sync_gps_time
    tmp_ext = gps_start_idx-1;
    tmp_len = Ng+tmp_ext;
    sync_start_idx = 1;
    sync_stop_idx = gps_start_idx - 1;
  else
    tmp_len = Ns; % gps_time is encompassed by sync_gps_time
    sync_start_idx = gps_stop_idx + 1;
    sync_stop_idx = gps_start_idx - 1;
  end
  
  tmp_gps_time = zeros(1,tmp_len);
  tmp_lat = zeros(1,tmp_len);
  tmp_lon = zeros(1,tmp_len);
  tmp_elev = zeros(1,tmp_len);
  tmp_pitch = zeros(1,tmp_len);
  tmp_roll = zeros(1,tmp_len);
  tmp_heading = zeros(1,tmp_len);
  
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    tmp_gps_time(1:Ng) = gps.gps_time;
    tmp_gps_time((Ng+1):end) = gps.sync_gps_time(sync_start_idx:end);
  elseif tmp_len > Ns % gps_time ends after sync_gps_time
    tmp_gps_time(gps_start_idx:gps_start_idx+Ng-1) = gps.gps_time;
    tmp_gps_time(1:sync_stop_idx) = gps.sync_gps_time(1:sync_stop_idx);
  else % gps_time is encompassed by sync_gps_time
    tmp_gps_time(gps_start_idx:gps_start_idx+Ng-1) = gps.gps_time;
    tmp_gps_time(1:sync_stop_idx) = gps.sync_gps_time(1:sync_stop_idx);
    tmp_gps_time(sync_start_idx:end) = gps.sync_gps_time(sync_start_idx:end);
  end
  
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    tmp_lat(1:Ng) = gps.lat;
    tmp_lat((Ng+1):end) = gps.sync_lat(sync_start_idx:end);
  elseif tmp_len > Ns % gps_time ends after sync_gps_time
    tmp_lat(gps_start_idx:gps_start_idx+Ng-1) = gps.lat;
    tmp_lat(1:sync_stop_idx) = gps.sync_lat(1:sync_stop_idx);
  else % gps_time is encompassed by sync_gps_time
    tmp_lat(gps_start_idx:gps_start_idx+Ng-1) = gps.lat;
    tmp_lat(1:sync_stop_idx) = gps.sync_lat(1:sync_stop_idx);
    tmp_lat(sync_start_idx:end) = gps.sync_lat(sync_start_idx:end);
  end
  
  if abs((gps.sync_lon(sync_stop_idx) - gps.lon(1))) > 1e-2
    gps.sync_lon = mod(gps.sync_lon,360);
  end
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    tmp_lon(1:Ng) = gps.lon;
    tmp_lon((Ng+1):end) = gps.sync_lon(sync_start_idx:end);
  elseif tmp_len > Ns % gps_time ends after sync_gps_time
    tmp_lon(gps_start_idx:gps_start_idx+Ng-1) = gps.lon;
    tmp_lon(1:sync_stop_idx) = gps.sync_lon(1:sync_stop_idx);
  else % gps_time is encompassed by sync_gps_time
    tmp_lon(gps_start_idx:gps_start_idx+Ng-1) = gps.lon;
    tmp_lon(1:sync_stop_idx) = gps.sync_lon(1:sync_stop_idx);
    tmp_lon(sync_start_idx:end) = gps.sync_lon(sync_start_idx:end);
  end
  
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    elev_shift = gps.elev(end) - gps.sync_elev(find(gps.sync_gps_time >= gps.gps_time(end),1))
    tmp_elev(1:Ng) = gps.elev;
    tmp_elev((Ng+1):end) = gps.sync_elev(sync_start_idx:end) + elev_shift;
  elseif tmp_len > Ns % gps_time ends after sync_gps_time
    elev_shift = gps.elev(1) - gps.sync_elev(find(gps.sync_gps_time >= gps.gps_time(1),1))
    tmp_elev(gps_start_idx:gps_start_idx+Ng-1) = gps.elev;
    tmp_elev(1:sync_stop_idx) = gps.sync_elev(1:sync_stop_idx) + elev_shift;
  else % gps_time is encompassed by sync_gps_time
    elev_shift = gps.elev(1) - gps.sync_elev(find(gps.sync_gps_time >= gps.gps_time(1),1))
    tmp_elev(gps_start_idx:gps_start_idx+Ng-1) = gps.elev;
    tmp_elev(1:sync_stop_idx) = gps.sync_elev(1:sync_stop_idx) + elev_shift;
    tmp_elev(sync_start_idx:end) = gps.sync_elev(sync_start_idx:end) + elev_shift;
  end
  
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    tmp_roll(1:Ng) = gps.roll;
    tmp_roll((Ng+1):end) = zeros(1,Ns-sync_start_idx+1);
  elseif tmp_len > Ns % gps_time ends after sync_gps_time
    tmp_roll(gps_start_idx:gps_start_idx+Ng-1) = gps.roll;
    tmp_roll(1:sync_stop_idx) = zeros(1,sync_stop_idx);
  else % gps_time is encompassed by sync_gps_time
    tmp_roll(gps_start_idx:gps_start_idx+Ng-1) = gps.roll;
    tmp_roll(1:sync_stop_idx) = zeros(1,sync_stop_idx);
    tmp_roll(sync_start_idx:end) = zeros(1,Ns-sync_start_idx+1);
  end
  
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    tmp_pitch(1:Ng) = gps.pitch;
    tmp_pitch((Ng+1):end) = zeros(1,Ns-sync_start_idx+1);
  elseif tmp_len > Ns % gps_time ends after sync_gps_time
    tmp_pitch(gps_start_idx:gps_start_idx+Ng-1) = gps.pitch;
    tmp_pitch(1:sync_stop_idx) = zeros(1,sync_stop_idx);
  else % gps_time is encompassed by sync_gps_time
    tmp_pitch(gps_start_idx:gps_start_idx+Ng-1) = gps.pitch;
    tmp_pitch(1:sync_stop_idx) = zeros(1,sync_stop_idx);
    tmp_pitch(sync_start_idx:end) = zeros(1,Ns-sync_start_idx+1);
  end
  
  if gps_start_idx == 1 % gps_time starts before sync_gps_time
    tmp_heading(1:Ng) = gps.heading;
    tmp_heading((Ng+1):end) = gps.sync_heading(sync_start_idx:end);
  elseif tmp_len > Ns % gps_time ends after sync_gps_time
    tmp_heading(gps_start_idx:gps_start_idx+Ng-1) = gps.heading;
    tmp_heading(1:sync_stop_idx) = gps.sync_heading(1:sync_stop_idx);
  else % gps_time is encompassed by sync_gps_time
    tmp_heading(gps_start_idx:gps_start_idx+Ng-1) = gps.heading;
    tmp_heading(1:sync_stop_idx) = gps.sync_heading(1:sync_stop_idx);
    tmp_heading(sync_start_idx:end) = gps.sync_heading(sync_start_idx:end);
  end
  
elseif merge_type == 2 % sync_gps time range is longer than gps time range
  Ng = length(gps.gps_time);
  Ns = length(gps.sync_gps_time);
  gps_start_idx = find(gps.sync_gps_time >= gps.gps_time(1),1);
  gps_stop_idx = find(gps.sync_gps_time <= gps.gps_time(end),1,'last');
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    tmp_ext = Ns-gps_stop_idx;
    tmp_len = Ng+tmp_ext;
    sync_start_idx = gps_stop_idx+1;
    sync_stop_idx = Ns;
  elseif gps_stop_idx == Ns % sync_gps_time starts before gps_time
    tmp_ext = gps_start_idx-1;
    tmp_len = Ng+tmp_ext;
    sync_start_idx = 1;
    sync_stop_idx = gps_start_idx - 1;
  end
  
  tmp_gps_time = zeros(1,tmp_len);
  tmp_lat = zeros(1,tmp_len);
  tmp_lon = zeros(1,tmp_len);
  tmp_elev = zeros(1,tmp_len);
  tmp_pitch = zeros(1,tmp_len);
  tmp_roll = zeros(1,tmp_len);
  tmp_heading = zeros(1,tmp_len);
  
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    tmp_gps_time(1:Ng) = gps.gps_time;
    tmp_gps_time((Ng+1):end) = gps.sync_gps_time(sync_start_idx:end);
  elseif gps_stop_idx == Ns % sync_gps_time starts before gps_time
    tmp_gps_time(gps_start_idx:gps_start_idx+Ng-1) = gps.gps_time;
    tmp_gps_time(1:sync_stop_idx) = gps.sync_gps_time(1:sync_stop_idx);
  end
  
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    tmp_lat(1:Ng) = gps.lat;
    tmp_lat((Ng+1):end) = gps.sync_lat(sync_start_idx:end);
  elseif tmp_len > Ns % sync_gps_time starts before gps_time
    tmp_lat(gps_start_idx:gps_start_idx+Ng-1) = gps.lat;
    tmp_lat(1:sync_stop_idx) = gps.sync_lat(1:sync_stop_idx);
  end
  
  if abs((gps.sync_lon(sync_stop_idx) - gps.lon(1))) > 1e-2
    gps.sync_lon = mod(gps.sync_lon,360);
  end
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    tmp_lon(1:Ng) = gps.lon;
    tmp_lon((Ng+1):end) = gps.sync_lon(sync_start_idx:end);
  elseif tmp_len > Ns % sync_gps_time starts before gps_time
    tmp_lon(gps_start_idx:gps_start_idx+Ng-1) = gps.lon;
    tmp_lon(1:sync_stop_idx) = gps.sync_lon(1:sync_stop_idx);
  end
  
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    elev_shift = gps.elev(end) - gps.sync_elev(find(gps.sync_gps_time >= gps.gps_time(end),1))
    tmp_elev(1:Ng) = gps.elev;
    tmp_elev((Ng+1):end) = gps.sync_elev(sync_start_idx:end) + elev_shift;
  elseif tmp_len > Ns % sync_gps_time starts before gps_time
    elev_shift = gps.elev(1) - gps.sync_elev(find(gps.sync_gps_time >= gps.gps_time(1),1))
    tmp_elev(gps_start_idx:gps_start_idx+Ng-1) = gps.elev;
    tmp_elev(1:sync_stop_idx) = gps.sync_elev(1:sync_stop_idx) + elev_shift;
  end
  
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    tmp_roll(1:Ng) = gps.roll;
    tmp_roll((Ng+1):end) = zeros(1,Ns-sync_start_idx+1);
  elseif tmp_len > Ns % sync_gps_time starts before gps_time
    tmp_roll(gps_start_idx:gps_start_idx+Ng-1) = gps.roll;
    tmp_roll(1:sync_stop_idx) = zeros(1,sync_stop_idx);
  end
  
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    tmp_pitch(1:Ng) = gps.pitch;
    tmp_pitch((Ng+1):end) = zeros(1,Ns-sync_start_idx+1);
  elseif tmp_len > Ns % sync_gps_time starts before gps_time
    tmp_pitch(gps_start_idx:gps_start_idx+Ng-1) = gps.pitch;
    tmp_pitch(1:sync_stop_idx) = zeros(1,sync_stop_idx);
  end
  
  if gps_start_idx == 1 % sync_gps_time extends beyond gps_time
    tmp_heading(1:Ng) = gps.heading;
    tmp_heading((Ng+1):end) = gps.sync_heading(sync_start_idx:end);
  elseif tmp_len > Ns % sync_gps_time starts before gps_time
    tmp_heading(gps_start_idx:gps_start_idx+Ng-1) = gps.heading;
    tmp_heading(1:sync_stop_idx) = gps.sync_heading(1:sync_stop_idx);
  end
end

  

fprintf('Check gps merging results now...\n')
keyboard
if 0
  plot(tmp_gps_time,'.')
  plot(diff(tmp_gps_time),'.')
  plot(tmp_elev,'.')
  plot(diff(tmp_elev),'.')
  plot(tmp_lon,tmp_lat,'.')
  plot(tmp_roll,'.')
  plot(tmp_pitch,'.')
  plot(tmp_heading,'.')
end

gps.gps_time = tmp_gps_time;
gps.lat = tmp_lat;
gps.lon = tmp_lon;
gps.elev = tmp_elev;
gps.roll = tmp_roll;
gps.pitch = tmp_pitch;
gps.heading = tmp_heading;
return