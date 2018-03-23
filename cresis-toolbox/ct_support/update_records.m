function update_records(param,param_override)
% update_records(param,param_override)
%
% This function updates the records according to the param struct. It
% really just resyncs the gps file to the records file. This is useful when
% the gps file is updated or the param.vectors.gps.time_offset is updated.
%
% params: parameter spreadsheet structure array
%
% Examples: See run_update_records.m
%
% Author: John Paden
%
% See also: run_update_records.m

param = merge_structs(param,param_override);

dbstack_info = dbstack;
fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', dbstack_info(1).name, param.day_seg, datestr(now,'HH:MM:SS'));
fprintf('=====================================================================\n');

%% Input checks

% save_changes: Logical, For debugging purposes, you can turn the file save on/off
save_changes = true;

%% Prep (load records and gps files)
records_fn = ct_filename_support(param,'','records');
records = load(records_fn);

gps_fn = ct_filename_support(param,'','gps',1);

fprintf('  Using gps file %s\n', gps_fn);
gps = load(gps_fn);

%% Update the records
tmp_records = records;
if any(strcmpi(param.radar_name,{'accum2','mcrds'}))
  % Determine time offset delta and apply to radar time
  delta_offset = param.vectors.gps.time_offset - records.param_records.vectors.gps.time_offset
  records.param_records.vectors.gps.time_offset = param.vectors.gps.time_offset;
  records.raw.radar_time = records.raw.radar_time + delta_offset;
  
  % Find the section of the GPS file that is to be used for this record
  guard_time = 5;
  good_idxs = find(gps.comp_time >= records.raw.comp_time(1)-guard_time ...
    & gps.comp_time <= records.raw.comp_time(end)+guard_time);
  good_radar_time = gps.radar_time(good_idxs);
  good_sync_gps_time = gps.sync_gps_time(good_idxs);
  % From these good indexes, remove any repeat radar times (usually caused
  % by there being more than one NMEA string every 1 PPS
  good_idxs = 1+find(diff(good_radar_time) ~= 0);
  good_radar_time = good_radar_time(good_idxs);
  good_sync_gps_time = good_sync_gps_time(good_idxs);
  % Interpolate gps.sync_gps_time to records.gps_time using gps.radar_time
  % and records.radar_time
  records.gps_time = interp1(good_radar_time, good_sync_gps_time, ...
    records.raw.radar_time,'linear','extrap');
  
  % Interpolate to find new trajectory
  records.lat = double(interp1(gps.gps_time,gps.lat,records.gps_time));
  records.lon = double(mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),records.gps_time)*180/pi+180, 360)-180);
  records.elev = double(interp1(gps.gps_time,gps.elev,records.gps_time));
  records.roll = double(interp1(gps.gps_time,gps.roll,records.gps_time));
  records.pitch = double(interp1(gps.gps_time,gps.pitch,records.gps_time));
  records.heading = double(mod(interp1(gps.gps_time,unwrap(gps.heading),records.gps_time)+pi,2*pi)-pi);
  records.gps_source = gps.gps_source;
else
  % Determine time offset delta and apply to radar time
  delta_offset = param.vectors.gps.time_offset - records.param_records.vectors.gps.time_offset
  records.param_records.vectors.gps.time_offset = param.vectors.gps.time_offset;
  records.gps_time = records.gps_time + delta_offset;
  
  % Interpolate to find new trajectory
  records.lat = interp1(gps.gps_time, gps.lat, records.gps_time);
  records.lon = mod(interp1(gps.gps_time,unwrap(gps.lon/180*pi),records.gps_time)*180/pi+180, 360)-180;
  records.elev = interp1(gps.gps_time, gps.elev, records.gps_time);
  records.roll = interp1(gps.gps_time, gps.roll, records.gps_time);
  records.pitch = interp1(gps.gps_time, gps.pitch, records.gps_time);
  records.heading = mod(interp1(gps.gps_time,unwrap(gps.heading),records.gps_time)+pi,2*pi)-pi;
  records.gps_source = gps.gps_source;
end

fprintf('  Delta GPS time offset is %.1f sec\n', delta_offset);

if save_changes
  % Save outputs
  fprintf('  Saving records %s\n', records_fn);
  save(records_fn,'-struct','records');
  create_records_aux_files(records_fn);
else
  fprintf('  Not saving information (TEST MODE)\n');
end
