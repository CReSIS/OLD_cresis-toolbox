if 0
  block_fn = 'E:\tmp\google_drive\radar\20180819-215243.mat';
  gps_fn = 'E:\tmp\google_drive\gps\180819.pos';
  param.day_seg = '20180819_01';
end
if 1
  block_fn = 'E:\tmp\google_drive\radar\20180523-233532.mat';
  gps_fn = 'E:\tmp\google_drive\gps\180523.pos';
  param.day_seg = '20180523_01';
end

[block_fn_dir,block_fn_name,block_fn_ext] = fileparts(block_fn);
load(block_fn);
param.season_name = '2018_Alaska_SO';
param.radar_name = 'rds';
% params = read_param_xls(ct_filename_param('rds_param_2018_Alaska_SO.xls'));

% data = block.ch0;
% data = bsxfun(@minus,data,mean(data,2));
% data = fir_dec(data,50);
% imagesc(lp(data,20));

records = [];
records.param_records.records.gps.time_offset = -40;
records.gps_time = datenum_to_epoch(block.time);
% Assume block.time is UTC so add leap seconds to get GPS time
records.gps_time = records.gps_time + utc_leap_seconds(records.gps_time(1)) + records.param_records.records.gps.time_offset;
if ~isempty(gps_fn)
  % University of Alaska Fairbanks, Chris Larsen, LIDAR CSV
  % TimeOfDay(UTC) PosLat(deg) PosLon(deg) PosHeight(m) AngleRoll(deg) AnglePitch(deg) Heading(deg)
  % 73652.044 59.51011043 -139.66689705 18.065 -0.853 8.523 206.700
  % OR 
  % 73667.680 59.50971898 -139.66655443 17.868 -0.618 8.647 125.636 0.026  0.020 0.035 
  
  gps_param = [];
  gps_param.format_str = '%f%f%f%f%f%f%f%f%f%f';
  gps_param.types = {'sec','lat_deg','lon_deg','elev_m','roll_deg','pitch_deg','heading_deg','f1','f2','f3'};
  gps_param.textscan = {};
  gps_param.headerlines = 1;
  gps_param.time_reference = 'utc';
  [year month day] = datevec(block.time(1));
  gps_param.year = year;
  gps_param.month = month;
  gps_param.day = day;
  gps = read_gps_general_ascii(gps_fn,gps_param);
  records.gps_source = 'UALidar-final';
else
  % NMEA data
  good_idxs = 1 + find(abs(diff(block.lat))>0 | abs(diff(block.lon))>0);
  gps.gps_time = records.gps_time(good_idxs);
  gps.lat = block.lat(good_idxs);
  gps.lon = block.lon(good_idxs);
  gps.elev = block.elev_air(good_idxs);
  keyboard
  gps.roll = zeros(size(gps.gps_time));
  gps.pitch = zeros(size(gps.gps_time));
  gps.heading = zeros(size(gps.gps_time));
  records.gps_source = 'nmea-field';
end

records.lat = interp1(gps.gps_time,gps.lat,records.gps_time);
records.lon = interp1(gps.gps_time,gps.lon,records.gps_time);
records.elev = interp1(gps.gps_time,gps.elev,records.gps_time);
records.roll = interp1(gps.gps_time,gps.roll,records.gps_time);
records.pitch = interp1(gps.gps_time,gps.pitch,records.gps_time);
records.heading = interp1(gps.gps_time,gps.heading,records.gps_time);

records.file_version = '1L';
records.notes = '';
records.offset = 1:size(records.gps_time,2);
records.radar_name = 'rds';
records.raw = [];
records.relative_filename = {{[block_fn_name,block_fn_ext]}};
records.relative_rec_num = {1};
records.settings.wfs_records = 1;
wf = 1;
records.settings.wfs(wf).Tpd = abs(block.chirp.len);
records.settings.wfs(wf).fc = block.chirp.cf;
records.settings.wfs(wf).fs = 1 / block.dt;
records.settings.wfs(wf).dt = block.dt;
records.settings.wfs(wf).BW = abs(block.chirp.bw);
records.settings.wfs(wf).f0 = block.chirp.cf - block.chirp.cf*block.chirp.bw/200;
records.settings.wfs(wf).f1 = block.chirp.cf + block.chirp.cf*block.chirp.bw/200;
records.settings.wfs(wf).t0 = block.twtt(1);
records.settings.wfs(wf).time = block.twtt;
records.settings.wfs(wf).presums = block.stack;
records.settings.wfs(wf).num_sam = size(block.ch0,1);
records.surface = block.twtt_surf;

records_fn = ct_filename_support(param,'','records');
fprintf('Saving %s\n', records_fn);
records_fn_dir = fileparts(records_fn);
if ~exist(records_fn_dir,'dir')
  mkdir(records_fn_dir);
end
save(records_fn,'-struct','records');
records_aux_files_create(records_fn);

global gRadar;
autogenerate_frames(param,gRadar)
