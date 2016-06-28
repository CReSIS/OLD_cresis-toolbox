% function gps = read_gps_arena()

fn = '/users/paden/tmp/19691231_180005_arena-ctu-ctu_NMEA.bin';
fn = '/users/paden/tmp/20160623_140746_arena-ctu-ctu_NMEA.bin';


fid = fopen(fn,'r');

A = fread(fid,14*2,'uint32');

fclose(fid);

finfo = frame_sync_info(fn,struct('sync','1DFCCF1A','file_mode','ieee-le'));

hdr_param = [];
hdr_param.frame_sync = uint32(hex2dec('1ACFFC1D'));
hdr_param.field_offsets = uint32([8:8:48]); % epri sec1 sec2 fractions
hdr_param.field_types = {uint64(1) uint64(1) uint64(1) uint64(1) uint64(1) uint64(1)};
hdr_param.file_mode = 'ieee-le';

[file_size offset rel_time_cntr rel_time_cntr_pps utctime lat lon elev] ...
  = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);


%%
rel_time_cntr = double(rel_time_cntr) / 10e6;

figure(1); clf;
set(1,'WindowStyle','docked');
subplot(2,1,1);
plot(rel_time_cntr)
title('Relative time');
subplot(2,1,2);
plot(diff(double(rel_time_cntr)),'.')
title('Relative time diff');
ylim([-3 3]);

%%
rel_time_cntr_pps = double(rel_time_cntr_pps) / 10e6;
figure(2); clf;
set(2,'WindowStyle','docked');
subplot(2,1,1);
plot(rel_time_cntr_pps)
title('PPS time');
subplot(2,1,2);
plot(diff(double(rel_time_cntr_pps)),'.')
title('PPS time diff');
ylim([-3 3]);

%%
year = 2000 + double(bitand(utctime,15*2^4))*10/2^4 + double(bitand(utctime,15));
month = double(bitand(utctime,15*2^12))*10/2^12 + double(bitand(utctime,15*2^8))/2^8;
day = double(bitand(utctime,15*2^20))*10/2^20 + double(bitand(utctime,15*2^16))/2^16;
hour = double(bitand(utctime,15*2^60))*10/2^60 + double(bitand(utctime,15*2^56))/2^56;
min = double(bitand(utctime,15*2^52))*10/2^52 + double(bitand(utctime,15*2^48))/2^48;
utctime = bitshift(utctime,-4*6);
sec = double(bitand(utctime,15));
utctime = bitshift(utctime,-4);
for idx = 1:5
  sec = sec/10 + double(bitand(utctime,15));
  utctime = bitshift(utctime,-4);
end
sec = sec * 10;

utc_time = datenum(year,month,day,hour,min,sec);
fprintf('%s to %s\n', datestr(utc_time(1)), datestr(utc_time(end)));
utc_time = datenum_to_epoch(utc_time);
utc_time_sod = epoch_to_sod(utc_time);

figure(3); clf;
set(3,'WindowStyle','docked');
subplot(2,1,1);
plot(utc_time_sod)
title('UTC time SOD');
subplot(2,1,2);
plot(diff(double(utc_time_sod)),'.');
title('UTC time SOD diff');
ylim([-3 3]);

%%
latsec = double(bitand(lat,15));
lat = bitshift(lat,-4);
for idx = 1:12
  latsec = latsec/10 + double(bitand(lat,15));
  lat = bitshift(lat,-4);
end
latsec = latsec * 10;

latdeg = double(bitand(lat,15));
lat = bitshift(lat,-4);
for idx = 1:1
  latdeg = latdeg/10 + double(bitand(lat,15));
  lat = bitshift(lat,-4);
end
latdeg = latdeg * 10;

if lat
  lat = -(latdeg + latsec/60);
else
  lat = latdeg + latsec/60;
end

figure(4); clf;
set(4,'WindowStyle','docked');
plot(lat);
title('Latitude');

%%
lonsec = double(bitand(lon,15));
lon = bitshift(lon,-4);
for idx = 1:11
  lonsec = lonsec/10 + double(bitand(lon,15));
  lon = bitshift(lon,-4);
end
lonsec = lonsec * 10;

londeg = double(bitand(lon,15));
lon = bitshift(lon,-4);
for idx = 1:2
  londeg = londeg/10 + double(bitand(lon,15));
  lon = bitshift(lon,-4);
end
londeg = londeg * 100;

if lon
  lon = -(londeg + lonsec/60);
else
  lon = londeg + lonsec/60;
end

figure(5); clf;
set(5,'WindowStyle','docked');
plot(lon);
title('Longitude');

%%
% 4 bits Geoid Sign
% 16 bits Geoid BCD whole number
% 4 bits Elevation Sign
% 20 bits Elevation BCD whole number
% 20 bits Elevation BCD fraction
elev_wgs = double(bitand(elev,15));
elev = bitshift(elev,-4);
for idx = 1:9
  elev_wgs = elev_wgs/10 + double(bitand(elev,15));
  elev = bitshift(elev,-4);
end
elev_wgs = elev_wgs * 1e4;
elev_sign = double(bitand(elev,15));
elev = bitshift(elev,-4);

geoid = double(bitand(elev,15));
elev = bitshift(elev,-4);
for idx = 1:3
  geoid = geoid/10 + double(bitand(elev,15));
  elev = bitshift(elev,-4);
end
geoid = geoid * 100;
if elev
  geoid = -geoid;
else
  geoid = geoid;
end

if elev_sign
  elev = -elev_wgs;
else
  elev = elev_wgs;
end

figure(6); clf;
set(6,'WindowStyle','docked');
subplot(2,1,1);
plot(elev);
title('Elevation');
subplot(2,1,2);
plot(geoid);
title('Geoid');

return

% end
