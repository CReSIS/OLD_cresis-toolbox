% FIND TIME DELAY FOR KU BAND RADAR
% Load qlook data and lidar data, find the correct time delay
% set time delay to zero, and run this with find_correction = 1
% then use time_delay in "Time Delay, Td" field of radar worksheet




% % % ---------------------------------------------------------------------
% COMPLETE THIS SECTION - ALL VARIABLES ARE STRINGS
tic;
qlook_base_dir = '/cresis/scratch2/mdce/kuband/2010_Greenland_DC8/CSARP_qlook/';
year = '2010';
month = '04';
day = '13';
seg = '01'; 
frm = '200'; % 3 numbers go here



% % % ---------------------------------------------------------------------
% % % ---------------------------------------------------------------------

fprintf('\n\n============================================================\n');
fprintf('find_time_delay.m \n');
fprintf('------------------------------------------------------------\n');
fprintf('Loading data (%1.1f sec)\n',toc);
atm_fns = get_filenames('/cresis/data1/NASA/2010_Greenland_DC8_ATM/',...
    strcat(year(3:4),month,day),'','');
lidar = read_lidar_atm(atm_fns);
temp_str = strcat(year,month,day,'_',seg);
kuband = load(strcat(qlook_base_dir,temp_str,'/Data_',temp_str,'_',frm,'.mat'));

fprintf('Calculating and plotting results (%1.1f sec)\n',toc);
new_kuband_surface = interp1(kuband.GPS_time,kuband.Surface,lidar.gps_time);
new_kuband_elevation = interp1(kuband.GPS_time,kuband.Elevation,lidar.gps_time);
new_kuband_lat = interp1(kuband.GPS_time,kuband.Latitude,lidar.gps_time);
new_kuband_lon = interp1(kuband.GPS_time,kuband.Longitude,lidar.gps_time);
physical_constants;
offset_m = lidar.surface - (new_kuband_elevation - new_kuband_surface*c/2);
offset_m = offset_m(~isnan(offset_m));
ave_diff = mean(offset_m);
time_correction = ave_diff*2/c;


% % % ---------------------------------------------------------------------
% PLOT/CHECK THE RESULTS
[a,good_idxs] = find(~isnan(new_kuband_surface));

figure(1);clf;
subplot(121);plot(lidar.gps_time(good_idxs),(new_kuband_elevation(good_idxs)...
    - new_kuband_surface(good_idxs)*c/2),'r');
title('plot of surface from kuband data'); ylabel('elevation (m)');
subplot(122);plot(lidar.gps_time(good_idxs),lidar.surface(good_idxs),'b');
title('plot of surface from lidar data'); ylabel('elevation (m)');
figure(2);clf;plot(lidar.lon(good_idxs),lidar.lat(good_idxs),'vr');hold on;...
    plot(new_kuband_lon(good_idxs),new_kuband_lat(good_idxs),'.b');
    title('Lat vs. Lon of both lidar and kuband');legend('lidar','kuband');

fprintf('Done (%1.1f sec)\n',toc);
fprintf('\n Time Delay for %s (frm: %s): Td = %13.5e\n',temp_str,frm,time_correction);
fprintf('============================================================\n');


