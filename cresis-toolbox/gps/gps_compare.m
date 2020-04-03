% script gps_compare
%
% Simple script to compare two GPS files
%
% Author: John Paden

day_str = '20190403';
gps = load(sprintf('/scratch/csarp_support/gps/2019_Greenland_P3/gps_%s.mat',day_str));
gps2 = load(sprintf('/scratch/csarp_support/gps/2019_Greenland_P3/gps_%s_novatel.mat',day_str));
out_dir = '~/';
out_types = {'fig','jpg'};

plot(angle(exp(j*interp1(gps.gps_time,gps.heading,gps2.gps_time))./exp(j*gps2.heading))*180/pi)
hfig = get_figures(1,1);
xlabel('Record');
ylabel('Heading delta (deg)');
grid on;
for idx=1:length(out_types)
  fn = fullfile(out_dir,sprintf('compare_gps_heading_%s.%s',day_str,out_types{idx}));
  fprintf('Saving %s\n', fn);
  saveas(hfig,fn);
end

plot(angle(exp(j*interp1(gps.gps_time,gps.roll,gps2.gps_time))./exp(j*gps2.roll))*180/pi)
xlabel('Record');
ylabel('Roll delta (deg)');
grid on;
for idx=1:length(out_types)
  fn = fullfile(out_dir,sprintf('compare_gps_roll_%s.%s',day_str,out_types{idx}));
  fprintf('Saving %s\n', fn);
  saveas(hfig,fn);
end

plot(angle(exp(j*interp1(gps.gps_time,gps.pitch,gps2.gps_time))./exp(j*gps2.pitch))*180/pi)
xlabel('Record');
ylabel('Pitch delta (deg)');
grid on;
for idx=1:length(out_types)
  fn = fullfile(out_dir,sprintf('compare_gps_pitch_%s.%s',day_str,out_types{idx}));
  fprintf('Saving %s\n', fn);
  saveas(hfig,fn);
end

plot(interp1(gps.gps_time,gps.elev,gps2.gps_time)-gps2.elev)
xlabel('Record');
ylabel('Elevation (m)');
grid on;
for idx=1:length(out_types)
  fn = fullfile(out_dir,sprintf('compare_gps_elev_%s.%s',day_str,out_types{idx}));
  fprintf('Saving %s\n', fn);
  saveas(hfig,fn);
end

plot(interp1(gps.gps_time,gps.lat,gps2.gps_time)-gps2.lat)
xlabel('Record');
ylabel('Latitude (deg)');
grid on;
for idx=1:length(out_types)
  fn = fullfile(out_dir,sprintf('compare_gps_lat_%s.%s',day_str,out_types{idx}));
  fprintf('Saving %s\n', fn);
  saveas(hfig,fn);
end

plot(interp1(gps.gps_time,gps.lon,gps2.gps_time)-gps2.lon)
xlabel('Record');
ylabel('Longitude (deg)');
grid on;
for idx=1:length(out_types)
  fn = fullfile(out_dir,sprintf('compare_gps_lon_%s.%s',day_str,out_types{idx}));
  fprintf('Saving %s\n', fn);
  saveas(hfig,fn);
end
