  
%% Flightline Data Interpolation - Petermann Line 4
% Years: 2010, 2011, 2013, 2014
% Author: Cody Barnett
%
% Load elevation and along track data from run_combine_passes.m that has
% been velocity corrected by velocity_coregister.m into a new structure
% for comparison. 
% Section 1 - load data into new struct and apply velocity correction to
% each along track profile from pass struct in previous functions.
% Section 2 - Mask out regions where there is no overlap and then clip 
% both the along track and elevation data to the correct sizes, with
% excess data being removed for non-overlap regions between years
% Section 3 - create query points for the both the beginning and ends of
% the profiles using points from both sides of the cut off profiles. Query
% points and elevation data is then interpolated to match the dimensions
% of the compared year data. Melt rates as then derived from this new
% resampled sections. 
 
%Make AT_data structure and save Bed, Surface, Along_track data
AT_data = struct('pass', [], 'vel', [], 'elevB', [], 'elevS', [],...
'elev_pad', [], 'elev_NC',[], 'Btrack', []);
%Save Bottom profiles
AT_data.elevB.P2010 = (pass(1).layers(2).layer_elev);
AT_data.elevB.P2011 = (pass(2).layers(2).layer_elev);
AT_data.elevB.P2013 = (pass(3).layers(2).layer_elev);
AT_data.elevB.P2014 = (pass(4).layers(2).layer_elev);

%Save Surface profiles
AT_data.elevS.P2010 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2011 = (pass(2).layers(1).layer_elev);
AT_data.elevS.P2013 = (pass(3).layers(1).layer_elev);
AT_data.elevS.P2014 = (pass(4).layers(1).layer_elev);

%Save annual alongtrack profile data
AT_data.pass.P2010 = pass(1).along_track;
AT_data.pass.P2011 = pass(2).along_track;
AT_data.pass.P2013 = pass(3).along_track;
AT_data.pass.P2014 = pass(4).along_track;

%Save annual velocity correction data
AT_data.vel.P2010 = pass(1).vel;
AT_data.vel.P2011 = pass(2).vel;
AT_data.vel.P2013 = pass(3).vel;
AT_data.vel.P2014 = pass(4).vel;

%Save Velocity Corrected Along_track data
AT_data.AT_vel.P2010 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2011 = pass(baseline_master_idx).along_track + pass(2).vel;
AT_data.AT_vel.P2013 = pass(baseline_master_idx).along_track + pass(3).vel;
AT_data.AT_vel.P2014 = pass(baseline_master_idx).along_track + pass(4).vel;

% Latitudes and longitudes (Adjust for the master pass)
AT_data.latitudes.P2010 = pass(1).lat;  
AT_data.latitudes.P2011 = interp1(pass(2).lat, pass(2).lat, pass(1).lat);
AT_data.latitudes.P2013 = interp1(pass(3).lat, pass(3).lat, pass(1).lat);
AT_data.latitudes.P2014 = interp1(pass(4).lat, pass(4).lat, pass(1).lat);
AT_data.longitudes.P2010 = pass(1).lon;
AT_data.longitudes.P2011 = interp1(pass(2).lon, pass(2).lon, pass(1).lon);
AT_data.longitudes.P2013 = interp1(pass(3).lon, pass(3).lon, pass(1).lon);
AT_data.longitudes.P2014 = interp1(pass(4).lon, pass(4).lon, pass(1).lon);

%% Section 2 - Mask, indexing, and clipping of along track
  
% Locate Along Track start element in each profile, save element ID as a
% variable for the clipping
% 2010
AT_data.Btrack.P10 = AT_data.AT_vel.P2010(AT_data.AT_vel.P2010 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P10 = find(AT_data.AT_vel.P2010 == ...
  AT_data.Btrack.P10(1));

% 2011  
AT_data.Btrack.P11 = AT_data.AT_vel.P2011(AT_data.AT_vel.P2011 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P11 = find(AT_data.AT_vel.P2011 == ...
  AT_data.Btrack.P11(1));

% 2013
AT_data.Btrack.P13 = AT_data.AT_vel.P2013(AT_data.AT_vel.P2013 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P13 = find(AT_data.AT_vel.P2013 == ...
  AT_data.Btrack.P13(1));
  
% 2014 
AT_data.Btrack.P14 = AT_data.AT_vel.P2014(AT_data.AT_vel.P2014 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P14 = find(AT_data.AT_vel.P2014 == ... 
  AT_data.Btrack.P14(1));

% Clipping from start point in each profile to the end of the profile
AT_data.Btrack_Beg_Clip.P10 = AT_data.AT_vel.P2010...
  (AT_data.find_AT_value.P10:end);
AT_data.Btrack_Beg_Clip.P11 = AT_data.AT_vel.P2011...
  (AT_data.find_AT_value.P11:end); 
AT_data.Btrack_Beg_Clip.P13 = AT_data.AT_vel.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.Btrack_Beg_Clip.P14 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end);

% Clipping from new start locations to a given value end element value 
AT_data.Btrack_End_Clip.P10 = AT_data.Btrack_Beg_Clip.P10...
  (AT_data.Btrack_Beg_Clip.P10 <= 5.79e+04);
AT_data.Btrack_End_Clip.P11 = AT_data.Btrack_Beg_Clip.P11...
  (AT_data.Btrack_Beg_Clip.P11 <= 5.79e+04);
AT_data.Btrack_End_Clip.P13 = AT_data.Btrack_Beg_Clip.P13...
  (AT_data.Btrack_Beg_Clip.P13 <= 5.79e+04);
AT_data.Btrack_End_Clip.P14 = AT_data.Btrack_Beg_Clip.P14...
  (AT_data.Btrack_Beg_Clip.P14 <= 5.79e+04); 

% Save along track data size as variable to see if there is any errors
AT_data.array_size.P10_AT = size(AT_data.Btrack_End_Clip.P10);
AT_data.array_size.P11_AT = size(AT_data.Btrack_End_Clip.P11);
AT_data.array_size.P13_AT = size(AT_data.Btrack_End_Clip.P13);
AT_data.array_size.P14_AT = size(AT_data.Btrack_End_Clip.P14);

%% Elevation data Clipping to Section size of Along Track files
% Elevation data beginning clipping from start element in Along Track BED
AT_data.elev_Beg_Clip.P2010 = AT_data.elevB.P2010...
  (AT_data.find_AT_value.P10:end); 
AT_data.elev_Beg_Clip.P2011 = AT_data.elevB.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip.P2013 = AT_data.elevB.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip.P2014 = AT_data.elevB.P2014...
  (AT_data.find_AT_value.P14:end);

% Elevation data end clipping from end of Along track data BED
AT_data.elev_End_Clip.P2010 = AT_data.elev_Beg_Clip.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip.P2011 = AT_data.elev_Beg_Clip.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip.P2013 = AT_data.elev_Beg_Clip.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip.P2014 = AT_data.elev_Beg_Clip.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));  

% Elevation data beginning clipping from start element in Along Track SURF
AT_data.elev_Beg_Clip_SURF.P2010 = AT_data.elevS.P2010...
  (AT_data.find_AT_value.P10:end);  
AT_data.elev_Beg_Clip_SURF.P2011 = AT_data.elevS.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip_SURF.P2013 = AT_data.elevS.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_SURF.P2014 = AT_data.elevS.P2014...
  (AT_data.find_AT_value.P14:end);
  
% Elevation data end clipping from end of Along track data SURF
AT_data.elev_End_Clip_SURF.P2010 = AT_data.elev_Beg_Clip_SURF.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip_SURF.P2011 = AT_data.elev_Beg_Clip_SURF.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip_SURF.P2013 = AT_data.elev_Beg_Clip_SURF.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_SURF.P2014 = AT_data.elev_Beg_Clip_SURF.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));

% Velocity data beginning clipping from start element in Along Track LAT
AT_data.elev_Beg_Clip_LAT.P2010 = AT_data.latitudes.P2010...
  (AT_data.find_AT_value.P10:end); 
AT_data.elev_Beg_Clip_LAT.P2011 = AT_data.latitudes.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip_LAT.P2013 = AT_data.latitudes.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_LAT.P2014 = AT_data.latitudes.P2014...
  (AT_data.find_AT_value.P14:end);
  
% Velocity data end clipping from end of Along track data LAT
AT_data.elev_End_Clip_LAT.P2010 = AT_data.elev_Beg_Clip_LAT.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip_LAT.P2011 = AT_data.elev_Beg_Clip_LAT.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip_LAT.P2013 = AT_data.elev_Beg_Clip_LAT.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_LAT.P2014 = AT_data.elev_Beg_Clip_LAT.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));

% Velocity data beginning clipping from start element in Along Track LON
AT_data.elev_Beg_Clip_LON.P2010 = AT_data.longitudes.P2010...
  (AT_data.find_AT_value.P10:end); 
AT_data.elev_Beg_Clip_LON.P2011 = AT_data.longitudes.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip_LON.P2013 = AT_data.longitudes.P2013...
  (AT_data.find_AT_value.P13:end); 
AT_data.elev_Beg_Clip_LON.P2014 = AT_data.longitudes.P2014...
  (AT_data.find_AT_value.P14:end);
  
% Velocity data end clipping from end of Along track data LON
AT_data.elev_End_Clip_LON.P2010 = AT_data.elev_Beg_Clip_LON.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip_LON.P2011 = AT_data.elev_Beg_Clip_LON.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip_LON.P2013 = AT_data.elev_Beg_Clip_LON.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_LON.P2014 = AT_data.elev_Beg_Clip_LON.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));

% Save along elevation data size as variable to see if there is any errors
% BED
AT_data.array_size.P10_elev = size(AT_data.elev_End_Clip.P2010);
AT_data.array_size.P11_elev = size(AT_data.elev_End_Clip.P2011);
AT_data.array_size.P13_elev = size(AT_data.elev_End_Clip.P2013);
AT_data.array_size.P14_elev = size(AT_data.elev_End_Clip.P2014);
% SURF
AT_data.array_size.P10_elevSURF = size(AT_data.elev_End_Clip_SURF.P2010);
AT_data.array_size.P11_elevSURF = size(AT_data.elev_End_Clip_SURF.P2011);
AT_data.array_size.P13_elevSURF = size(AT_data.elev_End_Clip_SURF.P2013);
AT_data.array_size.P14_elevSURF = size(AT_data.elev_End_Clip_SURF.P2014);
% LAT
AT_data.array_size.P10_LAT = size(AT_data.elev_End_Clip_LAT.P2010);
AT_data.array_size.P11_LAT = size(AT_data.elev_End_Clip_LAT.P2011);
AT_data.array_size.P13_LAT = size(AT_data.elev_End_Clip_LAT.P2013);
AT_data.array_size.P14_LAT = size(AT_data.elev_End_Clip_LAT.P2014);
% LON
AT_data.array_size.P10_LON = size(AT_data.elev_End_Clip_LON.P2010);
AT_data.array_size.P11_LON = size(AT_data.elev_End_Clip_LON.P2011);
AT_data.array_size.P13_LON = size(AT_data.elev_End_Clip_LON.P2013);
AT_data.array_size.P14_LON = size(AT_data.elev_End_Clip_LON.P2014);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
AT_data.query_array.P10 = (AT_data.AT_vel.P2010(1):10:...
  AT_data.AT_vel.P2010(end));
%AT_data.query_array.P11 = (AT_data.AT_vel.P2011(1):0.1:...
%AT_data.AT_vel.P2011(end));
%AT_data.query_array.P13 = (AT_data.AT_vel.P2013(1):0.1:...
%AT_data.AT_vel.P2013(end));
%AT_data.query_array.P14 = (AT_data.AT_vel.P2014(1):0.1:...
%AT_data.AT_vel.P2014(end));

% Apply interpolation to each profile using selected query array
% 2010 Bed
AT_data.interp_data.P10 = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip.P2010, AT_data.query_array.P10);
% 2011 Bed
AT_data.interp_data.P11 = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip.P2011, AT_data.query_array.P10);
% 2013 Bed
AT_data.interp_data.P13 = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip.P2013, AT_data.query_array.P10);
% 2014 Bed
AT_data.interp_data.P14 = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip.P2014, AT_data.query_array.P10);

% 2010 surf
AT_data.interp_data.P10_Surf = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip_SURF.P2010, AT_data.query_array.P10);
% 2011 surf
AT_data.interp_data.P11_Surf = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip_SURF.P2011, AT_data.query_array.P10);
% 2013 surf
AT_data.interp_data.P13_Surf = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_SURF.P2013, AT_data.query_array.P10);
% 2014 surf
AT_data.interp_data.P14_Surf = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_SURF.P2014, AT_data.query_array.P10);

% 2010 LAT
AT_data.interp_data.P10_LAT = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip_LAT.P2010, AT_data.query_array.P10);
% 2011 LAT
AT_data.interp_data.P11_LAT = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip_LAT.P2011, AT_data.query_array.P10);
% 2013 LAT
AT_data.interp_data.P13_LAT = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_LAT.P2013, AT_data.query_array.P10);
% 2014 LAT
AT_data.interp_data.P14_LAT = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LAT.P2014, AT_data.query_array.P10);

% 2010 LON
AT_data.interp_data.P10_LON = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip_LON.P2010, AT_data.query_array.P10);
% 2011 LON
AT_data.interp_data.P11_LON = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip_LON.P2011, AT_data.query_array.P10);
% 2013 LON
AT_data.interp_data.P13_LON = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_LON.P2013, AT_data.query_array.P10);
% 2014 LON
AT_data.interp_data.P14_LON = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LON.P2014, AT_data.query_array.P10);

% Calculate melt rates from interpolated profile pairings
% 2010-2011 melt (Vertical difference in Features)
AT_data.melt_rates.P10_P11 = AT_data.interp_data.P11 - ...
  AT_data.interp_data.P10;   
% 2011-2013 melt (Vertical difference in Features)
AT_data.melt_rates.P11_P13 = AT_data.interp_data.P13 - ...
  AT_data.interp_data.P11; 
% 2013-2014 melt (Vertical difference in Features)
AT_data.melt_rates.P13_P14 = AT_data.interp_data.P14 - ...
  AT_data.interp_data.P13;
% 2010-2014 melt (Vertical difference in Features)
AT_data.melt_rates.P10_P14 = AT_data.interp_data.P14 - ...
  AT_data.interp_data.P10;

%% Export data to csv
% Concatenate and take transpose of Lon, Lat, Surf, Bed fields. Concatenate horizontally for each year
AT_data.export.P10 = cat(2, AT_data.interp_data.P10_LON.', AT_data.interp_data.P10_LAT.', AT_data.interp_data.P10_Surf.', AT_data.interp_data.P10.' );
AT_data.export.P11 = cat(2, AT_data.interp_data.P11_LON.', AT_data.interp_data.P11_LAT.', AT_data.interp_data.P11_Surf.', AT_data.interp_data.P11.' );
AT_data.export.P13 = cat(2, AT_data.interp_data.P13_LON.', AT_data.interp_data.P13_LAT.', AT_data.interp_data.P13_Surf.', AT_data.interp_data.P13.' );
AT_data.export.P14 = cat(2, AT_data.interp_data.P14_LON.', AT_data.interp_data.P14_LAT.', AT_data.interp_data.P14_Surf.', AT_data.interp_data.P14.' );

%% Define Header Array of strings and vertically concatenate to data 
cheader = {'Lons', 'Lats', 'Surface', 'Depth'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

% change folder
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\CSV_export_files\'

%write header to file 2011
fid = fopen('P4_lat_lon_surf_bed_10.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P4_lat_lon_surf_bed_10.csv', AT_data.export.P10, '-append');

%write header to file 2014
fid = fopen('P4_lat_lon_surf_bed_11.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P4_lat_lon_surf_bed_11.csv', AT_data.export.P11, '-append');

%write header to file 2015
fid = fopen('P4_lat_lon_surf_bed_13.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P4_lat_lon_surf_bed_13.csv', AT_data.export.P13, '-append');

%write header to file 2017
fid = fopen('P4_lat_lon_surf_bed_14.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P4_lat_lon_surf_bed_14.csv', AT_data.export.P14, '-append');

%% Test Alignment Figures Section 
figure(1)
h1 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip.P2010);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
h3 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.elev_End_Clip.P2013);
h4 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
h5 = plot(AT_data.Btrack_End_Clip.P10(C(1))/1e3, AT_data.elev_End_Clip.P2010(C(1)),'+','Markersize',50);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2010', 'original 2011', 'original 2013', ...
  'original 2014', 'Location', 'southeast');
%%
figure(2)
subplot(2,1,1)
h1 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip.P2010, 'color','#0072BD');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011, 'color','#D95319');
h3 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.elev_End_Clip.P2013, 'color','#EDB120');
h4 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014, 'color','#77AC30');
h5 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip_SURF.P2010, 'color','#0072BD');
h6 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip_SURF.P2011, 'color','#D95319');
h7 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.elev_End_Clip_SURF.P2013, 'color','#EDB120');
h8 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip_SURF.P2014, 'color','#77AC30');
title('79N - Clipped aligned profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
xlim([0, 60]);
ylim([-600, 150]);
legend('original 2010', 'original 2011', 'original 2013', ...
    'original 2014', 'Location', 'southeast');
subplot(2,1,2)
h1 = plot(AT_data.AT_vel.P2010/1e3, AT_data.elevB.P2010, 'color','#0072BD');
hold on
h2 = plot(AT_data.AT_vel.P2011/1e3, AT_data.elevB.P2011, 'color','#D95319');
h3 = plot(AT_data.AT_vel.P2013/1e3, AT_data.elevB.P2013, 'color','#EDB120');
h4 = plot(AT_data.AT_vel.P2014/1e3, AT_data.elevB.P2014, 'color','#77AC30');
h5 = plot(AT_data.AT_vel.P2010/1e3, AT_data.elevS.P2010, 'color','#0072BD');
h6 = plot(AT_data.AT_vel.P2011/1e3, AT_data.elevS.P2011, 'color','#D95319');
h7 = plot(AT_data.AT_vel.P2013/1e3, AT_data.elevS.P2013, 'color','#EDB120');
h8 = plot(AT_data.AT_vel.P2014/1e3, AT_data.elevS.P2014, 'color','#77AC30');
title('79N - Unclipped alignment profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
xlim([0, 65]);
ylim([-600, 150]);
legend('original 2010', 'original 2011', 'original 2013', ...
    'original 2014', 'Location', 'southeast');

figure(3)
h1 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevB.P2010, 'color','#0072BD');
hold on
h2 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevB.P2011, 'color','#D95319');
h3 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevB.P2013, 'color','#EDB120');
h4 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevB.P2014, 'color','#77AC30');
h5 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevS.P2010, 'color','#0072BD');
h6 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevS.P2011, 'color','#D95319');
h7 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevS.P2013, 'color','#EDB120');
h8 = plot(pass(baseline_master_idx).along_track/1e3, AT_data.elevS.P2014, 'color','#77AC30');
title('79N - Clipped aligned profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
xlim([0, 60]);
ylim([-600, 150]);
legend('original 2010', 'original 2011', 'original 2013', ...
    'original 2014', 'Location', 'southeast');

figure(4)
h3 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip.P2010, 'color','#EDB120');
hold on 
h7 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip_SURF.P2010, 'color','#EDB120');
title('2010');

figure(5)
subplot(2,1,1)
h1 = plot(AT_data.AT_vel.P2011/1e3, AT_data.elevB.P2011, 'color','#77AC30');
hold on 
h2 = plot(AT_data.AT_vel.P2011/1e3, AT_data.elevS.P2011, 'color','#77AC30');
title('2011');
subplot(2,1,2)
h3 = plot(AT_data.AT_vel.P2010/1e3, AT_data.elevB.P2010, 'color','#77AC30');
hold on 
h4 = plot(AT_data.AT_vel.P2010/1e3, AT_data.elevS.P2010, 'color','#77AC30');
title('2010');


%% Interpolated data Figures
figure(2)
h1 = plot(AT_data.query_array.P10, AT_data.interp_data.P10);
hold on
h2 = plot(AT_data.query_array.P10, AT_data.interp_data.P11);
h3 = plot(AT_data.query_array.P10, AT_data.interp_data.P13);
h4 = plot(AT_data.query_array.P10, AT_data.interp_data.P14);
h5 = plot(AT_data.Btrack_End_Clip.P10, AT_data.elev_End_Clip.P2010);
h6 = plot(AT_data.Btrack_End_Clip.P11, AT_data.elev_End_Clip.P2011);
h7 = plot(AT_data.Btrack_End_Clip.P13, AT_data.elev_End_Clip.P2013);
h8 = plot(AT_data.Btrack_End_Clip.P14, AT_data.elev_End_Clip.P2014);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2010', 'interpolated 2011', 'interpolated 2013', ...
  'interpolated 2014', 'old 2010','old 2011', 'old 2013', 'old 2014',...
  'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P10, AT_data.interp_data.P10);
hold on
h2 = plot(AT_data.query_array.P10, AT_data.interp_data.P11);
h3 = plot(AT_data.query_array.P10, AT_data.interp_data.P13);
h4 = plot(AT_data.query_array.P10, AT_data.interp_data.P14);
h5 = plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P14);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2010', 'interpolated 2011', 'interpolated 2013',...
  'interpolated 2014', '2010-2014 melt', ...
  'Location', 'southeast');

%% Melt Rate Figures  
figure(4)
h1 = plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P11);
hold on
h2 = plot(AT_data.query_array.P10, AT_data.melt_rates.P11_P13);
h3 = plot(AT_data.query_array.P10, AT_data.melt_rates.P13_P14);
h4 = plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P14);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2010-2011 melt', '2011-2013 melt', '2013-2014 melt', ...
  '2010-2014 melt', 'Location', 'southeast');

% Subplot figure of melt rates
figure(5)
subplot(4,1,1);
plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P11);
title('Melt 2010-2011');

subplot(4,1,2);
plot(AT_data.query_array.P10, AT_data.melt_rates.P11_P13);
title('Melt 2011-2013');

subplot(4,1,3);
plot(AT_data.query_array.P10, AT_data.melt_rates.P13_P14);
title('Melt 2013-2014');

subplot(4,1,4);
plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P14);
title('Melt 2014-2017');


% Annual melt rate average based off 2010-2019 melt subtraction
% figure(6)
% plot(AT_data.query_array.P10/1e3, AT_data.melt_rates.P10_P14/4, 'color','g');
% title('Melt 2013-2017');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');
