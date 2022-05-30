  
%% Flightline Data Interpolation - 79N Line 1
% Years: 2010, 2014, 2016, 2017, 2018, 2019
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
AT_data.elevB.P2014 = (pass(2).layers(2).layer_elev);
AT_data.elevB.P2016 = (pass(3).layers(2).layer_elev);
AT_data.elevB.P2018 = (pass(4).layers(2).layer_elev);

%Save Surface profiles
AT_data.elevS.P2010 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2014 = (pass(2).layers(1).layer_elev);
AT_data.elevS.P2016 = (pass(3).layers(1).layer_elev);
AT_data.elevS.P2018 = (pass(4).layers(1).layer_elev);

%Save annual alongtrack profile data
AT_data.pass.P2010 = pass(1).along_track;
AT_data.pass.P2014 = pass(2).along_track;
AT_data.pass.P2016 = pass(3).along_track;
AT_data.pass.P2018 = pass(4).along_track;

%Save annual velocity correction data
AT_data.vel.P2010 = pass(1).vel;
AT_data.vel.P2014 = pass(2).vel;
AT_data.vel.P2016 = pass(3).vel;
AT_data.vel.P2018 = pass(4).vel;

%Save Velocity Corrected Along_track data
AT_data.AT_vel.P2010 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2014 = pass(baseline_master_idx).along_track + pass(2).vel;
AT_data.AT_vel.P2016 = pass(baseline_master_idx).along_track + pass(3).vel;
AT_data.AT_vel.P2018 = pass(baseline_master_idx).along_track + pass(4).vel;

% Latitudes and longitudes
AT_data.latitudes.P2010 = interp1(pass(1).lat, pass(1).lat, pass(2).lat);  
AT_data.latitudes.P2014 = pass(2).lat;
AT_data.latitudes.P2016 = interp1(pass(3).lat, pass(3).lat, pass(2).lat);
AT_data.latitudes.P2018 = interp1(pass(4).lat, pass(4).lat, pass(2).lat);
AT_data.longitudes.P2010 = interp1(pass(1).lon, pass(1).lon, pass(2).lon);
AT_data.longitudes.P2014 = pass(2).lon;
AT_data.longitudes.P2016 = interp1(pass(3).lon, pass(3).lon, pass(2).lon);
AT_data.longitudes.P2018 = interp1(pass(4).lon, pass(4).lon, pass(2).lon);

%% Section 2 - Mask, indexing, and clipping of along track
  
% Locate Along Track start element in each profile, save element ID as a
% variable for the clipping
% 2010
AT_data.Btrack.P10 = AT_data.AT_vel.P2010(AT_data.AT_vel.P2010 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P10 = find(AT_data.AT_vel.P2010 == ...
  AT_data.Btrack.P10(1));

% 2014  
AT_data.Btrack.P14 = AT_data.AT_vel.P2014(AT_data.AT_vel.P2014 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P14 = find(AT_data.AT_vel.P2014 == ...
  AT_data.Btrack.P14(1));

% 2016
AT_data.Btrack.P16 = AT_data.AT_vel.P2016(AT_data.AT_vel.P2016 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P16 = find(AT_data.AT_vel.P2016 == ...
  AT_data.Btrack.P16(1));
  
% 2018 
AT_data.Btrack.P18 = AT_data.AT_vel.P2018(AT_data.AT_vel.P2018 >= ...
  AT_data.AT_vel.P2010(1));
AT_data.find_AT_value.P18 = find(AT_data.AT_vel.P2018 == ... 
  AT_data.Btrack.P18(1));


% Clipping from start point in each profile to the end of the profile
AT_data.Btrack_Beg_Clip.P10 = AT_data.AT_vel.P2010...
  (AT_data.find_AT_value.P10:end);
AT_data.Btrack_Beg_Clip.P14 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end); 
AT_data.Btrack_Beg_Clip.P16 = AT_data.AT_vel.P2016...
  (AT_data.find_AT_value.P16:end);
AT_data.Btrack_Beg_Clip.P18 = AT_data.AT_vel.P2018...
  (AT_data.find_AT_value.P18:end);

% Clipping from new start locations to a given value end element value 
AT_data.Btrack_End_Clip.P10 = AT_data.Btrack_Beg_Clip.P10...
  (AT_data.Btrack_Beg_Clip.P10 <= 5.79e+04);
AT_data.Btrack_End_Clip.P14 = AT_data.Btrack_Beg_Clip.P14...
  (AT_data.Btrack_Beg_Clip.P14 <= 5.79e+04);
AT_data.Btrack_End_Clip.P16 = AT_data.Btrack_Beg_Clip.P16...
  (AT_data.Btrack_Beg_Clip.P16 <= 5.79e+04);
AT_data.Btrack_End_Clip.P18 = AT_data.Btrack_Beg_Clip.P18...
  (AT_data.Btrack_Beg_Clip.P18 <= 5.79e+04);

% Save along track data size as variable to see if there is any errors
AT_data.array_size.P10_AT = size(AT_data.Btrack_End_Clip.P10);
AT_data.array_size.P14_AT = size(AT_data.Btrack_End_Clip.P14);
AT_data.array_size.P16_AT = size(AT_data.Btrack_End_Clip.P16);
AT_data.array_size.P18_AT = size(AT_data.Btrack_End_Clip.P18);

%% Elevation data Clipping to Section size of Along Track files
% Elevation data beginning clipping from start element in Along Track
AT_data.elev_Beg_Clip.P2010 = AT_data.elevB.P2010...
  (AT_data.find_AT_value.P10:end); 
AT_data.elev_Beg_Clip.P2014 = AT_data.elevB.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip.P2016 = AT_data.elevB.P2016...
  (AT_data.find_AT_value.P16:end);
AT_data.elev_Beg_Clip.P2018 = AT_data.elevB.P2018...
  (AT_data.find_AT_value.P18:end);

% Elevation data end clipping from end of Along track data BED
AT_data.elev_End_Clip.P2010 = AT_data.elev_Beg_Clip.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip.P2014 = AT_data.elev_Beg_Clip.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip.P2016 = AT_data.elev_Beg_Clip.P2016...
  (1:length(AT_data.Btrack_End_Clip.P16));
AT_data.elev_End_Clip.P2018 = AT_data.elev_Beg_Clip.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18)); 

% Elevation data beginning clipping from start element in Along Track SURF
AT_data.elev_Beg_Clip_SURF.P2010 = AT_data.elevS.P2010...
  (AT_data.find_AT_value.P10:end);  
AT_data.elev_Beg_Clip_SURF.P2014 = AT_data.elevS.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_SURF.P2016 = AT_data.elevS.P2016...
  (AT_data.find_AT_value.P16:end);
AT_data.elev_Beg_Clip_SURF.P2018 = AT_data.elevS.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Elevation data end clipping from end of Along track data SURF
AT_data.elev_End_Clip_SURF.P2010 = AT_data.elev_Beg_Clip_SURF.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip_SURF.P2014 = AT_data.elev_Beg_Clip_SURF.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_SURF.P2016 = AT_data.elev_Beg_Clip_SURF.P2016...
  (1:length(AT_data.Btrack_End_Clip.P16));
AT_data.elev_End_Clip_SURF.P2018 = AT_data.elev_Beg_Clip_SURF.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Velocity data beginning clipping from start element in Along Track LAT
AT_data.elev_Beg_Clip_LAT.P2010 = AT_data.latitudes.P2010...
  (AT_data.find_AT_value.P10:end); 
AT_data.elev_Beg_Clip_LAT.P2014 = AT_data.latitudes.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_LAT.P2016 = AT_data.latitudes.P2016...
  (AT_data.find_AT_value.P16:end);
AT_data.elev_Beg_Clip_LAT.P2018 = AT_data.latitudes.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Velocity data end clipping from end of Along track data LAT
AT_data.elev_End_Clip_LAT.P2010 = AT_data.elev_Beg_Clip_LAT.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip_LAT.P2014 = AT_data.elev_Beg_Clip_LAT.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_LAT.P2016 = AT_data.elev_Beg_Clip_LAT.P2016...
  (1:length(AT_data.Btrack_End_Clip.P16));
AT_data.elev_End_Clip_LAT.P2018 = AT_data.elev_Beg_Clip_LAT.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Velocity data beginning clipping from start element in Along Track LON
AT_data.elev_Beg_Clip_LON.P2010 = AT_data.longitudes.P2010...
  (AT_data.find_AT_value.P10:end); 
AT_data.elev_Beg_Clip_LON.P2014 = AT_data.longitudes.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_LON.P2016 = AT_data.longitudes.P2016...
  (AT_data.find_AT_value.P16:end); 
AT_data.elev_Beg_Clip_LON.P2018 = AT_data.longitudes.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Velocity data end clipping from end of Along track data LON
AT_data.elev_End_Clip_LON.P2010 = AT_data.elev_Beg_Clip_LON.P2010...
  (1:length(AT_data.Btrack_End_Clip.P10));
AT_data.elev_End_Clip_LON.P2014 = AT_data.elev_Beg_Clip_LON.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_LON.P2016 = AT_data.elev_Beg_Clip_LON.P2016...
  (1:length(AT_data.Btrack_End_Clip.P16));
AT_data.elev_End_Clip_LON.P2018 = AT_data.elev_Beg_Clip_LON.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Save along elevation data size as variable to see if there is any errors
AT_data.array_size.P10_elev = size(AT_data.elev_End_Clip.P2010);
AT_data.array_size.P14_elev = size(AT_data.elev_End_Clip.P2014);
AT_data.array_size.P16_elev = size(AT_data.elev_End_Clip.P2016);
AT_data.array_size.P18_elev = size(AT_data.elev_End_Clip.P2018);
% SURF
AT_data.array_size.P10_elevSURF = size(AT_data.elev_End_Clip_SURF.P2010);
AT_data.array_size.P14_elevSURF = size(AT_data.elev_End_Clip_SURF.P2014);
AT_data.array_size.P16_elevSURF = size(AT_data.elev_End_Clip_SURF.P2016);
AT_data.array_size.P18_elevSURF = size(AT_data.elev_End_Clip_SURF.P2018);
% LAT
AT_data.array_size.P10_LAT = size(AT_data.elev_End_Clip_LAT.P2010);
AT_data.array_size.P14_LAT = size(AT_data.elev_End_Clip_LAT.P2014);
AT_data.array_size.P16_LAT = size(AT_data.elev_End_Clip_LAT.P2016);
AT_data.array_size.P18_LAT = size(AT_data.elev_End_Clip_LAT.P2018);
% LON
AT_data.array_size.P10_LON = size(AT_data.elev_End_Clip_LON.P2010);
AT_data.array_size.P14_LON = size(AT_data.elev_End_Clip_LON.P2014);
AT_data.array_size.P16_LON = size(AT_data.elev_End_Clip_LON.P2016);
AT_data.array_size.P18_LON = size(AT_data.elev_End_Clip_LON.P2018);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
AT_data.query_array.P10 = (AT_data.AT_vel.P2010(1):0.1:...
  AT_data.AT_vel.P2010(end));
%AT_data.query_array.P14 = (AT_data.AT_vel.P2014(1):0.1:...
%AT_data.AT_vel.P2014(end));
%AT_data.query_array.P16 = (AT_data.AT_vel.P2016(1):0.1:...
%AT_data.AT_vel.P2016(end));
%AT_data.query_array.P18 = (AT_data.AT_vel.P2018(1):0.1:...
%AT_data.AT_vel.P2018(end));

% Apply interpolation to each profile using selected query array
% 2010
AT_data.interp_data.P10 = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip.P2010, AT_data.query_array.P10);
% 2014
AT_data.interp_data.P14 = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip.P2014, AT_data.query_array.P10);
% 2016
AT_data.interp_data.P16 = interp1(AT_data.Btrack_End_Clip.P16, ...
  AT_data.elev_End_Clip.P2016, AT_data.query_array.P10);
% 2018
AT_data.interp_data.P18 = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip.P2018, AT_data.query_array.P10);

% 2010 surf
AT_data.interp_data.P10_Surf = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip_SURF.P2010, AT_data.query_array.P10);
% 2014 surf
AT_data.interp_data.P14_Surf = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_SURF.P2014, AT_data.query_array.P10);
% 2016 surf
AT_data.interp_data.P16_Surf = interp1(AT_data.Btrack_End_Clip.P16, ...
  AT_data.elev_End_Clip_SURF.P2016, AT_data.query_array.P10);
% 2018 surf
AT_data.interp_data.P18_Surf = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_SURF.P2018, AT_data.query_array.P10);

% 2010 LAT
AT_data.interp_data.P10_LAT = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip_LAT.P2010, AT_data.query_array.P10);
% 2014 LAT
AT_data.interp_data.P14_LAT = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LAT.P2014, AT_data.query_array.P10);
% 2016 LAT
AT_data.interp_data.P16_LAT = interp1(AT_data.Btrack_End_Clip.P16, ...
  AT_data.elev_End_Clip_LAT.P2016, AT_data.query_array.P10);
% 2018 LAT
AT_data.interp_data.P18_LAT = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_LAT.P2018, AT_data.query_array.P10);

% 2010 LON
AT_data.interp_data.P10_LON = interp1(AT_data.Btrack_End_Clip.P10, ...
  AT_data.elev_End_Clip_LON.P2010, AT_data.query_array.P10);
% 2014 LON
AT_data.interp_data.P14_LON = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LON.P2014, AT_data.query_array.P10);
% 2016 LON
AT_data.interp_data.P16_LON = interp1(AT_data.Btrack_End_Clip.P16, ...
  AT_data.elev_End_Clip_LON.P2016, AT_data.query_array.P10);
% 2018 LON
AT_data.interp_data.P18_LON = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_LON.P2018, AT_data.query_array.P10);

% Calculate melt rates from interpolated profile pairings
% 2010-2014 melt (Vertical difference in Features)
AT_data.melt_rates.P10_P14 = AT_data.interp_data.P14 - ...
  AT_data.interp_data.P10; 
  
% 2014-2016 melt (Vertical difference in Features)
AT_data.melt_rates.P14_P16 = AT_data.interp_data.P16 - ...
  AT_data.interp_data.P14;
  
% 2016-2018 melt (Vertical difference in Features)
AT_data.melt_rates.P16_P18 = AT_data.interp_data.P18 - ...
  AT_data.interp_data.P16;

% 2010-2018 melt (Vertical difference in Features)
AT_data.melt_rates.P10_P18 = AT_data.interp_data.P18 - ...
  AT_data.interp_data.P10;
  
%% Alignment of profiles test figures
figure(1)
h1 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip.P2010, 'color','#0072BD');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014, 'color','#D95319');
h3 = plot(AT_data.Btrack_End_Clip.P16/1e3, AT_data.elev_End_Clip.P2016, 'color','#EDB120');
h5 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018, 'color','#77AC30');
h6 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip_SURF.P2010, 'color','#0072BD');
h7 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip_SURF.P2014, 'color','#D95319');
h8 = plot(AT_data.Btrack_End_Clip.P16/1e3, AT_data.elev_End_Clip_SURF.P2016, 'color','#EDB120');
h9 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip_SURF.P2018, 'color','#77AC30');
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
xlim([10, 60]);
ylim([-600, 100]);
legend('original 2010', 'original 2014', 'original 2016', ...
    'original 2018', 'Location', 'southeast');

figure(2)
subplot(2,1,1)
h1 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip.P2010, 'color','#0072BD');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014, 'color','#D95319');
h3 = plot(AT_data.Btrack_End_Clip.P16/1e3, AT_data.elev_End_Clip.P2016, 'color','#EDB120');
h5 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018, 'color','#77AC30');
h6 = plot(AT_data.Btrack_End_Clip.P10/1e3, AT_data.elev_End_Clip_SURF.P2010, 'color','#0072BD');
h7 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip_SURF.P2014, 'color','#D95319');
h8 = plot(AT_data.Btrack_End_Clip.P16/1e3, AT_data.elev_End_Clip_SURF.P2016, 'color','#EDB120');
h9 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip_SURF.P2018, 'color','#77AC30');
title('79N - Clipped aligned profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
xlim([0, 60]);
ylim([-650, 150]);
legend('original 2010', 'original 2014', 'original 2016', ...
    'original 2018', 'Location', 'southeast');
subplot(2,1,2)
h1 = plot(AT_data.AT_vel.P2010/1e3, AT_data.elevB.P2010, 'color','#0072BD');
hold on
h2 = plot(AT_data.AT_vel.P2014/1e3, AT_data.elevB.P2014, 'color','#D95319');
h3 = plot(AT_data.AT_vel.P2016/1e3, AT_data.elevB.P2016, 'color','#EDB120');
h5 = plot(AT_data.AT_vel.P2018/1e3, AT_data.elevB.P2018, 'color','#77AC30');
h6 = plot(AT_data.AT_vel.P2010/1e3, AT_data.elevS.P2010, 'color','#0072BD');
h7 = plot(AT_data.AT_vel.P2014/1e3, AT_data.elevS.P2014, 'color','#D95319');
h8 = plot(AT_data.AT_vel.P2016/1e3, AT_data.elevS.P2016, 'color','#EDB120');
h9 = plot(AT_data.AT_vel.P2018/1e3, AT_data.elevS.P2018, 'color','#77AC30');
title('79N - Unclipped alignment profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
xlim([0, 75]);
ylim([-650, 150]);
legend('original 2010', 'original 2014', 'original 2016', ...
    'original 2018', 'Location', 'southeast');
%% Check interpolated data figures 
figure(3)
h1 = plot(AT_data.query_array.P10, AT_data.interp_data.P10);
hold on
h2 = plot(AT_data.query_array.P10, AT_data.interp_data.P14);
h3 = plot(AT_data.query_array.P10, AT_data.interp_data.P16);
h4 = plot(AT_data.query_array.P10, AT_data.interp_data.P18);
h5 = plot(AT_data.Btrack_End_Clip.P10, AT_data.elev_End_Clip.P2010);
h6 = plot(AT_data.Btrack_End_Clip.P14, AT_data.elev_End_Clip.P2014);
h7 = plot(AT_data.Btrack_End_Clip.P16, AT_data.elev_End_Clip.P2016);
h8 = plot(AT_data.Btrack_End_Clip.P18, AT_data.elev_End_Clip.P2018);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2010', 'interpolated 2014', 'interpolated 2016', ...
  'interpolated 2018', 'old 2010','old 2014','old 2016', 'old 2018', ...
  'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P10, AT_data.interp_data.P10);
hold on
h2 = plot(AT_data.query_array.P10, AT_data.interp_data.P14);
h3 = plot(AT_data.query_array.P10, AT_data.interp_data.P16);
h4 = plot(AT_data.query_array.P10, AT_data.interp_data.P18);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2010', 'interpolated 2014', 'interpolated 2016',...
  'interpolated 2018', 'Location', 'southeast');

%% Melt Rate Figures  
figure(4)
h1 = plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P14);
hold on
h2 = plot(AT_data.query_array.P10, AT_data.melt_rates.P14_P16);
h3 = plot(AT_data.query_array.P10, AT_data.melt_rates.P16_P18);
h4 = plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P18);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2010-2014 melt', '2014-2016 melt', '2016-2018 melt', ...
  '2010-2018 melt', 'Location', 'southeast');

% Subplot figure of melt rates
figure(5)
subplot(4,1,1);
plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P14);
title('79N Total Melt 2010-2014');

subplot(4,1,2);
plot(AT_data.query_array.P10, AT_data.melt_rates.P14_P16);
title('79N Total Melt 2014-2016');

subplot(4,1,3);
plot(AT_data.query_array.P10, AT_data.melt_rates.P16_P18);
title('79N Total Melt 2016-2018');

subplot(4,1,4);
plot(AT_data.query_array.P10, AT_data.melt_rates.P10_P18);
title('79N Total Melt 2010-2018');

% Annual melt rate average based off 2010-2019 melt subtraction
figure(6)
plot(AT_data.query_array.P10/1e3, AT_data.melt_rates.P10_P18/8, 'color','g');
title('79N Annual Average Melt 2010-2018');
ylim([-10, 50]);
xlabel('Along Track Distance (km)');
