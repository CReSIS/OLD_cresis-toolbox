%% Flightline Data Interpolation - Petermann Line 1 - 2014, 2017, 2018
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

% Make AT_data structure and save Bed, Surface, Along_track data
AT_data = struct('pass', [], 'vel', [], 'elevB', [], 'elevS', [],...
'elev_pad', [], 'elev_NC',[], 'Btrack', []);

% Save Bottom profiles
AT_data.elevB.P2014 = (pass(1).layers(2).layer_elev);
AT_data.elevB.P2017 = (pass(2).layers(2).layer_elev);
AT_data.elevB.P2018 = (pass(3).layers(2).layer_elev);

% Save Surface profiles
AT_data.elevS.P2014 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2017 = (pass(2).layers(1).layer_elev);
AT_data.elevS.P2018 = (pass(3).layers(1).layer_elev);

% Save annual alongtrack profile data
AT_data.pass.P2014 = pass(1).along_track;
AT_data.pass.P2017 = pass(2).along_track;
AT_data.pass.P2018 = pass(3).along_track;
 
% Save annual velocity correction data
AT_data.vel.P2014 = pass(1).vel;
AT_data.vel.P2017 = pass(2).vel;
AT_data.vel.P2018 = pass(3).vel;

% Save Velocity Corrected Along_track data
AT_data.AT_vel.P2014 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2017 = pass(baseline_master_idx).along_track + pass(2).vel;
AT_data.AT_vel.P2018 = pass(baseline_master_idx).along_track + pass(3).vel;

% Distance difference
AT_data.dist_14_17 = AT_data.AT_vel.P2014 - AT_data.AT_vel.P2017;
AT_data.dist_17_18 = AT_data.AT_vel.P2017 - AT_data.AT_vel.P2018;
AT_data.dist_14_18 = AT_data.AT_vel.P2014 - AT_data.AT_vel.P2018;

% Latitudes and longitudes
AT_data.latitudes.P2014 = interp1(pass(1).lat, pass(1).lat, pass(2).lat);  
AT_data.latitudes.P2017 = pass(2).lat;
AT_data.latitudes.P2018 = interp1(pass(3).lat, pass(3).lat, pass(2).lat);
AT_data.longitudes.P2014 = interp1(pass(1).lon, pass(1).lon, pass(2).lon);
AT_data.longitudes.P2017 = pass(2).lon;
AT_data.longitudes.P2018 = interp1(pass(3).lon, pass(3).lon, pass(2).lon);


figure(1)
h1 = plot(AT_data.pass.P2014/1e3, AT_data.elevB.P2014);
hold on
h2 = plot(AT_data.pass.P2017/1e3, AT_data.elevB.P2017);
h3 = plot(AT_data.pass.P2018/1e3, AT_data.elevB.P2018);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2014', 'original 2017', 'original 2018', ...
  'Location', 'southeast');
  

%% Section 2 - Mask, indexing, and clipping of along track
  
% Locate Along Track start element in each profile, save element ID as a
% variable for the clipping
% 2014
AT_data.Btrack.P14 = AT_data.AT_vel.P2014(AT_data.AT_vel.P2014 >= ...
  AT_data.AT_vel.P2014(1));
AT_data.find_AT_value.P14 = find(AT_data.AT_vel.P2014 == ...
  AT_data.Btrack.P14(1));
  
% 2017  
AT_data.Btrack.P17 = AT_data.AT_vel.P2017(AT_data.AT_vel.P2017 >= ...
  AT_data.AT_vel.P2014(1));
AT_data.find_AT_value.P17 = find(AT_data.AT_vel.P2017 == ...
  AT_data.Btrack.P17(1));
  
% 2018 
AT_data.Btrack.P18 = AT_data.AT_vel.P2018(AT_data.AT_vel.P2018 >= ...
  AT_data.AT_vel.P2014(1));
AT_data.find_AT_value.P18 = find(AT_data.AT_vel.P2018 == ...
  AT_data.Btrack.P18(1));
  
% Clipping from start point in each profile to the end of the profile
AT_data.Btrack_Beg_Clip.P14 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.Btrack_Beg_Clip.P17 = AT_data.AT_vel.P2017...
  (AT_data.find_AT_value.P17:end); 
AT_data.Btrack_Beg_Clip.P18 = AT_data.AT_vel.P2018...
  (AT_data.find_AT_value.P18:end);

% Clipping from new start locations to a given value end element value 
AT_data.Btrack_End_Clip.P14 = AT_data.Btrack_Beg_Clip.P14...
  (AT_data.Btrack_Beg_Clip.P14 <= 5.79e+04);
AT_data.Btrack_End_Clip.P17 = AT_data.Btrack_Beg_Clip.P17...
  (AT_data.Btrack_Beg_Clip.P17 <= 5.79e+04);
AT_data.Btrack_End_Clip.P18 = AT_data.Btrack_Beg_Clip.P18...
  (AT_data.Btrack_Beg_Clip.P18 <= 5.79e+04);
  
%Save along track data size as variable to see if there is any errors
AT_data.array_size.P14_AT = size(AT_data.Btrack_End_Clip.P14);
AT_data.array_size.P17_AT = size(AT_data.Btrack_End_Clip.P17);
AT_data.array_size.P18_AT = size(AT_data.Btrack_End_Clip.P18);
%% For terminus End element trimming (NEEDS WORK)
%AT_data.Btrack_E.P14 = AT_data.Btrack_C.P14(AT_data.Btrack_C.P14 <= ...
%AT_data.Btrack_C.P14(end));
%AT_data.find_AT_value_E.P14 = find(AT_data.Btrack_C.P14 == ...
%AT_data.Btrack_C.P14(end));
%     
%AT_data.Btrack_E.P17 = AT_data.Btrack_C.P17(AT_data.Btrack_C.P17 <= ...
%AT_data.Btrack_C.P14(end));
%AT_data.find_AT_value_E.P17 = find(AT_data.Btrack_C.P17 == ...
%AT_data.Btrack_C.P14(end));
%  
%AT_data.Btrack_E.P18 = AT_data.Btrack_C.P18(AT_data.Btrack_C.P18 <= ...
%AT_data.Btrack_C.P14(end));
%AT_data.find_AT_value_E.P18 = find(AT_data.Btrack_C.P18 == ...
%AT_data.Btrack_C.P14(end));
  
%AT_data.Btrack_CE.P14 = AT_data.Btrack_C.P14...
%(1:AT_data.find_AT_value_E.P14);
%AT_data.Btrack_CE.P17 = AT_data.Btrack_C.P17...
%(1:AT_data.find_AT_value_E.P17); 
%AT_data.Btrack_CE.P18 = AT_data.Btrack_C.P18...
%(1:AT_data.find_AT_value_E.P18);
%   
%AT_data.Etrack_C.P14 = AT_data.Btrack_C.P14(1:length(AT_data.Btrack_C.P18));
%AT_data.Etrack_C.P17 = AT_data.Btrack_C.P17(1:length(AT_data.Btrack_C.P18));
%AT_data.Etrack_C.P18 = AT_data.Btrack_C.P18(1:length(AT_data.Btrack_C.P18));

% Along track end clipping mask to the shortest profile (physically)
%AT_data.Etrack.P14 = AT_data.Btrack.P14(AT_data.Btrack.P14 ...
%<= AT_data.Btrack.P17(end));
%AT_data.Etrack.P17 = AT_data.Btrack.P17(AT_data.Btrack.P17 ...
%<= AT_data.Btrack.P17(end));
%AT_data.Etrack.P18 = AT_data.Btrack.P18(AT_data.Btrack.P18 ...
%<= AT_data.Btrack.P17(end));

%% Elevation data Clipping to Section size of Along Track files
% Elevation data beginning clipping from start element in Along Track BED
AT_data.elev_Beg_Clip.P2014 = AT_data.elevB.P2014...
  (AT_data.find_AT_value.P14:end); 
AT_data.elev_Beg_Clip.P2017 = AT_data.elevB.P2017...
  (AT_data.find_AT_value.P17:end);
AT_data.elev_Beg_Clip.P2018 = AT_data.elevB.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Elevation data end clipping from end of Along track data BED
AT_data.elev_End_Clip.P2014 = AT_data.elev_Beg_Clip.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip.P2017 = AT_data.elev_Beg_Clip.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));
AT_data.elev_End_Clip.P2018 = AT_data.elev_Beg_Clip.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Elevation data beginning clipping from start element in Along Track SURF
AT_data.elev_Beg_Clip_SURF.P2014 = AT_data.elevS.P2014...
  (AT_data.find_AT_value.P14:end); 
AT_data.elev_Beg_Clip_SURF.P2017 = AT_data.elevS.P2017...
  (AT_data.find_AT_value.P17:end);
AT_data.elev_Beg_Clip_SURF.P2018 = AT_data.elevS.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Elevation data end clipping from end of Along track data SURF
AT_data.elev_End_Clip_SURF.P2014 = AT_data.elev_Beg_Clip_SURF.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_SURF.P2017 = AT_data.elev_Beg_Clip_SURF.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));
AT_data.elev_End_Clip_SURF.P2018 = AT_data.elev_Beg_Clip_SURF.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Velocity data beginning clipping from start element in Along Track VEL
AT_data.elev_Beg_Clip_VEL.P2014 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end); 
AT_data.elev_Beg_Clip_VEL.P2017 = AT_data.AT_vel.P2017...
  (AT_data.find_AT_value.P17:end);
AT_data.elev_Beg_Clip_VEL.P2018 = AT_data.AT_vel.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Velocity data end clipping from end of Along track data VEL
AT_data.elev_End_Clip_VEL.P2014 = AT_data.elev_Beg_Clip_VEL.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_VEL.P2017 = AT_data.elev_Beg_Clip_VEL.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));
AT_data.elev_End_Clip_VEL.P2018 = AT_data.elev_Beg_Clip_VEL.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Velocity data beginning clipping from start element in Along Track LAT
AT_data.elev_Beg_Clip_LAT.P2014 = AT_data.latitudes.P2014...
  (AT_data.find_AT_value.P14:end); 
AT_data.elev_Beg_Clip_LAT.P2017 = AT_data.latitudes.P2017...
  (AT_data.find_AT_value.P17:end);
AT_data.elev_Beg_Clip_LAT.P2018 = AT_data.latitudes.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Velocity data end clipping from end of Along track data LAT
AT_data.elev_End_Clip_LAT.P2014 = AT_data.elev_Beg_Clip_LAT.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_LAT.P2017 = AT_data.elev_Beg_Clip_LAT.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));
AT_data.elev_End_Clip_LAT.P2018 = AT_data.elev_Beg_Clip_LAT.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Velocity data beginning clipping from start element in Along Track LON
AT_data.elev_Beg_Clip_LON.P2014 = AT_data.longitudes.P2014...
  (AT_data.find_AT_value.P14:end); 
AT_data.elev_Beg_Clip_LON.P2017 = AT_data.longitudes.P2017...
  (AT_data.find_AT_value.P17:end);
AT_data.elev_Beg_Clip_LON.P2018 = AT_data.longitudes.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Velocity data end clipping from end of Along track data LON
AT_data.elev_End_Clip_LON.P2014 = AT_data.elev_Beg_Clip_LON.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_LON.P2017 = AT_data.elev_Beg_Clip_LON.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));
AT_data.elev_End_Clip_LON.P2018 = AT_data.elev_Beg_Clip_LON.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Velocity data beginning clipping from start element in Along Track DIST
AT_data.elev_Beg_Clip_DIST.P2014 = AT_data.dist_14_17...
  (AT_data.find_AT_value.P14:end); 
AT_data.elev_Beg_Clip_DIST.P2017 = AT_data.dist_17_18...
  (AT_data.find_AT_value.P17:end);
AT_data.elev_Beg_Clip_DIST.P2018 = AT_data.dist_14_18...
  (AT_data.find_AT_value.P18:end);
  
% Velocity data end clipping from end of Along track data DIST
AT_data.elev_End_Clip_DIST.P2014 = AT_data.elev_Beg_Clip_DIST.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_DIST.P2017 = AT_data.elev_Beg_Clip_DIST.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));
AT_data.elev_End_Clip_DIST.P2018 = AT_data.elev_Beg_Clip_DIST.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));

% Save along elevation data size as variable to see if there is any errors
% BED
AT_data.array_size.P14_elev = size(AT_data.elev_End_Clip.P2014);
AT_data.array_size.P17_elev = size(AT_data.elev_End_Clip.P2017);
AT_data.array_size.P18_elev = size(AT_data.elev_End_Clip.P2018);
% SURF
AT_data.array_size.P14_elevSURF = size(AT_data.elev_End_Clip_SURF.P2014);
AT_data.array_size.P17_elevSURF = size(AT_data.elev_End_Clip_SURF.P2017);
AT_data.array_size.P18_elevSURF = size(AT_data.elev_End_Clip_SURF.P2018);
% VEL
AT_data.array_size.P14_VEL = size(AT_data.elev_End_Clip_VEL.P2014);
AT_data.array_size.P17_VEL = size(AT_data.elev_End_Clip_VEL.P2017);
AT_data.array_size.P18_VEL = size(AT_data.elev_End_Clip_VEL.P2018);
% LAT
AT_data.array_size.P14_LAT = size(AT_data.elev_End_Clip_LAT.P2014);
AT_data.array_size.P17_LAT = size(AT_data.elev_End_Clip_LAT.P2017);
AT_data.array_size.P18_LAT = size(AT_data.elev_End_Clip_LAT.P2018);
% LON
AT_data.array_size.P14_LON = size(AT_data.elev_End_Clip_LON.P2014);
AT_data.array_size.P17_LON = size(AT_data.elev_End_Clip_LON.P2017);
AT_data.array_size.P18_LON = size(AT_data.elev_End_Clip_LON.P2018);
% DIST
AT_data.array_size.P14_DIST = size(AT_data.elev_End_Clip_DIST.P2014);
AT_data.array_size.P17_DIST = size(AT_data.elev_End_Clip_DIST.P2017);
AT_data.array_size.P18_DIST = size(AT_data.elev_End_Clip_DIST.P2018);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
AT_data.query_array.P14 = (AT_data.AT_vel.P2014(1):10:...
  AT_data.AT_vel.P2014(end));
%AT_data.query_array.P17 = (AT_data.AT_vel.P2017(1):0.1:...
%AT_data.AT_vel.P2017(end));
%AT_data.query_array.P18 = (AT_data.AT_vel.P2018(1):0.1:...
%AT_data.AT_vel.P2018(end));
  
% Apply interpolation to each profile using selected query array
% 2014 bed
AT_data.interp_data.P14 = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip.P2014, AT_data.query_array.P14);
 
% 2017 bed
AT_data.interp_data.P17 = interp1(AT_data.Btrack_End_Clip.P17, ...
  AT_data.elev_End_Clip.P2017, AT_data.query_array.P14);
  
% 2018 bed
AT_data.interp_data.P18 = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip.P2018, AT_data.query_array.P14);

% 2014 surf
AT_data.interp_data.P14_Surf = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_SURF.P2014, AT_data.query_array.P14);
 
% 2017 surf
AT_data.interp_data.P17_Surf = interp1(AT_data.Btrack_End_Clip.P17, ...
  AT_data.elev_End_Clip_SURF.P2017, AT_data.query_array.P14);
  
% 2018 surf
AT_data.interp_data.P18_Surf = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_SURF.P2018, AT_data.query_array.P14);

% 2014 VEL
AT_data.interp_data.P14_VEL = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_VEL.P2014, AT_data.query_array.P14);
 
% 2017 VEL
AT_data.interp_data.P17_VEL = interp1(AT_data.Btrack_End_Clip.P17, ...
  AT_data.elev_End_Clip_VEL.P2017, AT_data.query_array.P14);
  
% 2018 VEL
AT_data.interp_data.P18_VEL = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_VEL.P2018, AT_data.query_array.P14);

% 2014 LAT
AT_data.interp_data.P14_LAT = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LAT.P2014, AT_data.query_array.P14);
 
% 2017 LAT
AT_data.interp_data.P17_LAT = interp1(AT_data.Btrack_End_Clip.P17, ...
  AT_data.elev_End_Clip_LAT.P2017, AT_data.query_array.P14);
  
% 2018 LAT
AT_data.interp_data.P18_LAT = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_LAT.P2018, AT_data.query_array.P14);

% 2014 LON
AT_data.interp_data.P14_LON = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LON.P2014, AT_data.query_array.P14);
 
% 2017 LON
AT_data.interp_data.P17_LON = interp1(AT_data.Btrack_End_Clip.P17, ...
  AT_data.elev_End_Clip_LON.P2017, AT_data.query_array.P14);
  
% 2018 LON
AT_data.interp_data.P18_LON = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_LON.P2018, AT_data.query_array.P14);

% 2014 DIST
AT_data.interp_data.P14_DIST = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_DIST.P2014, AT_data.query_array.P14);
 
% 2017 DIST
AT_data.interp_data.P17_DIST = interp1(AT_data.Btrack_End_Clip.P17, ...
  AT_data.elev_End_Clip_DIST.P2017, AT_data.query_array.P14);
  
% 2018 DIST
AT_data.interp_data.P18_DIST = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip_DIST.P2018, AT_data.query_array.P14);

% Find ice thickness based of difference of interp data
% 2014
AT_data.ice_thickness_interp.P14 = AT_data.interp_data.P14_Surf - AT_data.interp_data.P14;
% 2017
AT_data.ice_thickness_interp.P17 = AT_data.interp_data.P17_Surf - AT_data.interp_data.P17;
% 2018
AT_data.ice_thickness_interp.P18 = AT_data.interp_data.P18_Surf - AT_data.interp_data.P18;

% Calculate melt rates from interpolated profile pairings
% 2014-2017 melt (Vertical difference in Features)
AT_data.melt_rates.P14_P17 = AT_data.interp_data.P17 - ...
  AT_data.interp_data.P14; 
% 2017-2018 melt (Vertical difference in Features)
AT_data.melt_rates.P17_P18 = AT_data.interp_data.P18 - ...
  AT_data.interp_data.P17;
% 2014-2018 melt (Vertical difference in Features)
AT_data.melt_rates.P14_P18 = AT_data.interp_data.P18 - ...
  AT_data.interp_data.P14;
  
% THIS IS HOW TO DO THE INTERPOLATION
%MAKES ALONG TRACK SPACING IN THE PROFILE ALONG TRACK FILE EVERY 0.1M
%x_t14 = (AT_data.AT_vel.P2014(1):0.1:AT_data.AT_vel.P2014(end));
%INTERPOLATION BETWEEN THE SPACED LINE AND NEW YEARS
%FOR 2014
%elev_14_int = interp1(AT_data.Btrack_E.P14, AT_data.elev_EndC.P2014,x_t14);
%FOR 2017
%elev_17_int = interp1(AT_data.Btrack_E.P17, AT_data.elev_EndC.P2017,x_t14);
%FOR 2018
%elev_18_int = interp1(AT_data.Btrack_E.P18, AT_data.elev_EndC.P2018,x_t14);
%% Test Figures Section 

figure(1)
h1 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
h3 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2014', 'original 2017', 'original 2018', ...
  'Location', 'southeast');
  
figure(2)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P14);
hold on
h2 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P17);
h3 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P18);
h4 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
h5 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
h6 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2014', 'interpolated 2017', 'interpolated 2018', ...
  'old 2014','old 2017','old 2018', 'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P14);
hold on
h2 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P17);
h3 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P18);
h4 = plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P18);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2014', 'interpolated 2017', 'interpolated 2018',...
  '2014-2018 melt', 'Location', 'southeast');
  
figure(4)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P17);
hold on
h2 = plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P17_P18);
h3 = plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P18);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2014-2017 melt', '2017-2018 melt', '2014-2018 melt', ...
  'Location', 'southeast');

% Subplot of melt rates for paired years
figure(5)
subplot(3,1,1);
plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P17, 'color','m');
title('Melt 2014-2017');
ylim([-100, 400]);

subplot(3,1,2);
plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P17_P18, 'color','c');
title('Melt 2017-2018');
ylim([-100, 400]);
ylabel('Elevation Change (m)');

subplot(3,1,3);
plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P18, 'color','g');
title('Melt 2014-2018');
ylim([-100, 400]);
xlabel('Along Track Distance (km)');

% melt rate figures
figure(6)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P18/7, 'color','g');
hold on 
h2 = plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P17/3, 'color','c');
h3 = plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P17_P18/4, 'color','m');
title('Melt Average annual melt rates based off vertical differences');
ylim([-20, 100]);
xlabel('Along Track Distance (km)');
ylabel('Vertical Melt (m)');
legend('2014-2018 annual average', '2014-2017 annual average', '2017-2018 annual average', 'Location', 'northeast');

figure(7)
plot(AT_data.query_array.P14/1e3, AT_data.melt_rates.P14_P18/7, 'color','g');
ylim([0, 50]);
title('Average annual difference in heigh between 2014-2018');
xlabel('Along Track Distance (km)');
ylabel('Vertical Melt (m)');
legend('2014-2018 annual average');

%% Crevasse Framework
% Obtain local crevasse apices and low points
Crevasse_apex = islocalmax(AT_data.interp_data.P14);
Crevasse_base = islocalmin(AT_data.interp_data.P14);




figure(2)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P14);
hold on
h2 = plot(AT_data.query_array.P14/1e3, Crevasse_base);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2014', 'crevasse base', 'Location', 'southeast');




%% Basal Crevasse Apex Picking section:
% Crevasse Apex calculation and binary file formation
AT_data.crevasse.apex_bin.P14 = islocalmax(AT_data.interp_data.P14);
AT_data.crevasse.apex_bin.P17 = islocalmax(AT_data.interp_data.P17);
AT_data.crevasse.apex_bin.P18 = islocalmax(AT_data.interp_data.P18);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
AT_data.crevasse.apex_pic_P14 = AT_data.interp_data.P14.*AT_data.crevasse.apex_bin.P14;
AT_data.crevasse.apex_pic_P14(AT_data.crevasse.apex_pic_P14 == 0) = NaN;
AT_data.crevasse.apex_pic_P17 = AT_data.interp_data.P17.*AT_data.crevasse.apex_bin.P17;
AT_data.crevasse.apex_pic_P17(AT_data.crevasse.apex_pic_P17 == 0) = NaN;
AT_data.crevasse.apex_pic_P18 = AT_data.interp_data.P18.*AT_data.crevasse.apex_bin.P18;
AT_data.crevasse.apex_pic_P18(AT_data.crevasse.apex_pic_P18 == 0) = NaN;

%% Basal Crevasse Wall Base Picking:
% Crevasse wall base calculation and binary file formation
AT_data.crevasse.base_bin.P14 = islocalmin(AT_data.interp_data.P14);
AT_data.crevasse.base_bin.P17 = islocalmin(AT_data.interp_data.P17);
AT_data.crevasse.base_bin.P18 = islocalmin(AT_data.interp_data.P18);

% Take elevation values and multiple by local crevasse apex picks
% Turn all 0 values in binary format to NaN for exclusion in plotting
AT_data.crevasse.base_pic_P14 = AT_data.interp_data.P14.*AT_data.crevasse.base_bin.P14;
AT_data.crevasse.base_pic_P14(AT_data.crevasse.base_pic_P14 == 0) = NaN;
AT_data.crevasse.base_pic_P17 = AT_data.interp_data.P17.*AT_data.crevasse.base_bin.P17;
AT_data.crevasse.base_pic_P17(AT_data.crevasse.base_pic_P17 == 0) = NaN;
AT_data.crevasse.base_pic_P18 = AT_data.interp_data.P18.*AT_data.crevasse.base_bin.P18;
AT_data.crevasse.base_pic_P18(AT_data.crevasse.base_pic_P18 == 0) = NaN;
%% Basal Crevasse Plotting
figure(10)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P14,'r');
hold on;
h2 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P17,'b');
h3 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P18,'g');
h4 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.apex_pic_P14,'k*');
h5 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.apex_pic_P17,'k*');
h6 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.apex_pic_P18,'k*');
title('Basal Crevasse Apex Picks - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2014', 'corrected 2017', 'corrected 2018', ...
  'Apex Pick 2014', 'Apex Pick 2017', 'Apex Pick 2018', ...
  'Location', 'southeast');
hold off

figure(14)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P14,'r');
hold on;
h2 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P17,'b');
h3 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P18,'g');
h4 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.base_pic_P14,'k*');
h5 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.base_pic_P17,'k*');
h6 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.base_pic_P18,'k*');
title('Basal Crevasse Base Picks - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2014', 'corrected 2017', 'corrected 2018', ...
  'Base Pick 2014', 'Base Pick 2017', 'Base Pick 2018', ...
  'Location', 'southeast');
hold off;

figure(12)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P14,'r');
hold on
h2 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.base_pic_P14,'k*');
title('Basal Crevasse Base Picks (2014) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2014', 'Base Pick 2014', 'Location', 'southeast');

figure(13)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P17,'b');
hold on
h2 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.base_pic_P17,'k*');
title('Basal Crevasse Base Picks (2017) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2017', 'Base Pick 2017', 'Location', 'southeast');

figure(17)
h1 = plot(AT_data.query_array.P14/1e3, AT_data.interp_data.P18,'g');
hold on
h2 = plot(AT_data.query_array.P14/1e3, AT_data.crevasse.base_pic_P18,'k*');
title('Basal Crevasse Base Picks (2017) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('corrected 2018', 'Base Pick 2018', 'Location', 'southeast');

crevasse1 = AT_data.interp_data.P14(479:497);
along_track_crevasse1 = AT_data.query_array.P14(479:497)/1e3;
cross_area = cumtrapz(along_track_crevasse1, crevasse1);
figure(15)
h1 = plot(along_track_crevasse1, crevasse1);
cross_area1 = patch(along_track_crevasse1, crevasse1,'b','LineWidth',1.5);
title('Basal Crevasse Base Pick 1 (2014) - Profile 1 Petermann Corrected');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('basal crevasse 1', 'Location', 'southeast');

