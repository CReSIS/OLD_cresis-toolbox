  
%% Flightline Data Interpolation - Petermann Line 2
% Years: 2007, 2013, 2014, 2017
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
%AT_data.elevB.P2007 = (pass(1).layers(2).layer_elev);
AT_data.elevB.P2013 = (pass(1).layers(2).layer_elev);
AT_data.elevB.P2014 = (pass(2).layers(2).layer_elev);
AT_data.elevB.P2017 = (pass(3).layers(2).layer_elev);

%Save Surface profiles
%AT_data.elevS.P2007 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2013 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2014 = (pass(2).layers(1).layer_elev);
AT_data.elevS.P2017 = (pass(3).layers(1).layer_elev);

% Save Lidar Surface profiles
AT_data.elev_lidar.P2013 = (pass(1).layers(3).layer_elev);
AT_data.elev_lidar.P2014 = (pass(2).layers(3).layer_elev);
AT_data.elev_lidar.P2017 = (pass(3).layers(3).layer_elev);

%Save annual alongtrack profile data
%AT_data.pass.P2007 = pass(1).along_track;
AT_data.pass.P2013 = pass(1).along_track;
AT_data.pass.P2014 = pass(2).along_track;
AT_data.pass.P2017 = pass(3).along_track;

%Save annual velocity correction data
%AT_data.vel.P2007 = pass(1).vel;
AT_data.vel.P2013 = pass(1).vel;
AT_data.vel.P2014 = pass(2).vel;
AT_data.vel.P2017 = pass(3).vel;

%Save Velocity Corrected Along_track data
%AT_data.AT_vel.P2007 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2013 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2014 = pass(baseline_master_idx).along_track + pass(2).vel;
AT_data.AT_vel.P2017 = pass(baseline_master_idx).along_track + pass(3).vel;

% Latitudes and longitudes (Adjust for the master pass)
AT_data.latitudes.P2013 = pass(1).lat;
AT_data.latitudes.P2014 = interp1(pass(2).lat, pass(2).lat, pass(1).lat, 'nearest','extrap');
AT_data.latitudes.P2017 = interp1(pass(3).lat, pass(3).lat, pass(1).lat, 'nearest','extrap');

AT_data.longitudes.P2013 = pass(1).lon;
AT_data.longitudes.P2014 = interp1(pass(2).lon, pass(2).lon, pass(1).lon, 'nearest','extrap');
AT_data.longitudes.P2017 = interp1(pass(3).lon, pass(3).lon, pass(1).lon, 'nearest','extrap');

%% Section 2 - Mask, indexing, and clipping of along track
  
% Locate Along Track start element in each profile, save element ID as a
% variable for the clipping
% 2007
% AT_data.Btrack.P07 = AT_data.AT_vel.P2007(AT_data.AT_vel.P2007 >= ...
%   AT_data.AT_vel.P2007(1));
% AT_data.find_AT_value.P07 = find(AT_data.AT_vel.P2007 == ...
%   AT_data.Btrack.P07(1));

% 2013  
AT_data.Btrack.P13 = AT_data.AT_vel.P2013(AT_data.AT_vel.P2013 >= ...
  AT_data.AT_vel.P2013(1));
AT_data.find_AT_value.P13 = find(AT_data.AT_vel.P2013 == ...
  AT_data.Btrack.P13(1));

% 2014
AT_data.Btrack.P14 = AT_data.AT_vel.P2014(AT_data.AT_vel.P2014 >= ...
  AT_data.AT_vel.P2013(1));
AT_data.find_AT_value.P14 = find(AT_data.AT_vel.P2014 == ...
  AT_data.Btrack.P14(1));
  
% 2017 
AT_data.Btrack.P17 = AT_data.AT_vel.P2017(AT_data.AT_vel.P2017 >= ...
  AT_data.AT_vel.P2013(1));
AT_data.find_AT_value.P17 = find(AT_data.AT_vel.P2017 == ... 
  AT_data.Btrack.P17(1));

% Clipping from start point in each profile to the end of the profile
% AT_data.Btrack_Beg_Clip.P07 = AT_data.AT_vel.P2007...
%   (AT_data.find_AT_value.P07:end);
AT_data.Btrack_Beg_Clip.P13 = AT_data.AT_vel.P2013...
  (AT_data.find_AT_value.P13:end); 
AT_data.Btrack_Beg_Clip.P14 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.Btrack_Beg_Clip.P17 = AT_data.AT_vel.P2017...
 (AT_data.find_AT_value.P17:end);

% Clipping from new start locations to a given value end element value 
% AT_data.Btrack_End_Clip.P07 = AT_data.Btrack_Beg_Clip.P07...
%   (AT_data.Btrack_Beg_Clip.P07 <= 5.79e+04);
AT_data.Btrack_End_Clip.P13 = AT_data.Btrack_Beg_Clip.P13...
  (AT_data.Btrack_Beg_Clip.P13 <= 5.79e+04);
AT_data.Btrack_End_Clip.P14 = AT_data.Btrack_Beg_Clip.P14...
  (AT_data.Btrack_Beg_Clip.P14 <= 5.79e+04);
AT_data.Btrack_End_Clip.P17 = AT_data.Btrack_Beg_Clip.P17...
 (AT_data.Btrack_Beg_Clip.P17 <= 5.79e+04); 
 
% Save along track data size as variable to see if there is any errors
% AT_data.array_size.P07_AT = size(AT_data.Btrack_End_Clip.P07);
AT_data.array_size.P13_AT = size(AT_data.Btrack_End_Clip.P13);
AT_data.array_size.P14_AT = size(AT_data.Btrack_End_Clip.P14);
AT_data.array_size.P17_AT = size(AT_data.Btrack_End_Clip.P17);

%% Elevation data Clipping to Section size of Along Track files
% Elevation data beginning clipping from start element in Along Track BED
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip.P2013 = AT_data.elevB.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip.P2014 = AT_data.elevB.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip.P2017 = AT_data.elevB.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data BED
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip.P2013 = AT_data.elev_Beg_Clip.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip.P2014 = AT_data.elev_Beg_Clip.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip.P2017 = AT_data.elev_Beg_Clip.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));  

% Elevation data beginning clipping from start element in Along Track SURF
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip_SURF.P2013 = AT_data.elevS.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_SURF.P2014 = AT_data.elevS.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_SURF.P2017 = AT_data.elevS.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data SURF
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip_SURF.P2013 = AT_data.elev_Beg_Clip_SURF.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_SURF.P2014 = AT_data.elev_Beg_Clip_SURF.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_SURF.P2017 = AT_data.elev_Beg_Clip_SURF.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));  


% Elevation data beginning clipping from start element in Along Track LIDAR
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip_LIDAR.P2013 = AT_data.elev_lidar.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_LIDAR.P2014 = AT_data.elev_lidar.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_LIDAR.P2017 = AT_data.elev_lidar.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data LIDAR
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip_LIDAR.P2013 = AT_data.elev_Beg_Clip_LIDAR.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_LIDAR.P2014 = AT_data.elev_Beg_Clip_LIDAR.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_LIDAR.P2017 = AT_data.elev_Beg_Clip_LIDAR.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));  

% Elevation data beginning clipping from start element in Along Track LAT
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip_LAT.P2013 = AT_data.latitudes.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_LAT.P2014 = AT_data.latitudes.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_LAT.P2017 = AT_data.latitudes.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data LAT
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip_LAT.P2013 = AT_data.elev_Beg_Clip_LAT.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_LAT.P2014 = AT_data.elev_Beg_Clip_LAT.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_LAT.P2017 = AT_data.elev_Beg_Clip_LAT.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));  

% Elevation data beginning clipping from start element in Along Track LON
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip_LON.P2013 = AT_data.longitudes.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_LON.P2014 = AT_data.longitudes.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_LON.P2017 = AT_data.longitudes.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data LON
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip_LON.P2013 = AT_data.elev_Beg_Clip_LON.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_LON.P2014 = AT_data.elev_Beg_Clip_LON.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip_LON.P2017 = AT_data.elev_Beg_Clip_LON.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17)); 

% Elevation data beginning clipping from start element in Along Track ATVEL
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip_PASS.P2013 = AT_data.AT_vel.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_PASS.P2014 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip_PASS.P2017 = AT_data.AT_vel.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data ATVEL
% FOR 2014 & 2017 ARRAYS WERE MODIFIED BY ELEMENT SUBTRACTION 
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip_PASS.P2013 = AT_data.elev_Beg_Clip_PASS.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_PASS.P2014 = AT_data.elev_Beg_Clip_PASS.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14)-3);
AT_data.elev_End_Clip_PASS.P2017 = AT_data.elev_Beg_Clip_PASS.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17)-1); 

% Save along elevation data size as variable to see if there is any errors
% BED
AT_data.array_size.P13_elev = size(AT_data.elev_End_Clip.P2013);
AT_data.array_size.P14_elev = size(AT_data.elev_End_Clip.P2014);
AT_data.array_size.P17_elev = size(AT_data.elev_End_Clip.P2017);
% SURF
AT_data.array_size.P13_elevSURF = size(AT_data.elev_End_Clip_SURF.P2013);
AT_data.array_size.P14_elevSURF = size(AT_data.elev_End_Clip_SURF.P2014);
AT_data.array_size.P17_elevSURF = size(AT_data.elev_End_Clip_SURF.P2017);
% LAT
AT_data.array_size.P13_LAT = size(AT_data.elev_End_Clip_LAT.P2013);
AT_data.array_size.P14_LAT = size(AT_data.elev_End_Clip_LAT.P2014);
AT_data.array_size.P17_LAT = size(AT_data.elev_End_Clip_LAT.P2017);
% LON
AT_data.array_size.P13_LON = size(AT_data.elev_End_Clip_LON.P2013);
AT_data.array_size.P14_LON = size(AT_data.elev_End_Clip_LON.P2014);
AT_data.array_size.P17_LON = size(AT_data.elev_End_Clip_LON.P2017);
% ATVEL
AT_data.array_size.P13_PASS = size(AT_data.elev_End_Clip_PASS.P2013);
AT_data.array_size.P14_PASS = size(AT_data.elev_End_Clip_PASS.P2014);
AT_data.array_size.P17_PASS = size(AT_data.elev_End_Clip_PASS.P2017);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
% AT_data.query_array.P07 = (AT_data.AT_vel.P2007(1):10:...
%   AT_data.AT_vel.P2007(end));
AT_data.query_array.P13 = (AT_data.AT_vel.P2013(1):15:...
AT_data.AT_vel.P2013(end));
%AT_data.query_array.P14 = (AT_data.AT_vel.P2014(1):10:...
%AT_data.AT_vel.P2014(end));
%AT_data.query_array.P17 = (AT_data.AT_vel.P2017(1):10:...
%AT_data.AT_vel.P2017(end));

% Apply interpolation to each profile using selected query array
% 2007
% AT_data.interp_data.P07 = interp1(AT_data.Btrack_End_Clip.P07, ...
%   AT_data.elev_End_Clip.P2007, AT_data.query_array.P07);

% 2013 BED
AT_data.interp_data.P13 = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip.P2013, AT_data.query_array.P13,'nearest','extrap');
% 2014 BED
AT_data.interp_data.P14 = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip.P2014, AT_data.query_array.P13,'nearest','extrap');
% 2017 BED
AT_data.interp_data.P17 = interp1(AT_data.Btrack_End_Clip.P17, ...
 AT_data.elev_End_Clip.P2017, AT_data.query_array.P13,'nearest','extrap');

% 2013 SURF
AT_data.interp_data.P13_SURF = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_SURF.P2013, AT_data.query_array.P13,'nearest','extrap');
% 2014 SURF
AT_data.interp_data.P14_SURF = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_SURF.P2014, AT_data.query_array.P13,'nearest','extrap');
% 2017 SURF
AT_data.interp_data.P17_SURF = interp1(AT_data.Btrack_End_Clip.P17, ...
 AT_data.elev_End_Clip_SURF.P2017, AT_data.query_array.P13,'nearest','extrap');

% 2013 THICKNESS
AT_data.interp_data.P13_thickness = AT_data.interp_data.P13_SURF - AT_data.interp_data.P13;
% 2014 THICKNESS
AT_data.interp_data.P14_thickness = AT_data.interp_data.P14_SURF - AT_data.interp_data.P14;
% 2017 THICKNESS
AT_data.interp_data.P17_thickness = AT_data.interp_data.P17_SURF - AT_data.interp_data.P17;

% 2013 LIDAR
AT_data.interp_data.P13_LIDAR = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_LIDAR.P2013, AT_data.query_array.P13,'nearest','extrap');
% 2014 LIDAR
AT_data.interp_data.P14_LIDAR = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LIDAR.P2014, AT_data.query_array.P13,'nearest','extrap');
% 2017 LIDAR
AT_data.interp_data.P17_LIDAR = interp1(AT_data.Btrack_End_Clip.P17, ...
 AT_data.elev_End_Clip_LIDAR.P2017, AT_data.query_array.P13,'nearest','extrap');

% 2013 LAT
AT_data.interp_data.P13_LAT = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_LAT.P2013, AT_data.query_array.P13,'nearest','extrap');
% 2014 LAT
AT_data.interp_data.P14_LAT = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LAT.P2014, AT_data.query_array.P13,'nearest','extrap');
% 2017 LAT
AT_data.interp_data.P17_LAT = interp1(AT_data.Btrack_End_Clip.P17, ...
 AT_data.elev_End_Clip_LAT.P2017, AT_data.query_array.P13,'nearest','extrap');

% 2013 LON
AT_data.interp_data.P13_LON = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_LON.P2013, AT_data.query_array.P13,'nearest','extrap');
% 2014 LON
AT_data.interp_data.P14_LON = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip_LON.P2014, AT_data.query_array.P13,'nearest','extrap');
% 2017 LON
AT_data.interp_data.P17_LON = interp1(AT_data.Btrack_End_Clip.P17, ...
 AT_data.elev_End_Clip_LON.P2017, AT_data.query_array.P13,'nearest','extrap');

% 2013 ATVEL
AT_data.interp_data.P13_PASS = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_PASS.P2013, AT_data.query_array.P13,'nearest','extrap');
% 2014 ATVEL
AT_data.interp_data.P14_PASS = interp1(AT_data.Btrack_End_Clip.P14(1:end-3), ...
  AT_data.elev_End_Clip_PASS.P2014, AT_data.query_array.P13,'nearest','extrap');
% 2017 ATVEL
AT_data.interp_data.P17_PASS = interp1(AT_data.Btrack_End_Clip.P17(1:end-1), ...
 AT_data.elev_End_Clip_PASS.P2017, AT_data.query_array.P13,'nearest','extrap');


% Calculate melt rates from interpolated profile pairings
% 2007-2013 melt (Vertical difference in Features)
% AT_data.melt_rates.P07_P13 = AT_data.interp_data.P13 - ...
%   AT_data.interp_data.P07; 
  
% 2013-2014 melt (Vertical difference in Features)
AT_data.melt_rates.P13_P14 = AT_data.interp_data.P14 - ...
  AT_data.interp_data.P13;
  
% 2014-2017 melt (Vertical difference in Features)
AT_data.melt_rates.P14_P17 = AT_data.interp_data.P17 - ...
  AT_data.interp_data.P14;

% 2013-2017 melt (Vertical difference in Features)
AT_data.melt_rates.P13_P17 = AT_data.interp_data.P17 - ...
  AT_data.interp_data.P13;
  
%% Export data to csv
% Concatenate and take transpose of Lon, Lat, Surf, Bed fields. Concatenate horizontally for each year
AT_data.export.P13 = cat(2, AT_data.interp_data.P13_LON.', AT_data.interp_data.P13_LAT.', AT_data.interp_data.P13_thickness.',...
    AT_data.interp_data.P13_SURF.', AT_data.interp_data.P13_LIDAR.', AT_data.interp_data.P13.', AT_data.interp_data.P13_PASS.' );
AT_data.export.P14 = cat(2, AT_data.interp_data.P14_LON.', AT_data.interp_data.P14_LAT.', AT_data.interp_data.P14_thickness.', ...
    AT_data.interp_data.P14_SURF.', AT_data.interp_data.P14_LIDAR.', AT_data.interp_data.P14.', AT_data.interp_data.P14_PASS.');
AT_data.export.P17 = cat(2, AT_data.interp_data.P17_LON.', AT_data.interp_data.P17_LAT.', AT_data.interp_data.P17_thickness.', ...
    AT_data.interp_data.P17_SURF.', AT_data.interp_data.P17_LIDAR.', AT_data.interp_data.P17.', AT_data.interp_data.P17_PASS.' );

%% Define Header Array of strings and vertically concatenate to data 
cheader = {'Lons', 'Lats', 'Thickness','Surface', 'Lidar', 'Depth', 'AT_Velocity'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

% change folder
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\CSV_export_files\'

%write header to file 2013
fid = fopen('P2_lat_lon_surf_bed_13.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P2_lat_lon_surf_bed_13.csv', AT_data.export.P13, '-append');

%write header to file 2014
fid = fopen('P2_lat_lon_surf_bed_14.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P2_lat_lon_surf_bed_14.csv', AT_data.export.P14, '-append');

%write header to file 2017
fid = fopen('P2_lat_lon_surf_bed_17.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P2_lat_lon_surf_bed_17.csv', AT_data.export.P17, '-append');

%% Smoothing and Regression line
% Clip 2013 Profile terminus off (velocity correction too short)
AT_data.Btrack_End_Clip.P13_Clip = AT_data.Btrack_End_Clip.P13(1:3348);
AT_data.elev_End_Clip.P2013_Clip = AT_data.elev_End_Clip.P2013(1:...
    length(AT_data.Btrack_End_Clip.P13_Clip));
AT_data.elev_End_Clip_SURF.P2013_Clip = AT_data.elev_End_Clip_SURF.P2013(1:...
    length(AT_data.Btrack_End_Clip.P13_Clip));
% smooth data using gaussian weighted moving average and window size of 10
% First smoothing uses gaussian weight, quadratic linear regression (loess)
% and robust linear regression (rloess)
% 2013
AT_data.smoothed.P13 = smoothdata(AT_data.elev_End_Clip.P2013_Clip,'gaussian',100);
AT_data.smoothed.P13_loess = smoothdata(AT_data.elev_End_Clip.P2013_Clip,'loess',100);
AT_data.smoothed.P13_rloess = smoothdata(AT_data.elev_End_Clip.P2013_Clip,'rloess',100);
% 2014
AT_data.smoothed.P14 = smoothdata(AT_data.elev_End_Clip.P2014,'gaussian',100);
AT_data.smoothed.P14_loess = smoothdata(AT_data.elev_End_Clip.P2014,'loess',100);
AT_data.smoothed.P14_rloess = smoothdata(AT_data.elev_End_Clip.P2014,'rloess',100);
% 2017
AT_data.smoothed.P17 = smoothdata(AT_data.elev_End_Clip.P2017,'gaussian',100);
AT_data.smoothed.P17_loess = smoothdata(AT_data.elev_End_Clip.P2017,'loess',100);
AT_data.smoothed.P17_rloess = smoothdata(AT_data.elev_End_Clip.P2017,'rloess',100);

% Regression using Polyfit
% uses query array, along track clip array (Btrack_end_clip per year or 
% regression_query, the respective data points (smoothed data) and a 
% scalar degree value (6).
% 2013
AT_data.regression_query.P13 = linspace(AT_data.Btrack_End_Clip.P13_Clip(1), ...
    AT_data.Btrack_End_Clip.P13_Clip(end), length(AT_data.Btrack_End_Clip.P13_Clip));
AT_data.regression_points.P13 = polyfit(AT_data.Btrack_End_Clip.P13_Clip, ...
    AT_data.smoothed.P13, 6);
AT_data.regression_array.P13 = polyval(AT_data.regression_points.P13, ...
    AT_data.Btrack_End_Clip.P13_Clip);
% 2014
AT_data.regression_query.P14 = linspace(AT_data.Btrack_End_Clip.P14(1), ...
    AT_data.Btrack_End_Clip.P14(end), length(AT_data.Btrack_End_Clip.P14));
AT_data.regression_points.P14 = polyfit(AT_data.Btrack_End_Clip.P14, ...
    AT_data.smoothed.P14, 6);
AT_data.regression_array.P14 = polyval(AT_data.regression_points.P14, ...
    AT_data.Btrack_End_Clip.P14);
% 2017
AT_data.regression_query.P17 = linspace(AT_data.Btrack_End_Clip.P17(1), ...
    AT_data.Btrack_End_Clip.P17(end), length(AT_data.Btrack_End_Clip.P17));
AT_data.regression_points.P17 = polyfit(AT_data.Btrack_End_Clip.P17, ...
    AT_data.smoothed.P17, 6);
AT_data.regression_array.P17 = polyval(AT_data.regression_points.P17, ...
    AT_data.Btrack_End_Clip.P17);

% Hdyrostatic Equilibrium Profiles
% From Publications: 
% https://www.cambridge.org/core/journals/journal-of-glaciology/article/variability-of-basal-melt-beneath-the-pine-island-glacier-ice-shelf-west-antarctica/F26EDA3C49D3A5140F702C21434FBAA7
% https://www.cambridge.org/core/journals/journal-of-glaciology/article/interannual-changes-of-the-floating-ice-shelf-of-petermann-gletscher-north-greenland-from-2000-to-2012/5DA6EADBD7E2C73FC0B79505689D4D3D#R4
% Cameron Lewis Scott Dissertation
% Define constants
rho_ice = 917;
rho_water = 1026;
Firn_air_coeff = 0;
% Bindschadler equation for Hydrostatic Equilibtrium
% Derive Thickness H
% Z_s = (1 - (rho_i/rho_w))*H +(rho_i/rho_w)*firn_air_coeff

% 2013
AT_data.thickness.P13_Clip = AT_data.elev_End_Clip_SURF.P2013_Clip; %- AT_data.elev_End_Clip.P2013_Clip;
AT_data.hydrostatic.P13 = -(((AT_data.thickness.P13_Clip - Firn_air_coeff)*rho_water)/(rho_water - rho_ice)) +Firn_air_coeff;
% 2014
AT_data.thickness.P14 = AT_data.elev_End_Clip_SURF.P2014; %- AT_data.elev_End_Clip.P2014;
AT_data.hydrostatic.P14 = -(((AT_data.thickness.P14 - Firn_air_coeff)*rho_water)/(rho_water - rho_ice)) +Firn_air_coeff;
% 2017
AT_data.thickness.P17 = AT_data.elev_End_Clip_SURF.P2017; %- AT_data.elev_End_Clip.P2017;
AT_data.hydrostatic.P17 = - (((AT_data.thickness.P17 - Firn_air_coeff)*rho_water)/(rho_water - rho_ice)) +Firn_air_coeff;

%(((rho_ice/rho_water)*AT_data.thickness.P17) + ((rho_ice/rho_water)*Firn_air_coeff));

% Get number of fields in structure then automate
fn = fieldnames(AT_data.thickness);
fn1 = fieldnames(AT_data.elev_End_Clip_SURF);
for i = 1:numel(fn1)
    AT_data.thickness1.(fn1{i}) = AT_data.elev_End_Clip_SURF.(fn1{i}) - AT_data.elev_End_Clip.(fn1{i});
end
    
% Surface Derived Bed Profile
% Use hydrostatic equilibrium assumption that 1/10th of the shelf thickness
% is above the waterline, while 9/10 of the thickness is below the water
% multiply each surface profile by -9/10 to produce negative elvation value
AT_data.surf_derived.P13 = (AT_data.elev_End_Clip_SURF.P2013_Clip - AT_data.elev_End_Clip.P2013_Clip)*(-9/10);
AT_data.surf_derived.P14 = (AT_data.elev_End_Clip_SURF.P2014 - AT_data.elev_End_Clip.P2014)*(-9/10);
AT_data.surf_derived.P17 = (AT_data.elev_End_Clip_SURF.P2017 - AT_data.elev_End_Clip.P2017)*(-9/10);

%% 1 through 5 order regressions for 2013
% 5-order
AT_data.regression_points.P13_5 = polyfit(AT_data.Btrack_End_Clip.P13_Clip, ...
    AT_data.smoothed.P13, 5);
AT_data.regression_array.P13_5 = polyval(AT_data.regression_points.P13_5, ...
    AT_data.Btrack_End_Clip.P13_Clip);
% 4-order
AT_data.regression_points.P13_4 = polyfit(AT_data.Btrack_End_Clip.P13_Clip, ...
    AT_data.smoothed.P13, 4);
AT_data.regression_array.P13_4 = polyval(AT_data.regression_points.P13_4, ...
    AT_data.Btrack_End_Clip.P13_Clip);
% 3-order
AT_data.regression_points.P13_3 = polyfit(AT_data.Btrack_End_Clip.P13_Clip, ...
    AT_data.smoothed.P13, 3);
AT_data.regression_array.P13_3 = polyval(AT_data.regression_points.P13_3, ...
    AT_data.Btrack_End_Clip.P13_Clip);
% 2-order
AT_data.regression_points.P13_2 = polyfit(AT_data.Btrack_End_Clip.P13_Clip, ...
    AT_data.smoothed.P13, 2);
AT_data.regression_array.P13_2 = polyval(AT_data.regression_points.P13_2, ...
    AT_data.Btrack_End_Clip.P13_Clip);
% 1-order
AT_data.regression_points.P13_1 = polyfit(AT_data.Btrack_End_Clip.P13_Clip, ...
    AT_data.smoothed.P13, 1);
AT_data.regression_array.P13_1 = polyval(AT_data.regression_points.P13_1, ...
    AT_data.Btrack_End_Clip.P13_Clip);
% 2014 1-order
AT_data.regression_points.P14_1 = polyfit(AT_data.Btrack_End_Clip.P14, ...
    AT_data.smoothed.P14, 1);
AT_data.regression_array.P14_1 = polyval(AT_data.regression_points.P14_1, ...
    AT_data.Btrack_End_Clip.P14);
% 2017 1-order
AT_data.regression_points.P17_1 = polyfit(AT_data.Btrack_End_Clip.P17, ...
    AT_data.smoothed.P17, 1);
AT_data.regression_array.P17_1 = polyval(AT_data.regression_points.P17_1, ...
    AT_data.Btrack_End_Clip.P17);

% TEST FIGURE FOR Multi-Order Regression
figure(105)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13_1);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13_2);
h3 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13_3);
h4 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13_4);
h5 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13_5);
h6 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13);
title('P13 1-6 Order regression');
xlabel('along track distance (km)');
ylabel('elevation (m)');
legend('1-order', '2-order', '3-order', '4-order', '5-order', '6-order',...
    'Location', 'southeast');

% TESt FIGURE 1-Order Regression all years
figure(106)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13_1, 'color', 'r');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.regression_array.P14_1, 'color', 'b');
h3 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.regression_array.P17_1, 'color', 'g');
title('1-Order Regression 2013-2017');
xlabel('along track distance (km)');
ylabel('elevation (m)');
legend('2013', '2014', '2017', 'Location', 'southeast');

% TESt FIGURE Surface Derived Depths for all years
figure(109)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.surf_derived.P13, 'color', 'r');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.surf_derived.P14, 'color', 'b');
h3 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.surf_derived.P17, 'color', 'g');
h4 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.elev_End_Clip_SURF.P2013_Clip, 'color', 'r');
h5 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip_SURF.P2014, 'color', 'b');
h6 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip_SURF.P2017, 'color', 'g');
title('Surface Derived Basal Profiles 2013-2017');
xlabel('along track distance (km)');
ylabel('elevation (m)');
legend('2013', '2014', '2017', 'Location', 'southeast');

%% TEST FIGURE 3 x 3 original profile, Bindschadler-derived, and Regression
f1 = figure(110);
subplot(3,3,1)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.elev_End_Clip.P2013_Clip, 'color','r');
title('Original 2013');
subplot(3,3,2)
h2 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.hydrostatic.P13, 'color', 'r');
title('Hydrostatic 2013');
subplot(3,3,3)
h3 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13, '--', 'color','r');
title('6-Order Regression 2013');
subplot(3,3,4)
h4 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014, 'color','b');
title('Original 2014');
subplot(3,3,5)
h5 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.hydrostatic.P14, 'color', 'b');
title('Hydrostatic 2014');
subplot(3,3,6)
h6 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.regression_array.P14, '--', 'color','b');
title('6-Order Regression 2014');
subplot(3,3,7)
h7 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017, 'color','g');
title('Original 2017');
subplot(3,3,8)
h8 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.hydrostatic.P17, 'color', 'g');
title('Hydrostatic 2017');
subplot(3,3,9)
h9 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.regression_array.P17, '--', 'color','g');
title('6-Order Regression 2017');

f1_ax = findobj(f1,'Type','Axes');
for i=1:length(f1_ax)
    xlabel(f1_ax(i),{'AT Distance (km)'});
    ylabel(f1_ax(i),{'Elevation (m)'});
    xlim([0 60]);
    ylim([-500 -100]);
end
% original_array = [9 6 3];
% hydro_array = [2 5 8];
% regress_array = [1 4 7];
% num_array = [3 4 7];
% for i= original_array
%     for j = num_array
%     title_Label = sprintf('Original 201%d',j);
%     title(f1_ax(i),{title_Label}, 'FontSize',10);
%     end
% end
%%
f2 = figure(111);
subplot(3,1,1);
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.elev_End_Clip.P2013_Clip, 'color','r');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.hydrostatic.P13, 'color', 'g');
h3 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13, '--', 'color','b');
%h4 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.surf_derived.P13, '--', 'color', 'k');
title('2013 Original, Hydostatic, & Regression');
subplot(3,1,2);
h5 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014, 'color','r');
hold on
h6 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.hydrostatic.P14, 'color', 'g');
h7 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.regression_array.P14, '--', 'color','b');
%h8 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.surf_derived.P14, '--', 'color', 'k');
title('2014 Original, Hydostatic, & Regression');
subplot(3,1,3);
h9 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017, 'color','r');
hold on
h10 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.hydrostatic.P17, 'color', 'g');
h11 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.regression_array.P17, '--', 'color','b');
%h12 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.surf_derived.P17, '--', 'color', 'k');
title('2017 Original, Hydostatic, & Regression');

f2_ax = findobj(f2,'Type','Axes');
num_array = [3 4 7];
for j = num_array
    for i=1:length(f2_ax)
        xlabel(f2_ax(i),{'AT Distance (km)'});
        ylabel(f2_ax(i),{'Elevation (m)'});
        xlim([0 60]);
        ylim([-500 0]);
        legend(f2_ax(i),{'Original', 'Hydrostatic','Regression','9/10th method'}, 'Location','southeast');
        %title_Label = sprintf('Original 201%d',j);
        %title(f2_ax(i),{title_Label}, 'FontSize',10);
    end
end

%% Find melt between 10 & 15 km
% 10km Point location AT
AT_data.point_10km.P13 = AT_data.Btrack_End_Clip.P13_Clip(AT_data.Btrack_End_Clip.P13_Clip >= 10000);
AT_data.point_10km.P14 = AT_data.Btrack_End_Clip.P14(AT_data.Btrack_End_Clip.P14 >= 10000);
AT_data.point_10km.P17 = AT_data.Btrack_End_Clip.P17(AT_data.Btrack_End_Clip.P17 >= 10000);


AT_data.find_AT_value.P17 = find(AT_data.AT_vel.P2017 == ...
  AT_data.Btrack.P17(1));

%%
% Obtain mean value per year and make into horizontal line same
% 2013
AT_data.mean_elev_bed.P13 = mean(AT_data.smoothed.P13);
AT_data.mean_elev_array.P13 = AT_data.mean_elev_bed.P13.* ones(length(AT_data.Btrack_End_Clip.P13),1);
% 2014
AT_data.mean_elev_bed.P14 = mean(AT_data.smoothed.P14);
AT_data.mean_elev_array.P14 = AT_data.mean_elev_bed.P14.* ones(length(AT_data.Btrack_End_Clip.P14),1);
% 2017
AT_data.mean_elev_bed.P17 = mean(AT_data.smoothed.P17);
AT_data.mean_elev_array.P17 = AT_data.mean_elev_bed.P17.* ones(length(AT_data.Btrack_End_Clip.P17),1);

% TEST FIGURE FOR SMOOTHING
figure(100)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.elev_End_Clip.P2013_Clip);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.smoothed.P13);
h3 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.smoothed.P13_loess);
h4 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.smoothed.P13_rloess);
title('P13 original & Smoothed');
xlabel('along track distance (km)');
ylabel('elevation (m)');
legend('original', 'guassian', 'loess', 'rloess', 'Location', 'southeast');

% TEST FIGURE FOR REGRESSION
figure(101)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.elev_End_Clip.P2013_Clip);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.smoothed.P13);
h3 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.smoothed.P13_loess);
h4 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.smoothed.P13_rloess);
h5 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13);
title('P13 original & Smoothed');
xlabel('along track distance (km)');
ylabel('elevation (m)');
legend('original', 'gaussian', 'loess', 'rloess', 'regression', 'Location', 'southeast');

%% figure(108)
figure(108)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13, '--r');
h4 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.elev_End_Clip.P2013_Clip, 'color','r');
title('6th Order Ice Shelf Regression Lines 2013-2017');
xlabel('Along Track Distance (km)');
ylabel('Elevation (m)');
ylim([-500 -100]);
%% Regression Figure with 6th order fit
figure(102)
subplot(3,1,1)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13, '--', 'color','r');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.regression_array.P14, '--', 'color', 'b');
h3 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.regression_array.P17, '--', 'color', 'g');
x1 = xline(5.2351,':r',{'2013 GL'});
x1.LabelHorizontalAlignment = 'left';
x2 = xline(5.0549,':b',{'2014 GL'});
x2.LabelVerticalAlignment = 'middle';
x2.LabelHorizontalAlignment = 'left';
x3 = xline(5.2499,':g',{'2017 GL'});
xline(10);
xline(20);
xline(30);
xline(40);
xline(50);
title('6th Order Ice Shelf Regression Lines 2013-2017');
xlabel('Along Track Distance (km)');
ylabel('Elevation (m)');
ylim([-500 -100]);
legend('regression 2013', 'regression 2014', 'regression 2017', 'Location', 'southeast');
subplot(3,1,2)
h4 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.elev_End_Clip.P2013_Clip, 'color','r');
hold on
h5 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014, 'color','b');
h6 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017, 'color','g');
title('Original Ice Shelf Profiles 2013-2017');
xlabel('Along Track Distance (km)');
ylabel('Elevation (m)');
ylim([-500 -100]);
legend('original 2013', 'original 2014', 'original 2017', 'Location', 'southeast');
subplot(3,1,3)
h1 = plot(AT_data.Btrack_End_Clip.P13_Clip/1e3, AT_data.regression_array.P13, '--', 'color','r');
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.regression_array.P14, '--', 'color', 'b');
h3 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.regression_array.P17, '--', 'color', 'g');
h7 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.mean_elev_array.P13, 'color','r');
h8 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.mean_elev_array.P14, 'color','b');
h9 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.mean_elev_array.P17, 'color','g');
title('Regression lines and Mean Profle Elevations 2013-2017');
xlabel('Along Track Distance (km)');
ylabel('Elevation (m)');
ylim([-500 -100]);
legend('regression 2013', 'regression 2014', 'regression 2017', ...
    'mean 2013', 'mean 2014', 'mean 2017', 'Location', 'southeast', 'NumColumns', 2);

% subplot(4,1,3)
% h7 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.mean_elev_array.P13, 'color','r');
% hold on 
% h8 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.mean_elev_array.P14, 'color','b');
% h9 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.mean_elev_array.P17, 'color','g');
% title('Mean Profle Elevations 2013-17');
% xlabel('Along Track Distance (km)');
% ylabel('Elevation (m)');
% ylim([-500 -100]);
% legend('mean 2013', 'mean 2014', 'mean 2017', 'Location', 'southeast');

%% define along track segments 
% 0 - 5 km
% 54.9142
AT_segment_1_P13 = AT_data.Btrack_End_Clip.P13(AT_data.Btrack_End_Clip.P13 < 10000);
AT_segment_1_P14 = AT_data.Btrack_End_Clip.P14(AT_data.Btrack_End_Clip.P14 < 10000);
AT_segment_1_P17 = AT_data.Btrack_End_Clip.P17(AT_data.Btrack_End_Clip.P17 < 10000);

first_segment_13 = AT_data.smoothed.P13(1:length(AT_segment_1_P13));
first_segment_14 = AT_data.smoothed.P14(1:length(AT_segment_1_P14));
first_segment_17 = AT_data.smoothed.P17(1:length(AT_segment_1_P17));

figure(101)
h1 = plot(AT_segment_1_P13/1e3, first_segment_13);
hold on 
h2 = plot(AT_segment_1_P14/1e3, first_segment_14);
h3 = plot(AT_segment_1_P17/1e3, first_segment_17);
xlabel('along track distance (km)');
ylabel('elevation (m)');
legend('original', 'smoothed', 'Location', 'southeast');

% 2013  
% AT_data.Btrack.P13 = AT_data.AT_vel.P2013(AT_data.AT_vel.P2013 >= ...
%   AT_data.AT_vel.P2013(1));
% AT_data.find_AT_value.P13 = find(AT_data.AT_vel.P2013 == ...
%   AT_data.Btrack.P13(1));

% AT_data.elev_Beg_Clip.P2013 = AT_data.elevB.P2013...
%   (AT_data.find_AT_value.P13:end);

%%
x = linspace(0,4*pi,10);
y = sin(x);
p = polyfit(x,y,7);
x1 = linspace(0,4*pi);
y1 = polyval(p,x1);

%% Test Figures Section 

figure(1)
%h1 = plot(AT_data.Btrack_End_Clip.P07/1e3, AT_data.elev_End_Clip.P2007);
h2 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.elev_End_Clip.P2013);
hold on
h3 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
h4 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2013', 'original 2014', ...
  'original 2017', 'Location', 'southeast');
%%  
figure(2)
h1 = plot(AT_data.query_array.P07, AT_data.interp_data.P07);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.interp_data.P13);
h3 = plot(AT_data.query_array.P07, AT_data.interp_data.P14);
h4 = plot(AT_data.query_array.P07, AT_data.interp_data.P17);
h5 = plot(AT_data.Btrack_End_Clip.P07, AT_data.elev_End_Clip.P2007);
h6 = plot(AT_data.Btrack_End_Clip.P13, AT_data.elev_End_Clip.P2013);
h7 = plot(AT_data.Btrack_End_Clip.P14, AT_data.elev_End_Clip.P2014);
h8 = plot(AT_data.Btrack_End_Clip.P17, AT_data.elev_End_Clip.P2017);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2007', 'interpolated 2013', 'interpolated 2014', ...
  'interpolated 2017', 'old 2007','old 2013','old 2014', 'old 2017', ...
  'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P07, AT_data.interp_data.P07);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.interp_data.P13);
h3 = plot(AT_data.query_array.P07, AT_data.interp_data.P14);
h4 = plot(AT_data.query_array.P07, AT_data.interp_data.P17);
h5 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P17);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2007', 'interpolated 2013', 'interpolated 2014',...
  'interpolated 2017', '2007-2017 melt', 'Location', 'southeast');
  
figure(4)
h1 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P13);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.melt_rates.P13_P14);
h3 = plot(AT_data.query_array.P07, AT_data.melt_rates.P14_P17);
h4 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P17);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2007-2013 melt', '2013-2014 melt', '2014-2017 melt', ...
  '2007-2017 melt', 'Location', 'southeast');

% Subplot figure of melt rates
figure(5)
subplot(4,1,1);
plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P13);
title('Melt 2007-2013');

subplot(4,1,2);
plot(AT_data.query_array.P07, AT_data.melt_rates.P13_P14);
title('Melt 2013-2014');

subplot(4,1,3);
plot(AT_data.query_array.P07, AT_data.melt_rates.P14_P17);
title('Melt 2014-2017');

subplot(4,1,4);
plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P17);
title('Melt 2007-2014');

% Annual melt rate average based off 2007-2019 melt subtraction
% figure(6)
% plot(AT_data.query_array.P07/1e3, AT_data.melt_rates.P07_P17/10, 'color','g');
% title('Melt 2007-2017');
% ylim([-070, 400]);
% xlabel('Along Track Distance (km)');
