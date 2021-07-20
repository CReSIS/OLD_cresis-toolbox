  
%% Flightline Data Interpolation - Ryder Line 1
% Years: 2007, 2011, 2013, 2015, 2019
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
AT_data.elevB.P2011 = (pass(1).layers(2).layer_elev);
AT_data.elevB.P2013 = (pass(2).layers(2).layer_elev);
AT_data.elevB.P2015 = (pass(3).layers(2).layer_elev);
%AT_data.elevB.P2019 = (pass(4).layers(2).layer_elev);

%Save Surface profiles
%AT_data.elevS.P2007 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2011 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2013 = (pass(2).layers(1).layer_elev);
AT_data.elevS.P2015 = (pass(3).layers(1).layer_elev);
%AT_data.elevS.P2019 = (pass(4).layers(1).layer_elev);

%Save annual alongtrack profile data
%AT_data.pass.P2007 = pass(1).along_track;
AT_data.pass.P2011 = pass(1).along_track;
AT_data.pass.P2013 = pass(2).along_track;
AT_data.pass.P2015 = pass(3).along_track;
%AT_data.pass.P2019 = pass(4).along_track;

%Save annual velocity correction data
%AT_data.vel.P2007 = pass(1).vel;
AT_data.vel.P2011 = pass(1).vel;
AT_data.vel.P2013 = pass(2).vel;
AT_data.vel.P2015 = pass(3).vel;
%AT_data.vel.P2019 = pass(4).vel;

%Save Velocity Corrected Along_track data
%AT_data.AT_vel.P2007 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2011 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2013 = pass(baseline_master_idx).along_track + pass(2).vel;
AT_data.AT_vel.P2015 = pass(baseline_master_idx).along_track + pass(3).vel;
%AT_data.AT_vel.P2019 = pass(baseline_master_idx).along_track + pass(4).vel;

% Latitudes and longitudes (Adjust for the master pass)  
AT_data.latitudes.P2011 = interp1(pass(1).lat, pass(1).lat, pass(4).lat);
AT_data.latitudes.P2013 = interp1(pass(2).lat, pass(2).lat, pass(4).lat);
AT_data.latitudes.P2015 = interp1(pass(3).lat, pass(3).lat, pass(4).lat);
%AT_data.latitudes.P2019 = pass(4).lat;
AT_data.longitudes.P2011 = interp1(pass(1).lon, pass(1).lon, pass(4).lon);
AT_data.longitudes.P2013 = interp1(pass(2).lon, pass(2).lon, pass(4).lon);
AT_data.longitudes.P2015 = interp1(pass(3).lon, pass(3).lon, pass(4).lon);
%AT_data.longitudes.P2019 = pass(4).lon;
%% Section 2 - Mask, indexing, and clipping of along track
  
% Locate Along Track start element in each profile, save element ID as a
% variable for the clipping
% 2007
% AT_data.Btrack.P07 = AT_data.AT_vel.P2007(AT_data.AT_vel.P2007 >= ...
%   AT_data.AT_vel.P2007(1));
% AT_data.find_AT_value.P07 = find(AT_data.AT_vel.P2007 == ...
%   AT_data.Btrack.P07(1));

% 2011  
AT_data.Btrack.P11 = AT_data.AT_vel.P2011(AT_data.AT_vel.P2011 >= ...
  AT_data.AT_vel.P2011(1));
AT_data.find_AT_value.P11 = find(AT_data.AT_vel.P2011 == ...
  AT_data.Btrack.P11(1));

% 2013
AT_data.Btrack.P13 = AT_data.AT_vel.P2013(AT_data.AT_vel.P2013 >= ...
  AT_data.AT_vel.P2011(1));
AT_data.find_AT_value.P13 = find(AT_data.AT_vel.P2013 == ...
  AT_data.Btrack.P13(1));
  
% 2015 
AT_data.Btrack.P15 = AT_data.AT_vel.P2015(AT_data.AT_vel.P2015 >= ...
  AT_data.AT_vel.P2011(1));
AT_data.find_AT_value.P15 = find(AT_data.AT_vel.P2015 == ... 
  AT_data.Btrack.P15(1));

% 2019 
% AT_data.Btrack.P19 = AT_data.AT_vel.P2019(AT_data.AT_vel.P2019 >= ...
%   AT_data.AT_vel.P2011(1));
% AT_data.find_AT_value.P19 = find(AT_data.AT_vel.P2019 == ... 
%   AT_data.Btrack.P19(1));

% Clipping from start point in each profile to the end of the profile
% AT_data.Btrack_Beg_Clip.P07 = AT_data.AT_vel.P2007...
%   (AT_data.find_AT_value.P07:end);
AT_data.Btrack_Beg_Clip.P11 = AT_data.AT_vel.P2011...
  (AT_data.find_AT_value.P11:end); 
AT_data.Btrack_Beg_Clip.P13 = AT_data.AT_vel.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.Btrack_Beg_Clip.P15 = AT_data.AT_vel.P2015...
  (AT_data.find_AT_value.P15:end);
% AT_data.Btrack_Beg_Clip.P19 = AT_data.AT_vel.P2019...
%   (AT_data.find_AT_value.P19:end);

% Clipping from new start locations to a given value end element value 
% AT_data.Btrack_End_Clip.P07 = AT_data.Btrack_Beg_Clip.P07...
%   (AT_data.Btrack_Beg_Clip.P07 <= 5.79e+04);
AT_data.Btrack_End_Clip.P11 = AT_data.Btrack_Beg_Clip.P11...
  (AT_data.Btrack_Beg_Clip.P11 <= 5.79e+04);
AT_data.Btrack_End_Clip.P13 = AT_data.Btrack_Beg_Clip.P13...
  (AT_data.Btrack_Beg_Clip.P13 <= 5.79e+04);
AT_data.Btrack_End_Clip.P15 = AT_data.Btrack_Beg_Clip.P15...
  (AT_data.Btrack_Beg_Clip.P15 <= 5.79e+04); 
% AT_data.Btrack_End_Clip.P19 = AT_data.Btrack_Beg_Clip.P19...
%   (AT_data.Btrack_Beg_Clip.P19 <= 5.79e+04);

% Save along track data size as variable to see if there is any errors
%AT_data.array_size.P07_AT = size(AT_data.Btrack_End_Clip.P07);
AT_data.array_size.P11_AT = size(AT_data.Btrack_End_Clip.P11);
AT_data.array_size.P13_AT = size(AT_data.Btrack_End_Clip.P13);
AT_data.array_size.P15_AT = size(AT_data.Btrack_End_Clip.P15);
% AT_data.array_size.P19_AT = size(AT_data.Btrack_End_Clip.P19);

%% Elevation data Clipping to Section size of Along Track files
% Elevation data beginning clipping from start element in Along Track BED
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip.P2011 = AT_data.elevB.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip.P2013 = AT_data.elevB.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip.P2015 = AT_data.elevB.P2015...
  (AT_data.find_AT_value.P15:end);
% AT_data.elev_Beg_Clip.P2019 = AT_data.elevB.P2019...
%   (AT_data.find_AT_value.P19:end);

% Elevation data end clipping from end of Along track data BED
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip.P2011 = AT_data.elev_Beg_Clip.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip.P2013 = AT_data.elev_Beg_Clip.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip.P2015 = AT_data.elev_Beg_Clip.P2015...
  (1:length(AT_data.Btrack_End_Clip.P15));  
% AT_data.elev_End_Clip.P2019 = AT_data.elev_Beg_Clip.P2019...
%   (1:length(AT_data.Btrack_End_Clip.P19));

% Elevation data beginning clipping from start element in Along Track SURF  
AT_data.elev_Beg_Clip_SURF.P2011 = AT_data.elevS.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip_SURF.P2013 = AT_data.elevS.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_SURF.P2015 = AT_data.elevS.P2015...
  (AT_data.find_AT_value.P15:end);
% AT_data.elev_Beg_Clip_SURF.P2019 = AT_data.elevS.P2019...
%   (AT_data.find_AT_value.P19:end);
  
% Elevation data end clipping from end of Along track data SURF
AT_data.elev_End_Clip_SURF.P2011 = AT_data.elev_Beg_Clip_SURF.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip_SURF.P2013 = AT_data.elev_Beg_Clip_SURF.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_SURF.P2015 = AT_data.elev_Beg_Clip_SURF.P2015...
  (1:length(AT_data.Btrack_End_Clip.P15));
% AT_data.elev_End_Clip_SURF.P2019 = AT_data.elev_Beg_Clip_SURF.P2019...
%   (1:length(AT_data.Btrack_End_Clip.P19));

% Velocity data beginning clipping from start element in Along Track LAT 
AT_data.elev_Beg_Clip_LAT.P2011 = AT_data.latitudes.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip_LAT.P2013 = AT_data.latitudes.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip_LAT.P2015 = AT_data.latitudes.P2015...
  (AT_data.find_AT_value.P15:end);
% AT_data.elev_Beg_Clip_LAT.P2019 = AT_data.latitudes.P2019...
%   (AT_data.find_AT_value.P19:end);
  
% Velocity data end clipping from end of Along track data LAT
AT_data.elev_End_Clip_LAT.P2011 = AT_data.elev_Beg_Clip_LAT.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip_LAT.P2013 = AT_data.elev_Beg_Clip_LAT.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_LAT.P2015 = AT_data.elev_Beg_Clip_LAT.P2015...
  (1:length(AT_data.Btrack_End_Clip.P15));
% AT_data.elev_End_Clip_LAT.P2019 = AT_data.elev_Beg_Clip_LAT.P2019...
%   (1:length(AT_data.Btrack_End_Clip.P19));

% Velocity data beginning clipping from start element in Along Track LON 
AT_data.elev_Beg_Clip_LON.P2011 = AT_data.longitudes.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip_LON.P2013 = AT_data.longitudes.P2013...
  (AT_data.find_AT_value.P13:end); 
AT_data.elev_Beg_Clip_LON.P2015 = AT_data.longitudes.P2015...
  (AT_data.find_AT_value.P15:end);
% AT_data.elev_Beg_Clip_LON.P2019 = AT_data.longitudes.P2019...
%   (AT_data.find_AT_value.P19:end);
  
% Velocity data end clipping from end of Along track data LON
AT_data.elev_End_Clip_LON.P2011 = AT_data.elev_Beg_Clip_LON.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip_LON.P2013 = AT_data.elev_Beg_Clip_LON.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_LON.P2015 = AT_data.elev_Beg_Clip_LON.P2015...
  (1:length(AT_data.Btrack_End_Clip.P15));
% AT_data.elev_End_Clip_LON.P2019 = AT_data.elev_Beg_Clip_LON.P2019...
%   (1:length(AT_data.Btrack_End_Clip.P19));

% Velocity data beginning clipping from start element in Along Track PASS 
AT_data.elev_Beg_Clip_PASS.P2011 = AT_data.pass.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.elev_Beg_Clip_PASS.P2013 = AT_data.pass.P2013...
  (AT_data.find_AT_value.P13:end); 
AT_data.elev_Beg_Clip_PASS.P2015 = AT_data.pass.P2015...
  (AT_data.find_AT_value.P15:end);
% AT_data.elev_Beg_Clip_PASS.P2019 = AT_data.pass.P2019...
%   (AT_data.find_AT_value.P19:end);
  
% Velocity data end clipping from end of Along track data PASS
AT_data.elev_End_Clip_PASS.P2011 = AT_data.elev_Beg_Clip_PASS.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip_PASS.P2013 = AT_data.elev_Beg_Clip_PASS.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip_PASS.P2015 = AT_data.elev_Beg_Clip_PASS.P2015...
  (1:length(AT_data.Btrack_End_Clip.P15));
% AT_data.elev_End_Clip_PASS.P2019 = AT_data.elev_Beg_Clip_PASS.P2019...
%   (1:length(AT_data.Btrack_End_Clip.P19));

% Save along elevation data size as variable to see if there is any errors
% BED
%AT_data.array_size.P07_elev = size(AT_data.elev_End_Clip.P2007);
AT_data.array_size.P11_elev = size(AT_data.elev_End_Clip.P2011);
AT_data.array_size.P13_elev = size(AT_data.elev_End_Clip.P2013);
AT_data.array_size.P15_elev = size(AT_data.elev_End_Clip.P2015);
% AT_data.array_size.P19_elev = size(AT_data.elev_End_Clip.P2019);
% SURF
AT_data.array_size.P11_elevSURF = size(AT_data.elev_End_Clip_SURF.P2011);
AT_data.array_size.P13_elevSURF = size(AT_data.elev_End_Clip_SURF.P2013);
AT_data.array_size.P15_elevSURF = size(AT_data.elev_End_Clip_SURF.P2015);
% AT_data.array_size.P19_elevSURF = size(AT_data.elev_End_Clip_SURF.P2019);
% LAT
AT_data.array_size.P11_LAT = size(AT_data.elev_End_Clip_LAT.P2011);
AT_data.array_size.P13_LAT = size(AT_data.elev_End_Clip_LAT.P2013);
AT_data.array_size.P15_LAT = size(AT_data.elev_End_Clip_LAT.P2015);
% AT_data.array_size.P19_LAT = size(AT_data.elev_End_Clip_LAT.P2019);
% LON
AT_data.array_size.P11_LON = size(AT_data.elev_End_Clip_LON.P2011);
AT_data.array_size.P13_LON = size(AT_data.elev_End_Clip_LON.P2013);
AT_data.array_size.P15_LON = size(AT_data.elev_End_Clip_LON.P2015);
% AT_data.array_size.P19_LON = size(AT_data.elev_End_Clip_LON.P2019);
% PASS
AT_data.array_size.P11_PASS = size(AT_data.elev_End_Clip_PASS.P2011);
AT_data.array_size.P13_PASS = size(AT_data.elev_End_Clip_PASS.P2013);
AT_data.array_size.P15_PASS = size(AT_data.elev_End_Clip_PASS.P2015);
% AT_data.array_size.P19_PASS = size(AT_data.elev_End_Clip_PASS.P2019);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
% AT_data.query_array.P07 = (AT_data.AT_vel.P2007(1):0.1:...
%   AT_data.AT_vel.P2007(end));
AT_data.query_array.P11 = (AT_data.AT_vel.P2011(1):0.1:...
AT_data.AT_vel.P2011(end));
%AT_data.query_array.P13 = (AT_data.AT_vel.P2013(1):0.1:...
%AT_data.AT_vel.P2013(end));
%AT_data.query_array.P15 = (AT_data.AT_vel.P2015(1):0.1:...
%AT_data.AT_vel.P2015(end));
%AT_data.query_array.P19 = (AT_data.AT_vel.P2019(1):0.1:...
%AT_data.AT_vel.P2019(end));

% Apply interpolation to each profile using selected query array
% 2007
% AT_data.interp_data.P07 = interp1(AT_data.Btrack_End_Clip.P07, ...
%   AT_data.elev_End_Clip.P2007, AT_data.query_array.P07);

% 2011 BED
AT_data.interp_data.P11 = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip.P2011, AT_data.query_array.P11);
% 2013 BED
AT_data.interp_data.P13 = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip.P2013, AT_data.query_array.P11);
% 2015 BED
AT_data.interp_data.P15 = interp1(AT_data.Btrack_End_Clip.P15, ...
  AT_data.elev_End_Clip.P2015, AT_data.query_array.P11);
% 2019 BED
% AT_data.interp_data.P19 = interp1(AT_data.Btrack_End_Clip.P19, ...
%   AT_data.elev_End_Clip.P2019, AT_data.query_array.P11);

% 2011 SURF
AT_data.interp_data.P11_Surf = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip_SURF.P2011, AT_data.query_array.P11);
% 2013 SURF
AT_data.interp_data.P13_Surf = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_SURF.P2013, AT_data.query_array.P11);
% 2015 SURF
AT_data.interp_data.P15_Surf = interp1(AT_data.Btrack_End_Clip.P15, ...
  AT_data.elev_End_Clip_SURF.P2015, AT_data.query_array.P11);
% 2019 SURF
% AT_data.interp_data.P19_Surf = interp1(AT_data.Btrack_End_Clip.P19, ...
%   AT_data.elev_End_Clip_SURF.P2019, AT_data.query_array.P11);

% 2011 LAT
AT_data.interp_data.P11_LAT = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip_LAT.P2011, AT_data.query_array.P11);
% 2013 LAT
AT_data.interp_data.P13_LAT = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_LAT.P2013, AT_data.query_array.P11);
% 2015 LAT
AT_data.interp_data.P15_LAT = interp1(AT_data.Btrack_End_Clip.P15, ...
  AT_data.elev_End_Clip_LAT.P2015, AT_data.query_array.P11);
% 2019 LAT
% AT_data.interp_data.P19_LAT = interp1(AT_data.Btrack_End_Clip.P19, ...
%   AT_data.elev_End_Clip_LAT.P2019, AT_data.query_array.P11);

% 2011 LON
AT_data.interp_data.P11_LON = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip_LON.P2011, AT_data.query_array.P11);
% 2013 LON
AT_data.interp_data.P13_LON = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip_LON.P2013, AT_data.query_array.P11);
% 2015 LON
AT_data.interp_data.P15_LON = interp1(AT_data.Btrack_End_Clip.P15, ...
  AT_data.elev_End_Clip_LON.P2015, AT_data.query_array.P11);
% 2019 LON
% AT_data.interp_data.P19_LON = interp1(AT_data.Btrack_End_Clip.P19, ...
%   AT_data.elev_End_Clip_LON.P2019, AT_data.query_array.P11);

% Calculate melt rates from interpolated profile pairings
% 2007-2011 melt (Vertical difference in Features)
% AT_data.melt_rates.P07_P11 = AT_data.interp_data.P11 - ...
%   AT_data.interp_data.P07; 
  
% 2011-2013 melt (Vertical difference in Features)
AT_data.melt_rates.P11_P13 = AT_data.interp_data.P13 - ...
  AT_data.interp_data.P11;
% 2013-2015 melt (Vertical difference in Features)
AT_data.melt_rates.P13_P15 = AT_data.interp_data.P15 - ...
  AT_data.interp_data.P13;
% 2011-2015 melt (Vertical difference in Features)
AT_data.melt_rates.P11_P15 = AT_data.interp_data.P15 - ...
  AT_data.interp_data.P11;
% 2015-2019 melt (Vertical difference in Features)
% AT_data.melt_rates.P15_P19 = AT_data.interp_data.P19 - ...
%   AT_data.interp_data.P15;
% 2011-2019 melt (Vertical difference in Features)
% AT_data.melt_rates.P11_P19 = AT_data.interp_data.P19 - ...
%   AT_data.interp_data.P11;

%% Export data to csv
% Concatenate and take transpose of Lon, Lat, Surf, Bed fields. Concatenate horizontally for each year
AT_data.export.P11 = cat(2, AT_data.interp_data.P11_LON.', AT_data.interp_data.P11_LAT.', AT_data.interp_data.P11_SURF.', AT_data.interp_data.P11.' );
AT_data.export.P13 = cat(2, AT_data.interp_data.P13_LON.', AT_data.interp_data.P13_LAT.', AT_data.interp_data.P13_SURF.', AT_data.interp_data.P13.' );
AT_data.export.P15 = cat(2, AT_data.interp_data.P15_LON.', AT_data.interp_data.P15_LAT.', AT_data.interp_data.P15_SURF.', AT_data.interp_data.P15.' );
AT_data.export.P19 = cat(2, AT_data.interp_data.P19_LON.', AT_data.interp_data.P19_LAT.', AT_data.interp_data.P19_SURF.', AT_data.interp_data.P19.' );

%% Define Header Array of strings and vertically concatenate to data 
cheader = {'Lons', 'Lats', 'Surface', 'Depth'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

% change folder
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\CSV_export_files\'

%write header to file 2011
fid = fopen('Ryder1_lat_lon_surf_bed_11.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('Ryder1_lat_lon_surf_bed_11.csv', AT_data.export.P11, '-append');

%write header to file 2013
fid = fopen('Ryder1_lat_lon_surf_bed_13.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('Ryder1_lat_lon_surf_bed_13.csv', AT_data.export.P13, '-append');

%write header to file 2015
fid = fopen('Ryder1_lat_lon_surf_bed_15.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('Ryder1_lat_lon_surf_bed_15.csv', AT_data.export.P15, '-append');

%write header to file 2019
fid = fopen('Ryder1_lat_lon_surf_bed_19.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('Ryder1_lat_lon_surf_bed_19.csv', AT_data.export.P19, '-append');

%% Test Alignment Figures Section 
figure(1)
hold on
h1 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
h2 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.elev_End_Clip.P2013);
h3 = plot(AT_data.Btrack_End_Clip.P15/1e3, AT_data.elev_End_Clip.P2015);
%h4 = plot(AT_data.Btrack_End_Clip.P19/1e3, AT_data.elev_End_Clip.P2019);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2011', 'original 2013', 'original 2015', ...
  'original 2019', 'Location', 'southeast');

figure(2)
h1 = plot(AT_data.Btrack_End_Clip.P07/1e3, AT_data.elev_End_Clip.P2007);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
h3 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.elev_End_Clip.P2013);
h4 = plot(AT_data.Btrack_End_Clip.P15/1e3, AT_data.elev_End_Clip.P2015);
h5 = plot(AT_data.Btrack_End_Clip.P19/1e3, AT_data.elev_End_Clip.P2019);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2007', 'original 2011', 'original 2013', ...
  'original 2015', 'original 2019', 'Location', 'southeast');
 
%% Interpolated Figures
figure(2)
h1 = plot(AT_data.query_array.P11, AT_data.interp_data.P07);
hold on
h2 = plot(AT_data.query_array.P11, AT_data.interp_data.P11);
h3 = plot(AT_data.query_array.P11, AT_data.interp_data.P13);
h4 = plot(AT_data.query_array.P11, AT_data.interp_data.P15);
%h5 = plot(AT_data.query_array.P11, AT_data.interp_data.P19);
h7 = plot(AT_data.Btrack_End_Clip.P11, AT_data.elev_End_Clip.P2011);
h8 = plot(AT_data.Btrack_End_Clip.P13, AT_data.elev_End_Clip.P2013);
h9 = plot(AT_data.Btrack_End_Clip.P15, AT_data.elev_End_Clip.P2015);
%h10 = plot(AT_data.Btrack_End_Clip.P19, AT_data.elev_End_Clip.P2019);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2007', 'interpolated 2011', 'interpolated 2013', ...
  'interpolated 2015', 'interpolated 2019', 'old 2007','old 2011', ...
  'old 2013', 'old 2015', 'old 2019', 'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P07, AT_data.interp_data.P07);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.interp_data.P11);
h3 = plot(AT_data.query_array.P07, AT_data.interp_data.P13);
h4 = plot(AT_data.query_array.P07, AT_data.interp_data.P15);
h5 = plot(AT_data.query_array.P07, AT_data.interp_data.P19);
h6 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P19);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2007', 'interpolated 2011', 'interpolated 2013',...
  'interpolated 2015', 'interpolated 2019', '2007-2019 melt', ...
  'Location', 'southeast');

%% Melt Figures
figure(4)
h1 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P11);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.melt_rates.P11_P13);
h3 = plot(AT_data.query_array.P07, AT_data.melt_rates.P13_P15);
h4 = plot(AT_data.query_array.P07, AT_data.melt_rates.P15_P19);
h5 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P19);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2007-2011 melt', '2011-2013 melt', '2013-2015 melt', ...
  '2015-2019 melt', '2007-2019 melt', 'Location', 'southeast');

% Subplot figure of melt rates
figure(5)
subplot(5,1,1);
plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P11);
title('Melt 2007-2011');

subplot(5,1,2);
plot(AT_data.query_array.P07, AT_data.melt_rates.P11_P13);
title('Melt 2011-2013');

subplot(5,1,3);
plot(AT_data.query_array.P07, AT_data.melt_rates.P13_P15);
title('Melt 2013-2015');

subplot(5,1,4);
plot(AT_data.query_array.P07, AT_data.melt_rates.P15_P19);
title('Melt 2015-2019');

subplot(5,1,5);
plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P19);
title('Melt 2007-2019');

% Annual melt rate average based off 2007-2019 melt subtraction
% figure(6)
% plot(AT_data.query_array.P07/1e3, AT_data.melt_rates.P07_P19/12, 'color','g');
% title('Melt 2013-2019');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');
