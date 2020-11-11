  
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
  
2017 
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
% Elevation data beginning clipping from start element in Along Track
% AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
%   (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip.P2013 = AT_data.elevB.P2013...
  (AT_data.find_AT_value.P13:end);
AT_data.elev_Beg_Clip.P2014 = AT_data.elevB.P2014...
  (AT_data.find_AT_value.P14:end);
AT_data.elev_Beg_Clip.P2017 = AT_data.elevB.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data
% AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
%   (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip.P2013 = AT_data.elev_Beg_Clip.P2013...
  (1:length(AT_data.Btrack_End_Clip.P13));
AT_data.elev_End_Clip.P2014 = AT_data.elev_Beg_Clip.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip.P2017 = AT_data.elev_Beg_Clip.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));  

% Save along elevation data size as variable to see if there is any errors
AT_data.array_size.P07_elev = size(AT_data.elev_End_Clip.P2007);
AT_data.array_size.P13_elev = size(AT_data.elev_End_Clip.P2013);
AT_data.array_size.P14_elev = size(AT_data.elev_End_Clip.P2014);
AT_data.array_size.P17_elev = size(AT_data.elev_End_Clip.P2017);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
% AT_data.query_array.P07 = (AT_data.AT_vel.P2007(1):0.1:...
%   AT_data.AT_vel.P2007(end));
AT_data.query_array.P13 = (AT_data.AT_vel.P2013(1):0.1:...
AT_data.AT_vel.P2013(end));
%AT_data.query_array.P14 = (AT_data.AT_vel.P2014(1):0.1:...
%AT_data.AT_vel.P2014(end));
%AT_data.query_array.P17 = (AT_data.AT_vel.P2017(1):0.1:...
%AT_data.AT_vel.P2017(end));

% Apply interpolation to each profile using selected query array
% 2007
% AT_data.interp_data.P07 = interp1(AT_data.Btrack_End_Clip.P07, ...
%   AT_data.elev_End_Clip.P2007, AT_data.query_array.P07);

% 2013
AT_data.interp_data.P13 = interp1(AT_data.Btrack_End_Clip.P13, ...
  AT_data.elev_End_Clip.P2013, AT_data.query_array.P13);

% 2014
AT_data.interp_data.P14 = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip.P2014, AT_data.query_array.P13);

% 2017
AT_data.interp_data.P17 = interp1(AT_data.Btrack_End_Clip.P17, ...
 AT_data.elev_End_Clip.P2017, AT_data.query_array.P13);

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
  
%% Test Figures Section 

figure(1)
h1 = plot(AT_data.Btrack_End_Clip.P07/1e3, AT_data.elev_End_Clip.P2007);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P13/1e3, AT_data.elev_End_Clip.P2013);
h3 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
h4 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2007', 'original 2013', 'original 2014', ...
  'original 2017', 'Location', 'southeast');
  
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
