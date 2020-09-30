  
%% Flightline Data Interpolation - Petermann Line 3
% Years: 2007, 2010A, 2010B, 2017
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
AT_data.elevB.P2007 = (pass(1).layers(2).layer_elev);
AT_data.elevB.P2010A = (pass(2).layers(2).layer_elev);
AT_data.elevB.P2010B = (pass(3).layers(2).layer_elev);
AT_data.elevB.P2017 = (pass(4).layers(2).layer_elev);

%Save Surface profiles
AT_data.elevS.P2007 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2010A = (pass(2).layers(1).layer_elev);
AT_data.elevS.P2010B = (pass(3).layers(1).layer_elev);
AT_data.elevS.P2017 = (pass(4).layers(1).layer_elev);

%Save annual alongtrack profile data
AT_data.pass.P2007 = pass(1).along_track;
AT_data.pass.P2010A = pass(2).along_track;
AT_data.pass.P2010B = pass(3).along_track;
AT_data.pass.P2017 = pass(4).along_track;

%Save annual velocity correction data
AT_data.vel.P2007 = pass(1).vel;
AT_data.vel.P2010A = pass(2).vel;
AT_data.vel.P2010B = pass(3).vel;
AT_data.vel.P2017 = pass(4).vel;

%Save Velocity Corrected Along_track data
AT_data.AT_vel.P2007 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2010A = pass(baseline_master_idx).along_track + pass(2).vel;
AT_data.AT_vel.P2010B = pass(baseline_master_idx).along_track + pass(3).vel;
AT_data.AT_vel.P2017 = pass(baseline_master_idx).along_track + pass(4).vel;

%% Section 2 - Mask, indexing, and clipping of along track
  
% Locate Along Track start element in each profile, save element ID as a
% variable for the clipping
% 2007
AT_data.Btrack.P07 = AT_data.AT_vel.P2007(AT_data.AT_vel.P2007 >= ...
  AT_data.AT_vel.P2007(1));
AT_data.find_AT_value.P07 = find(AT_data.AT_vel.P2007 == ...
  AT_data.Btrack.P07(1));

% 2010A  
AT_data.Btrack.P10A = AT_data.AT_vel.P2010A(AT_data.AT_vel.P2010A >= ...
  AT_data.AT_vel.P2007(1));
AT_data.find_AT_value.P10A = find(AT_data.AT_vel.P2010A == ...
  AT_data.Btrack.P10A(1));

% 2010B
AT_data.Btrack.P10B = AT_data.AT_vel.P2010B(AT_data.AT_vel.P2010B >= ...
  AT_data.AT_vel.P2007(1));
AT_data.find_AT_value.P10B = find(AT_data.AT_vel.P2010B == ...
  AT_data.Btrack.P10B(1));
  
 2017 
 AT_data.Btrack.P17 = AT_data.AT_vel.P2017(AT_data.AT_vel.P2017 >= ...
   AT_data.AT_vel.P2007(1));
 AT_data.find_AT_value.P17 = find(AT_data.AT_vel.P2017 == ... 
   AT_data.Btrack.P17(1));

% Clipping from start point in each profile to the end of the profile
AT_data.Btrack_Beg_Clip.P07 = AT_data.AT_vel.P2007...
  (AT_data.find_AT_value.P07:end);
AT_data.Btrack_Beg_Clip.P10A = AT_data.AT_vel.P2010A...
  (AT_data.find_AT_value.P10A:end); 
AT_data.Btrack_Beg_Clip.P10B = AT_data.AT_vel.P2010B...
  (AT_data.find_AT_value.P10B:end);
AT_data.Btrack_Beg_Clip.P17 = AT_data.AT_vel.P2017...
  (AT_data.find_AT_value.P17:end);

% Clipping from new start locations to a given value end element value 
AT_data.Btrack_End_Clip.P07 = AT_data.Btrack_Beg_Clip.P07...
  (AT_data.Btrack_Beg_Clip.P07 <= 5.79e+04);
AT_data.Btrack_End_Clip.P10A = AT_data.Btrack_Beg_Clip.P10A...
  (AT_data.Btrack_Beg_Clip.P10A <= 5.79e+04);
AT_data.Btrack_End_Clip.P10B = AT_data.Btrack_Beg_Clip.P10B...
  (AT_data.Btrack_Beg_Clip.P10B <= 5.79e+04);
AT_data.Btrack_End_Clip.P17 = AT_data.Btrack_Beg_Clip.P17...
  (AT_data.Btrack_Beg_Clip.P17 <= 5.79e+04); 

% Save along track data size as variable to see if there is any errors
AT_data.array_size.P07_AT = size(AT_data.Btrack_End_Clip.P07);
AT_data.array_size.P10A_AT = size(AT_data.Btrack_End_Clip.P10A);
AT_data.array_size.P10B_AT = size(AT_data.Btrack_End_Clip.P10B);
AT_data.array_size.P17_AT = size(AT_data.Btrack_End_Clip.P17);

%% Elevation data Clipping to Section size of Along Track files
% Elevation data beginning clipping from start element in Along Track
AT_data.elev_Beg_Clip.P2007 = AT_data.elevB.P2007...
  (AT_data.find_AT_value.P07:end); 
AT_data.elev_Beg_Clip.P2010A = AT_data.elevB.P2010A...
  (AT_data.find_AT_value.P10A:end);
AT_data.elev_Beg_Clip.P2010B = AT_data.elevB.P2010B...
  (AT_data.find_AT_value.P10B:end);
AT_data.elev_Beg_Clip.P2017 = AT_data.elevB.P2017...
  (AT_data.find_AT_value.P17:end);

% Elevation data end clipping from end of Along track data
AT_data.elev_End_Clip.P2007 = AT_data.elev_Beg_Clip.P2007...
  (1:length(AT_data.Btrack_End_Clip.P07));
AT_data.elev_End_Clip.P2010A = AT_data.elev_Beg_Clip.P2010A...
  (1:length(AT_data.Btrack_End_Clip.P10A));
AT_data.elev_End_Clip.P2010B = AT_data.elev_Beg_Clip.P2010B...
  (1:length(AT_data.Btrack_End_Clip.P10B));
AT_data.elev_End_Clip.P2017 = AT_data.elev_Beg_Clip.P2017...
  (1:length(AT_data.Btrack_End_Clip.P17));  

% Save along elevation data size as variable to see if there is any errors
AT_data.array_size.P07_elev = size(AT_data.elev_End_Clip.P2007);
AT_data.array_size.P10A_elev = size(AT_data.elev_End_Clip.P2010A);
AT_data.array_size.P10B_elev = size(AT_data.elev_End_Clip.P2010B);
AT_data.array_size.P17_elev = size(AT_data.elev_End_Clip.P2017);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
AT_data.query_array.P07 = (AT_data.AT_vel.P2007(1):0.1:...
  AT_data.AT_vel.P2007(end));
%AT_data.query_array.P10A = (AT_data.AT_vel.P2010A(1):0.1:...
%AT_data.AT_vel.P2010A(end));
%AT_data.query_array.P10B = (AT_data.AT_vel.P2010B(1):0.1:...
%AT_data.AT_vel.P2010B(end));
%AT_data.query_array.P17 = (AT_data.AT_vel.P2017(1):0.1:...
%AT_data.AT_vel.P2017(end));

% Apply interpolation to each profile using selected query array
% 2007
AT_data.interp_data.P07 = interp1(AT_data.Btrack_End_Clip.P07, ...
  AT_data.elev_End_Clip.P2007, AT_data.query_array.P07);

% 2010A
AT_data.interp_data.P10A = interp1(AT_data.Btrack_End_Clip.P10A, ...
  AT_data.elev_End_Clip.P2010A, AT_data.query_array.P07);

% 2010B
AT_data.interp_data.P10B = interp1(AT_data.Btrack_End_Clip.P10B, ...
  AT_data.elev_End_Clip.P2010B, AT_data.query_array.P07);

% 2017
AT_data.interp_data.P17 = interp1(AT_data.Btrack_End_Clip.P17, ...
  AT_data.elev_End_Clip.P2017, AT_data.query_array.P07);

% Calculate melt rates from interpolated profile pairings
% 2007-2010A melt (Vertical difference in Features)
AT_data.melt_rates.P07_P10A = AT_data.interp_data.P10A - ...
  AT_data.interp_data.P07; 
  
% 2010A-2010B melt (Vertical difference in Features)
AT_data.melt_rates.P10A_P10B = AT_data.interp_data.P10B - ...
  AT_data.interp_data.P10A;
  
% 2010B-2017 melt (Vertical difference in Features)
AT_data.melt_rates.P10B_P17 = AT_data.interp_data.P17 - ...
  AT_data.interp_data.P10B;

% 2007-2017 melt (Vertical difference in Features)
AT_data.melt_rates.P07_P17 = AT_data.interp_data.P17 - ...
  AT_data.interp_data.P07;
  
%% Test Figures Section 

figure(1)
h1 = plot(AT_data.Btrack_End_Clip.P07/1e3, AT_data.elev_End_Clip.P2007);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P10A/1e3, AT_data.elev_End_Clip.P2010A);
h3 = plot(AT_data.Btrack_End_Clip.P10B/1e3, AT_data.elev_End_Clip.P2010B);
h4 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2007', 'original 2010A', 'original 2010B', ...
  'original 2017', 'Location', 'southeast');
  
figure(2)
h1 = plot(AT_data.query_array.P07, AT_data.interp_data.P07);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.interp_data.P10A);
h3 = plot(AT_data.query_array.P07, AT_data.interp_data.P10B);
h4 = plot(AT_data.query_array.P07, AT_data.interp_data.P17);
h5 = plot(AT_data.Btrack_End_Clip.P07, AT_data.elev_End_Clip.P2007);
h6 = plot(AT_data.Btrack_End_Clip.P10A, AT_data.elev_End_Clip.P2010A);
h7 = plot(AT_data.Btrack_End_Clip.P10B, AT_data.elev_End_Clip.P2010B);
h8 = plot(AT_data.Btrack_End_Clip.P17, AT_data.elev_End_Clip.P2017);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2007', 'interpolated 2010A', 'interpolated 2010B', ...
  'interpolated 2017', 'old 2007','old 2010A','old 2010B', 'old 2017',  ...
  'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P07, AT_data.interp_data.P07);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.interp_data.P10A);
h3 = plot(AT_data.query_array.P07, AT_data.interp_data.P10B);
h4 = plot(AT_data.query_array.P07, AT_data.interp_data.P17);
h5 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P17);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2007', 'interpolated 2010A', 'interpolated 2010B',...
  'interpolated 2017', '2007-2017 melt', 'Location', 'southeast');
  
figure(4)
h1 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P10A);
hold on
h2 = plot(AT_data.query_array.P07, AT_data.melt_rates.P10A_P10B);
h3 = plot(AT_data.query_array.P07, AT_data.melt_rates.P10B_P17);
h4 = plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P17);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2007-2010A melt', '2010A-2010B melt', '2010B-2017 melt', ...
  '2007-2017 melt', 'Location', 'southeast');

% Subplot figure of melt rates
figure(5)
subplot(4,1,1);
plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P10A);
title('Melt 2007-2010A');

subplot(4,1,2);
plot(AT_data.query_array.P07, AT_data.melt_rates.P10A_P10B);
title('Melt 2010A-2010B');

subplot(4,1,3);
plot(AT_data.query_array.P07, AT_data.melt_rates.P10B_P17);
title('Melt 2010B-2018');

subplot(4,1,4);
plot(AT_data.query_array.P07, AT_data.melt_rates.P07_P17);
title('Melt 2007-2017');

% Annual melt rate average based off 2007-2019 melt subtraction
% figure(6)
% plot(AT_data.query_array.P07/1e3, AT_data.melt_rates.P07_P17/10, 'color','g');
% title('Melt 2007-2017');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');
