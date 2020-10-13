%% Flightline Data Interpolation - Petermann Line 1 
% Years: 2011, 2014, 2015, 2017, 2018, 2019
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
AT_data.elevB.P2011 = (pass(1).layers(2).layer_elev);
AT_data.elevB.P2014 = (pass(2).layers(2).layer_elev);
%AT_data.elevB.P2015 = pass(3).layers(2).layer_elev);
%AT_data.elevB.P2017 = pass(4).layers(2).layer_elev);
AT_data.elevB.P2018 = (pass(5).layers(2).layer_elev);
%AT_data.elevB.P2019 = pass(6).layers(2).layer_elev);

% Save Surface profiles
AT_data.elevS.P2011 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2014 = (pass(2).layers(1).layer_elev);
%AT_data.elevS.P2015 = (pass(3).layers(1).layer_elev);
%AT_data.elevS.P2017 = (pass(4).layers(1).layer_elev);
AT_data.elevS.P2018 = (pass(5).layers(1).layer_elev);
%AT_data.elevS.P2019 = (pass(6).layers(1).layer_elev);

% Save annual alongtrack profile data
AT_data.pass.P2011 = pass(1).along_track;
AT_data.pass.P2014 = pass(2).along_track;
%AT_data.pass.P2015 = pass(3).along_track;
%AT_data.pass.P2017 = pass(4).along_track;
AT_data.pass.P2018= pass(5).along_track;
%AT_data.pass.P2019 = pass(6).along_track;

% Save annual velocity correction data
AT_data.vel.P2011 = pass(1).vel;
AT_data.vel.P2014 = pass(2).vel;
%AT_data.vel.P2015 = pass(3).vel;
%AT_data.vel.P2017 = pass(4).vel;
AT_data.vel.P2018 = pass(5).vel;
%AT_data.vel.P2016 = pass(6).vel;

% Save Velocity Corrected Along_track data
AT_data.AT_vel.P2011 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2014 = pass(baseline_master_idx).along_track + pass(2).vel;
%AT_data.AT_vel.P2015 = pass(baseline_master_idx).along_track + pass(3).vel;
%AT_data.AT_vel.P2017 = pass(baseline_master_idx).along_track + pass(4).vel;
AT_data.AT_vel.P2018 = pass(baseline_master_idx).along_track + pass(5).vel;
%AT_data.AT_vel.P2019 = pass(baseline_master_idx).along_track + pass(6).vel;

%% Section 2 - Mask, indexing, and clipping of along track
  
% Locate Along Track start element in each profile, save element ID as a
% variable for the clipping
% 2011
AT_data.Btrack.P11 = AT_data.AT_vel.P2011(AT_data.AT_vel.P2011 >= ...
  AT_data.AT_vel.P2011(1));
AT_data.find_AT_value.P11 = find(AT_data.AT_vel.P2011 == ...
  AT_data.Btrack.P11(1));
  
% 2014  
AT_data.Btrack.P14 = AT_data.AT_vel.P2014(AT_data.AT_vel.P2014 >= ...
  AT_data.AT_vel.P2011(1));
AT_data.find_AT_value.P14 = find(AT_data.AT_vel.P2014 == ...
  AT_data.Btrack.P14(1));

% 2015 
% AT_data.Btrack.P15 = AT_data.AT_vel.P2015(AT_data.AT_vel.P2015 >= ...
%   AT_data.AT_vel.P2011(1));
% AT_data.find_AT_value.P15 = find(AT_data.AT_vel.P2015 == ...
%   AT_data.Btrack.P15(1));

% 2017 
% AT_data.Btrack.P17 = AT_data.AT_vel.P2017(AT_data.AT_vel.P2017 >= ...
%   AT_data.AT_vel.P2011(1));
% AT_data.find_AT_value.P17 = find(AT_data.AT_vel.P2017 == ...
%   AT_data.Btrack.P17(1));

% 2018 
AT_data.Btrack.P18 = AT_data.AT_vel.P2018(AT_data.AT_vel.P2018 >= ...
  AT_data.AT_vel.P2011(1));
AT_data.find_AT_value.P18 = find(AT_data.AT_vel.P2018 == ...
  AT_data.Btrack.P18(1));
  
% 2019 
% AT_data.Btrack.P19 = AT_data.AT_vel.P2019(AT_data.AT_vel.P2019 >= ...
%   AT_data.AT_vel.P2011(1));
% AT_data.find_AT_value.P19 = find(AT_data.AT_vel.P2019 == ...
%   AT_data.Btrack.P19(1));

% Clipping from start point in each profile to the end of the profile
AT_data.Btrack_Beg_Clip.P11 = AT_data.AT_vel.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.Btrack_Beg_Clip.P14 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end);
%AT_data.Btrack_Beg_Clip.P15 = AT_data.AT_vel.P2015...
%  (AT_data.find_AT_value.P15:end);
%AT_data.Btrack_Beg_Clip.P17 = AT_data.AT_vel.P2017...
%  (AT_data.find_AT_value.P17:end);
AT_data.Btrack_Beg_Clip.P18 = AT_data.AT_vel.P2018...
  (AT_data.find_AT_value.P18:end);
%AT_data.Btrack_Beg_Clip.P19 = AT_data.AT_vel.P2019...
%  (AT_data.find_AT_value.P19:end);

% Clipping from new start locations to a given value end element value 
AT_data.Btrack_End_Clip.P11 = AT_data.Btrack_Beg_Clip.P11...
  (AT_data.Btrack_Beg_Clip.P11 <= 5.79e+04);
AT_data.Btrack_End_Clip.P14 = AT_data.Btrack_Beg_Clip.P14...
  (AT_data.Btrack_Beg_Clip.P14 <= 5.79e+04);
%AT_data.Btrack_End_Clip.P15 = AT_data.Btrack_Beg_Clip.P15...
%  (AT_data.Btrack_Beg_Clip.P15 <= 5.79e+04);
%AT_data.Btrack_End_Clip.P17 = AT_data.Btrack_Beg_Clip.P17...
%  (AT_data.Btrack_Beg_Clip.P17 <= 5.79e+04);
AT_data.Btrack_End_Clip.P18 = AT_data.Btrack_Beg_Clip.P18...
  (AT_data.Btrack_Beg_Clip.P18 <= 5.79e+04);
%AT_data.Btrack_End_Clip.P19 = AT_data.Btrack_Beg_Clip.P19...
%  (AT_data.Btrack_Beg_Clip.P19 <= 5.79e+04);
  
%Save along track data size as variable to see if there is any errors
AT_data.array_size.P11_AT = size(AT_data.Btrack_End_Clip.P11);
AT_data.array_size.P14_AT = size(AT_data.Btrack_End_Clip.P14);
%AT_data.array_size.P15_AT = size(AT_data.Btrack_End_Clip.P15);
%AT_data.array_size.P17_AT = size(AT_data.Btrack_End_Clip.P17);
AT_data.array_size.P18_AT = size(AT_data.Btrack_End_Clip.P18);
%AT_data.array_size.P19_AT = size(AT_data.Btrack_End_Clip.P19);

%% For terminus End element trimming (NEEDS WORK)
%AT_data.Btrack_E.P11 = AT_data.Btrack_C.P11(AT_data.Btrack_C.P11 <= ...
%AT_data.Btrack_C.P11(end));
%AT_data.find_AT_value_E.P11 = find(AT_data.Btrack_C.P11 == ...
%AT_data.Btrack_C.P11(end));
%     
%AT_data.Btrack_E.P14 = AT_data.Btrack_C.P14(AT_data.Btrack_C.P14 <= ...
%AT_data.Btrack_C.P11(end));
%AT_data.find_AT_value_E.P14 = find(AT_data.Btrack_C.P14 == ...
%AT_data.Btrack_C.P11(end));
%  
%AT_data.Btrack_E.P18 = AT_data.Btrack_C.P18(AT_data.Btrack_C.P18 <= ...
%AT_data.Btrack_C.P11(end));
%AT_data.find_AT_value_E.P18 = find(AT_data.Btrack_C.P18 == ...
%AT_data.Btrack_C.P11(end));
  
%AT_data.Btrack_CE.P11 = AT_data.Btrack_C.P11...
%(1:AT_data.find_AT_value_E.P11);
%AT_data.Btrack_CE.P14 = AT_data.Btrack_C.P14...
%(1:AT_data.find_AT_value_E.P14); 
%AT_data.Btrack_CE.P18 = AT_data.Btrack_C.P18...
%(1:AT_data.find_AT_value_E.P18);
%   
%AT_data.Etrack_C.P11 = AT_data.Btrack_C.P11(1:length(AT_data.Btrack_C.P18));
%AT_data.Etrack_C.P14 = AT_data.Btrack_C.P14(1:length(AT_data.Btrack_C.P18));
%AT_data.Etrack_C.P18 = AT_data.Btrack_C.P18(1:length(AT_data.Btrack_C.P18));

% Along track end clipping mask to the shortest profile (physically)
%AT_data.Etrack.P11 = AT_data.Btrack.P11(AT_data.Btrack.P11 ...
%<= AT_data.Btrack.P14(end));
%AT_data.Etrack.P14 = AT_data.Btrack.P14(AT_data.Btrack.P14 ...
%<= AT_data.Btrack.P14(end));
%AT_data.Etrack.P18 = AT_data.Btrack.P18(AT_data.Btrack.P18 ...
%<= AT_data.Btrack.P14(end));

%% Elevation data Clipping to Section size of Along Track files
% Elevation data beginning clipping from start element in Along Track
AT_data.elev_Beg_Clip.P2011 = AT_data.elevB.P2011...
  (AT_data.find_AT_value.P11:end); 
AT_data.elev_Beg_Clip.P2014 = AT_data.elevB.P2014...
  (AT_data.find_AT_value.P14:end);
%AT_data.elev_Beg_Clip.P2015 = AT_data.elevB.P2015...
%  (AT_data.find_AT_value.P15:end);
%AT_data.elev_Beg_Clip.P2017 = AT_data.elevB.P2017...
%  (AT_data.find_AT_value.P17:end);
AT_data.elev_Beg_Clip.P2018 = AT_data.elevB.P2018...
  (AT_data.find_AT_value.P18:end);
%AT_data.elev_Beg_Clip.P2019 = AT_data.elevB.P2019...
%  (AT_data.find_AT_value.P19:end);  

% Elevation data end clipping from end of Along track data
AT_data.elev_End_Clip.P2011 = AT_data.elev_Beg_Clip.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip.P2014 = AT_data.elev_Beg_Clip.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
%AT_data.elev_End_Clip.P2015 = AT_data.elev_Beg_Clip.P2015...
%  (1:length(AT_data.Btrack_End_Clip.P15));
%AT_data.elev_End_Clip.P2017 = AT_data.elev_Beg_Clip.P2017...
%  (1:length(AT_data.Btrack_End_Clip.P17));
AT_data.elev_End_Clip.P2018 = AT_data.elev_Beg_Clip.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));
%AT_data.elev_End_Clip.P2019 = AT_data.elev_Beg_Clip.P2019...
%  (1:length(AT_data.Btrack_End_Clip.P19));

% Save along elevation data size as variable to see if there is any errors
AT_data.array_size.P11_elev = size(AT_data.elev_End_Clip.P2011);
AT_data.array_size.P14_elev = size(AT_data.elev_End_Clip.P2014);
%AT_data.array_size.P15_elev = size(AT_data.elev_End_Clip.P2015);
%AT_data.array_size.P17_elev = size(AT_data.elev_End_Clip.P2017);
AT_data.array_size.P18_elev = size(AT_data.elev_End_Clip.P2018);
%AT_data.array_size.P19_elev = size(AT_data.elev_End_Clip.P2019);

%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
AT_data.query_array.P11 = (AT_data.AT_vel.P2011(1):0.1:...
  AT_data.AT_vel.P2011(end));
%AT_data.query_array.P14 = (AT_data.AT_vel.P2014(1):0.1:...
%AT_data.AT_vel.P2014(end));
%AT_data.query_array.P15 = (AT_data.AT_vel.P2015(1):0.1:...
%AT_data.AT_vel.P2015(end));
%AT_data.query_array.P17 = (AT_data.AT_vel.P2017(1):0.1:...
%AT_data.AT_vel.P2017(end));
%AT_data.query_array.P18 = (AT_data.AT_vel.P2018(1):0.1:...
%AT_data.AT_vel.P2018(end));
%AT_data.query_array.P19 = (AT_data.AT_vel.P2019(1):0.1:...
%AT_data.AT_vel.P2019(end));

% Apply interpolation to each profile using selected query array
% 2011
AT_data.interp_data.P11 = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip.P2011, AT_data.query_array.P11);
 
% 2014
AT_data.interp_data.P14 = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip.P2014, AT_data.query_array.P11);

% 2015
%AT_data.interp_data.P15 = interp1(AT_data.Btrack_End_Clip.P15, ...
%  AT_data.elev_End_Clip.P2015, AT_data.query_array.P11);

% 2017
%AT_data.interp_data.P17 = interp1(AT_data.Btrack_End_Clip.P17, ...
%  AT_data.elev_End_Clip.P2017, AT_data.query_array.P11);

% 2018
AT_data.interp_data.P18 = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip.P2018, AT_data.query_array.P11);

% 2019
%AT_data.interp_data.P19 = interp1(AT_data.Btrack_End_Clip.P19, ...
%  AT_data.elev_End_Clip.P2019, AT_data.query_array.P11);

% Calculate melt rates from interpolated profile pairings
% 2011-2014 melt (Vertical difference in Features)
AT_data.melt_rates.P11_P14 = AT_data.interp_data.P14 - ...
  AT_data.interp_data.P11; 
% 2014-2015 melt (Vertical difference in Features)
% AT_data.melt_rates.P14_P15 = AT_data.interp_data.P15 - ...
%   AT_data.interp_data.P14;
% 2015-2017 melt (Vertical difference in Features)
% AT_data.melt_rates.P15_P17 = AT_data.interp_data.P17 - ...
%   AT_data.interp_data.P17;
% 2017-2018 melt (Vertical difference in Features)
% AT_data.melt_rates.P17_P18 = AT_data.interp_data.P18 - ...
%   AT_data.interp_data.P17;
% 2018-2019 melt (Vertical difference in Features)
% AT_data.melt_rates.P18_P19 = AT_data.interp_data.P19 - ...
%   AT_data.interp_data.P18;
% 2011-2018 melt (Vertical difference in Features)
% AT_data.melt_rates.P11_P19 = AT_data.interp_data.P19 - ...
%   AT_data.interp_data.P11;
  
% THIS IS HOW TO DO THE INTERPOLATION
%MAKES ALONG TRACK SPACING IN THE PROFILE ALONG TRACK FILE EVERY 0.1M
%x_t11 = (AT_data.AT_vel.P2011(1):0.1:AT_data.AT_vel.P2011(end));
%INTERPOLATION BETWEEN THE SPACED LINE AND NEW YEARS
%FOR 2011
%elev_11_int = interp1(AT_data.Btrack_E.P11, AT_data.elev_EndC.P2011,x_t11);
%FOR 2014
%elev_14_int = interp1(AT_data.Btrack_E.P14, AT_data.elev_EndC.P2014,x_t11);
%FOR 2018
%elev_18_int = interp1(AT_data.Btrack_E.P18, AT_data.elev_EndC.P2018,x_t11);
%% Test Figures Section 

figure(1)
h1 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
%h3 = plot(AT_data.Btrack_End_Clip.P15/1e3, AT_data.elev_End_Clip.P2015);
%h4 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
h5 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
%h6 = plot(AT_data.Btrack_End_Clip.P19/1e3, AT_data.elev_End_Clip.P2019);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2011', 'original 2014', 'original 2015', ...
  'original 2017', 'original 2018', 'original 2019', ...
  'Location', 'southeast');
  
figure(2)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11);
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14);
%h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P15);
%h4 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P17);
h5 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18);
%h6 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P19);
h7 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
h8 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
%h9 = plot(AT_data.Btrack_End_Clip.P15/1e3, AT_data.elev_End_Clip.P2015);
%h10 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
h11 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
%h12 = plot(AT_data.Btrack_End_Clip.P19/1e3, AT_data.elev_End_Clip.P2019);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2011', 'interpolated 2014', 'interpolated 2015', ...
  'interpolated 2017', 'interpolated 2018', 'interpolated 2019', ...
  'old 2011','old 2014', 'old 2015', 'old 2017', 'old 2018', 'old 2019', ...
  'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11);
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14);
%h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P15)
%h4 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P17)
h5 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18);
%h6 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P19)
%h7 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P19);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2011', 'interpolated 2014', 'interpolated 2015',...
  'interpolated 2017', 'interpolated 2018', 'interpolated 2019', ...
  '2011-2019 melt', 'Location', 'southeast');
  
figure(4)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P14);
hold on
%h2 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P14_P15);
%h3 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P15_P17);
%h4 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P17_P18);
%h5 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P18_P19);
%h6 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P19);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2011-2014 melt', '2014-2018 melt', '2011-2018 melt', ...
  'Location', 'southeast');

% Subplot of melt rates for paired years
figure(5)
subplot(6,1,1);
plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P14, 'color','m');
title('Melt 2011-2014');
ylim([-100, 400]);

% subplot(6,1,2);
% plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P14_P15, 'color','c');
% title('Melt 2014-2015');
% ylim([-100, 400]);
% ylabel('Elevation Change (m)');
% 
% subplot(6,1,3);
% plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P15_P17, 'color','g');
% title('Melt 2015-2017');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');
% 
% subplot(6,1,4);
% plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P17_P18, 'color','g');
% title('Melt 2017-2018');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');
% 
% subplot(6,1,5);
% plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P18_P19, 'color','g');
% title('Melt 2018-2019');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');
%
% subplot(6,1,6);
% plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P19, 'color','g');
% title('Melt 2011-2019');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');


% Annual melt rate average based off 2011-2019 melt subtraction
% figure(6)
% plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P19/8, 'color','g');
% title('Melt 2016-2018');
% ylim([-100, 400]);
% xlabel('Along Track Distance (km)');
