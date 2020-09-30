%% Flightline Data Interpolation - Petermann Line 1 - 2011, 2014, 2018
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
AT_data.elevB.P2018 = (pass(3).layers(2).layer_elev);

% Save Surface profiles
AT_data.elevS.P2011 = (pass(1).layers(1).layer_elev);
AT_data.elevS.P2014 = (pass(2).layers(1).layer_elev);
AT_data.elevS.P2018 = (pass(3).layers(1).layer_elev);

% Save annual alongtrack profile data
AT_data.pass.P2011 = pass(1).along_track;
AT_data.pass.P2014 = pass(2).along_track;
AT_data.pass.P2018= pass(3).along_track;
 
% Save annual velocity correction data
AT_data.vel.P2011 = pass(1).vel;
AT_data.vel.P2014 = pass(2).vel;
AT_data.vel.P2018 = pass(3).vel;

% Save Velocity Corrected Along_track data
AT_data.AT_vel.P2011 = pass(baseline_master_idx).along_track + pass(1).vel;
AT_data.AT_vel.P2014 = pass(baseline_master_idx).along_track + pass(2).vel;
AT_data.AT_vel.P2018 = pass(baseline_master_idx).along_track + pass(3).vel;
  
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
  
% 2018 
AT_data.Btrack.P18 = AT_data.AT_vel.P2018(AT_data.AT_vel.P2018 >= ...
  AT_data.AT_vel.P2011(1));
AT_data.find_AT_value.P18 = find(AT_data.AT_vel.P2018 == ...
  AT_data.Btrack.P18(1));
  
% Clipping from start point in each profile to the end of the profile
AT_data.Btrack_Beg_Clip.P11 = AT_data.AT_vel.P2011...
  (AT_data.find_AT_value.P11:end);
AT_data.Btrack_Beg_Clip.P14 = AT_data.AT_vel.P2014...
  (AT_data.find_AT_value.P14:end); 
AT_data.Btrack_Beg_Clip.P18 = AT_data.AT_vel.P2018...
  (AT_data.find_AT_value.P18:end);

% Clipping from new start locations to a given value end element value 
AT_data.Btrack_End_Clip.P11 = AT_data.Btrack_Beg_Clip.P11...
  (AT_data.Btrack_Beg_Clip.P11 <= 5.79e+04);
AT_data.Btrack_End_Clip.P14 = AT_data.Btrack_Beg_Clip.P14...
  (AT_data.Btrack_Beg_Clip.P14 <= 5.79e+04);
AT_data.Btrack_End_Clip.P18 = AT_data.Btrack_Beg_Clip.P18...
  (AT_data.Btrack_Beg_Clip.P18 <= 5.79e+04);
  
%Save along track data size as variable to see if there is any errors
AT_data.array_size.P11_AT = size(AT_data.Btrack_End_Clip.P11);
AT_data.array_size.P14_AT = size(AT_data.Btrack_End_Clip.P14);
AT_data.array_size.P18_AT = size(AT_data.Btrack_End_Clip.P18);
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
AT_data.elev_Beg_Clip.P2018 = AT_data.elevB.P2018...
  (AT_data.find_AT_value.P18:end);
  
% Elevation data end clipping from end of Along track data
AT_data.elev_End_Clip.P2011 = AT_data.elev_Beg_Clip.P2011...
  (1:length(AT_data.Btrack_End_Clip.P11));
AT_data.elev_End_Clip.P2014 = AT_data.elev_Beg_Clip.P2014...
  (1:length(AT_data.Btrack_End_Clip.P14));
AT_data.elev_End_Clip.P2018 = AT_data.elev_Beg_Clip.P2018...
  (1:length(AT_data.Btrack_End_Clip.P18));
  
% Save along elevation data size as variable to see if there is any errors
AT_data.array_size.P11_elev = size(AT_data.elev_End_Clip.P2011);
AT_data.array_size.P14_elev = size(AT_data.elev_End_Clip.P2014);
AT_data.array_size.P18_elev = size(AT_data.elev_End_Clip.P2018);
%% Derive Along track spacing array, interpolated profiles and melt by Year
  
% Creates along track array query points for Interp1. Makes for all years
% but only use 1 array as interp1. Multiple years are so you can
% interpolate by which ever line you choose. Sample spacing is every 0.1m
% Switch middle value to change sample step
  
AT_data.query_array.P11 = (AT_data.AT_vel.P2011(1):0.1:...
  AT_data.AT_vel.P2011(end));
%AT_data.query_array.P14 = (AT_data.AT_vel.P2014(1):0.1:...
%AT_data.AT_vel.P2014(end));
%AT_data.query_array.P18 = (AT_data.AT_vel.P2018(1):0.1:...
%AT_data.AT_vel.P2018(end));
  
% Apply interpolation to each profile using selected query array
% 2011
AT_data.interp_data.P11 = interp1(AT_data.Btrack_End_Clip.P11, ...
  AT_data.elev_End_Clip.P2011, AT_data.query_array.P11);
 
% 2014
AT_data.interp_data.P14 = interp1(AT_data.Btrack_End_Clip.P14, ...
  AT_data.elev_End_Clip.P2014, AT_data.query_array.P11);
  
% 2018
AT_data.interp_data.P18 = interp1(AT_data.Btrack_End_Clip.P18, ...
  AT_data.elev_End_Clip.P2018, AT_data.query_array.P11);
  
% Calculate melt rates from interpolated profile pairings
% 2011-2014 melt (Vertical difference in Features)
AT_data.melt_rates.P11_P14 = AT_data.interp_data.P14 - ...
  AT_data.interp_data.P11; 
% 2014-2018 melt (Vertical difference in Features)
AT_data.melt_rates.P14_P18 = AT_data.interp_data.P18 - ...
  AT_data.interp_data.P14;
% 2011-2018 melt (Vertical difference in Features)
AT_data.melt_rates.P11_P18 = AT_data.interp_data.P18 - ...
  AT_data.interp_data.P11;
  
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
h3 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2011', 'original 2014', 'original 2018', ...
  'Location', 'southeast');
  
figure(2)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11);
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14);
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18);
h4 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
h5 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
h6 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
title('Test Plot 2 - Interpolated files same location as Originals');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2011', 'interpolated 2014', 'interpolated 2018', ...
  'old 2011','old 2014','old 2018', 'Location', 'southeast');
  
figure(3)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P11);
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P14);
h3 = plot(AT_data.query_array.P11/1e3, AT_data.interp_data.P18);
h4 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P18);
title('Test Plot 3 - Interpolation Profiles and Melt Profiles');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('interpolated 2011', 'interpolated 2014', 'interpolated 2018',...
  '2011-2018 melt', 'Location', 'southeast');
  
figure(4)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P14);
hold on
h2 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P14_P18);
h3 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P18);
title('Test Plot 4 - Melt Profiles for all years');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('2011-2014 melt', '2014-2018 melt', '2011-2018 melt', ...
  'Location', 'southeast');

% Subplot of melt rates for paired years
figure(5)
subplot(3,1,1);
plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P14, 'color','m');
title('Melt 2011-2014');
ylim([-100, 400]);

subplot(3,1,2);
plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P14_P18, 'color','c');
title('Melt 2014-2018');
ylim([-100, 400]);
ylabel('Elevation Change (m)');

subplot(3,1,3);
plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P18, 'color','g');
title('Melt 2011-2018');
ylim([-100, 400]);
xlabel('Along Track Distance (km)');

% melt rate figures
figure(6)
h1 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P18/7, 'color','g');
hold on 
h2 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P14/3, 'color','c');
h3 = plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P14_P18/4, 'color','m');
title('Melt Average annual melt rates based off vertical differences');
ylim([-20, 100]);
xlabel('Along Track Distance (km)');
ylabel('Vertical Melt (m)');
legend('2011-2018 annual average', '2011-2014 annual average', '2014-2018 annual average', 'Location', 'northeast');

figure(7)
plot(AT_data.query_array.P11/1e3, AT_data.melt_rates.P11_P18/7, 'color','g');
ylim([0, 50]);
title('Average annual difference in heigh between 2011-2018');
xlabel('Along Track Distance (km)');
ylabel('Vertical Melt (m)');
legend('2011-2018 annual average');
  %% NaN array padding section
  %makes NAN arrays and the transposes for reshapeing of along_track data 
  nan_array_11 = NaN(449, 1); %2011
  nan_array_11T = nan_array_11';
  
  nan_array_14 = NaN(261, 1); %2014
  nan_array_14T = nan_array_14';
  k= 449-261;
  nan_array_14B = NaN(k, 1); %2014 end nan array
  nan_array_14BT = nan_array_14B';
  
  %Concatenates to beginning/end of respective yearly profiles
  AT_data.elev_pad.P2011 = cat(1, nan_array_11, Bprof_2011);
  AT_data.elev_pad.P2014 = cat(1, nan_array_14, Bprof_2014, nan_array_14B);
  AT_data.elev_pad.P2018 = cat(1, Bprof_2018, nan_array_11);

  %Save along_track and vel values as individual variables (different sizes) in
  %structure
  
  %add yearly velocities to along_track master profile
  AT_data.CP.P2011 = AT_data.pass.P2018(1,2:end) + AT_data.vel.P2011;
  AT_data.CP.P2014 = AT_data.pass.P2018(1,2:end) + AT_data.vel.P2014;
  AT_data.CP.P2018 = AT_data.pass.P2018(1,2:end) + AT_data.vel.P2018;
  
  % Calculate melt by differencing annual profiles
  AT_data.melt.m11_18 = AT_data.elev_pad.P2018 - AT_data.elev_pad.P2011;
  AT_data.melt.m11_14 = AT_data.elev_pad.P2018 - AT_data.elev_pad.P2011;
  AT_data.melt.m14_18 = AT_data.elev_pad.P2018 - AT_data.elev_pad.P2014;
  
  %interpolate profiles to make them the same size as the along_track file
  AT_data.interp.P2011 = interp1(AT_data.elev_pad.P2011, AT_data.CP.P2011);
  AT_data.interp.P2014 = interp1(AT_data.elev_pad.P2014, AT_data.CP.P2014);
  AT_data.interp.P2018 = interp1(AT_data.elev_pad.P2018, AT_data.CP.P2018);
  
  %Test figures
  figure(1)
  h1 = plot(AT_data.CP.P2011/1e3, AT_data.elev.P2011);
  hold on
  h2 = plot(AT_data.CP.P2014/1e3, AT_data.elev.P2014);
  h3 = plot(AT_data.CP.P2018/1e3, AT_data.elev.P2018);
  xlabel = ('Along track distance (km)');
  ylabel = ('Elevation (m)');
  legend('2011', '2014', '2018','Location', 'southeast');
  
  figure(2)
  h1 = plot(AT_data.CP.P2011/1e3, AT_data.elev.P2011);
  hold on
  h2 = plot(AT_data.CP.P2014/1e3, AT_data.elev.P2014);
  h3 = plot(AT_data.CP.P2018/1e3, AT_data.elev.P2018);
  h4 = plot(AT_data.CP.P2018/1e3, AT_data.elev.P2018-AT_data.elev.P2011);
  xlabel = ('Along track distance (km)');
  ylabel = ('Elevation (m)');
  legend('new 2011', 'new 2014', 'new 2018','Location', 'southeast');
  
  figure(3)
  h1 = plot(AT_data.CP.P2011/1e3, AT_data.interp.P2011);
  hold on
  h2 = plot(AT_data.CP.P2014/1e3, AT_data.interp.P2014);
  h3 = plot(AT_data.CP.P2018/1e3, AT_data.interp.P2018);
  %h4 = plot(AT_data.CP.P2018/1e3, AT_data.interp.P2018-AT_data.interp.P2011);
  xlabel = ('Along track distance (km)');
  ylabel = ('Elevation (m)');
  legend('new 2011', 'new 2014', 'new 2018','Location', 'southeast');

%%
Crevasse_apex = localmax(AT_data.elevB.P2011);
Crevasse_base = localmin(AT_data.elevB.P2011);


%% Correlation Window approach
maxlag = 250; % minimum bin offset used by cameron
corr_window = 4* max_lag; 
% width of correlation window, real corr_window is +1

% make nan arrays the length of each profile for surface and bed
% surface
AT_data.window_S_11 = NaN(1,length(AT_data.elevS.P2011));
AT_data.window_S_14 = NaN(1,length(AT_data.elevS.P2014));
AT_data.window_S_18 = NaN(1,length(AT_data.elevS.P2018));
% Bed
AT_data.window_B_11 = NaN(1,length(AT_data.elevB.P2011));
AT_data.window_B_14 = NaN(1,length(AT_data.elevB.P2014));
AT_data.window_B_18 = NaN(1,length(AT_data.elevB.P2018));

% temporary correlation window file
temp_corr = NaN(length(AT_data.elevS.P2011),2*max_lag + 1);
temp_shift = [];
%%
%2011-2014 surface adjustment
for cidx = ((corr_window/2)+1):(length(AT_data.window_S_11) - (corr_window/2) - 1)
  % create tmp_signals 1 and 2
  tmp_signal_1 = AT_data.elevS.P2011(cidx - (corr_window/2):cidx+(corr_window/2));
  tmp_signal_2 = AT_data.elevS.P2014(cidx - (corr_window/2):cidx+(corr_window/2)); 
  
  [temp_corr(cidx,:),temp_lag] = xcorr(tmp_signal_1 - mean(tmp_signal_1),(tmp_signal_2 - mean(tmp_signal_2)), maxlag, 'coeff');
  vel_pick_window = (-10:10) + max_lag + 1 + AT_data.shift_data_14;
  
  
  if sum(temp_corr(cidx, max_lag + 1:end) == max(temp_corr(cidx, maxlag + 1:end))) ==1
    temp_shift = temp_lag(temp_corr(cidx, vel_pick_window) == max(temp_corr(cidx, vel_pick_window))) + vel_pick_window(1) -1;
    if temp_shift > vel_pick_window(1) - max_lag - 1
      AT_data.window_S_14(cidx) = temp_shift;
    end
  end
end
 %% Paired profiles 
  figure(1)
  h_plot = plot(alongtrack/1e3, Bprof_2011);
  grid on;
  hold on;
  h_plot2 = plot(alongtrack/1e3, Bprof_2014);
  legend({'2011 Profile', '2014 Profile'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Basal Reflector Elevation (m)');
  title(' 2011, 2014, 2018 Petermann Basal Reflectors');
  
  figure(2)
  h_plot = plot(alongtrack/1e3, Bprof_2014);
  grid on;
  hold on;
  h_plot2 = plot(alongtrack/1e3, Bprof_2018);
  legend({'2014 Profile', '2018 Profile'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title(' 2014-2018 Basal Melt');
  
  figure(3)
  h_plot = plot(alongtrack/1e3, Bprof_2011);
  grid on;
  hold on;
  h_plot2 = plot(alongtrack/1e3, Bprof_2018);
  legend({'2011 Profile', '2018 Profile'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title(' 2011-2018 Basal Melt');
  
  %% 
  figure(4)
  plot(alongtrack/1e3, Bmelt_14_11);
  grid on;
  xlabel('Along-track (km)');
  ylabel('\Delta height (M');
  title(' 2014-2011 Basal ');
  
  %%   
%   surf1 = cumtrapz(pass(1).layers(1).layer_elev);
%   base1 = cumtrapz(pass(1).layers(2).layer_elev);
%   surf2 = cumtrapz(pass(2).layers(1).layer_elev);
%   base2 = cumtrapz(pass(2).layers(2).layer_elev);
%   
%   prof2011 = trapz(pass(1).layers(1).layer_elev, pass(1).layers(2).layer_elev);
%   prof2014 = trapz(pass(2).layers(1).layer_elev, pass(2).layers(2).layer_elev); 
%   prof2018 = trapz(pass(3).layers(1).layer_elev, pass(3).layers(2).layer_elev);
  
  %% difference in Heigh of Profiles 
  fig = figure(5); 
  h1 = subplot(3,1,1);
  plot(alongtrack/1e3, Bmelt_14_11);
  grid on;
  title('Vertical Difference between 2014 and 2011 Basal Reflectors');
  
  h2 =  subplot(3,1,2);
  plot(alongtrack/1e3, Bmelt_18_14);
  grid on;
  title('Vertical Difference between 2018 and 2014 Basal Reflectors');

  h3 = subplot(3,1,3);
  plot(alongtrack/1e3, Bmelt_18_11);
  grid on;
  title('Vertical Difference between 2018 and 2011 Basal Reflectors');
 
  h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  
%   han = axes(fig, 'visible', 'off');
%   han.Title.Visible = 'on';
%   han.XLabel.Visible = 'on';
%   han.YLabel.Visible = 'on';
%   xlabel(han, 'Along-track (km)');
%   ylabel(han, 'Vertical change in Height (m)');
  
  %% Basal Reflectors of all 3 years WORKING SPOT
  figure(26)
  h_plot = plot(AT_data.AT_vel.P2011/1e3, Bprof_2011);
  grid on;
  hold on;
  h_plot2 = plot(AT_data.AT_vel.P2014/1e3, Bprof_2014);
  h_plot3 = plot(AT_data.AT_vel.P2018/1e3, Bprof_2018);
  legend({'2011 Profile', '2014 Profile', '2018 Profile'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title(' 2011, 2014, 2018 Petermann Basal Reflectors');
  
  figure(27)
  h_plot = plot(AT_data.AT_vel.P2011/1e3, Bprof_2011);
  grid on;
  hold on;
  %h_plot2 = plot(AT_data.AT_vel.P2014/1e3, Bprof_2014);
  h_plot3 = plot(AT_data.AT_vel.P2018/1e3, Bprof_2018);
  h_plot4 = plot(AT_data.AT_vel.P2011/1e3, Bprof_2018 - Bprof_2011);
  legend({'2011 Profile', '2014 Profile', '2018 Profile'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title(' 2011, 2014, 2018 Petermann Basal Reflectors');
  
  %%
  lat11 = pass(1).layers(2).lat;
  lat18 = pass(3).layers(2).lat;
  lon11 = pass(1).layers(2).lon;
  lon18 = pass(3).layers(2).lon;
  
  figure(31)
  h_plot1 = plot(lon11, lat11);
  hold on;
  h_plot2 = plot(lon18, lat18);
  grid on; 
  title('flight trajectories 2011/2018');
  legend({'2011', '2018'}, 'Location', 'northeast');
  xlabel('Longitude (degrees)');
  ylabel('Latitude (degrees)');
  
  
  
%   AT_data.CP = [];
%   AT_data.CP(1) = (pass(3).along_track(1,2:end) + pass(1).vel); 
%   AT_data.CP(2) = (pass(3).along_track(1,2:end) + pass(2).vel); 
%   AT_data.CP(3) = (pass(3).along_track(1,2:end) + pass(3).vel); 
  
%   for pass_idx = 1:length(pass)
%     AT_data.CP(pass_idx) = (pass(3).along_track(1,2:end)/1e3 + pass(pass_idx).vel);
%   end
    %% Interpolations for each profile
  for pass_idx = 1:length(pass)
    AT_data.interp_data(pass_idx) = interpr1();
  end
  
  
  

 
  
 %% 
    %Resample to same size
  AT_data.CP.P2011 = (AT_data.pass.P2011(1), length(AT_data.pass.P2011), length(AT_data.pass.P2018));
%%
  along_track

% AT_11 = linspace(along_track_11(1), along_track_11(end), along_track_18); 
%      AT_14 = linspace(along_track_14(1), along_track_14(end), along_track_18);
%      AT_18 = along_track_18;
    
  %Save Velocity correction data as individual variables (all same size)
  
  
  %Make new interpolated array using vel data for new along track
%   alongtrack_int = interp1(along);
  
  
  
  % Linspace data Testing for new along track interpolation
  lsTest = linspace(Bprof_2018(1),Bprof_2018(end),length(Bprof_2011));
  %alongtrack_int = interp1(lsTest, newprof11,  
  
  
  %alongtrack_int1 = interp1(AT_data.AT_vel.P2011, newprof11 ,AT_data.AT_vel.P2018);
  %alongtrack_int2 = interp1(AT_data.AT_vel.P2018, newprof18 ,AT_data.AT_vel.P2011);
  %alongtrack_int3 = interp1(newprof11, AT_data.AT_vel.P2018, newprof18);
  
  %alongtrack_new = cat(2, AT_data.AT_vel.P2011, nan_arrary_tran);
  %%
  figure(28)
  h1 = plot(tmp/1e3, newprof11);
  hold on;
  h2 = plot(tmp/1e3, newprof18);
  h3 = plot(tmp/1e3, (newprof18 - newprof11));
  h4 = plot(AT_data.AT_vel.P2011/1e3, Bprof_2011);
  h5 = plot(AT_data.AT_vel.P2018/1e3, Bprof_2018);
  legend({'new 2011', 'new 2018', 'melt difference', '2011 original', '2018 original'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title('Spatially correct, incorrect, and melt profiles 2011/2018');
  grid on;
  
  figure(29)
  h1 = plot(tmp/1e3, newprof11);
  hold on;
  h2 = plot(tmp/1e3, newprof18);
  h3 = plot(tmp/1e3, (newprof18 - newprof11));
  legend({'new 2011', 'new 2018', 'melt difference'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title('Spatially Incorrect 2011/2018 profiles');
  grid on;
  
  figure(30)
  h4 = plot(AT_data.AT_vel.P2011/1e3, Bprof_2011);
  hold on;
  h5 = plot(AT_data.AT_vel.P2018/1e3, Bprof_2018);
  legend({'original 2011', 'original 2018'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title('Spatially Correct 2011/2018 profiles');
  grid on;
  
  
  figure(32)
  h4 = plot(tmp/1e3, newprof11);
  hold on;
  h5 = plot(tmp/1e3, newprof18);
  legend({'new 2011', 'new 2018'}, 'Location', 'southeast');
  xlabel('Along-track (km)');
  ylabel('Cross-sectional Area(m^2)');
  title('Spatially Incorrect 2011/2018 profiles');
  grid on;
  %%
  figure(33)
  plot(tmp/1e3, (newprof18 - newprof11);
  legend({'2018-2011 melt'}, 'Location', 'southeast');
  title('Melt 2018-2011');
  xlabel('Cross-sectional area (km)');
  %% Cumulative loss
  
  fig = figure(8); 
  h1 = subplot(3,1,1);
  plot(alongtrack/1e3, area_11_14);
  grid on;
  title('Vertical Difference between 2014 and 2011 Basal Reflectors');
  
  h2 =  subplot(3,1,2);
  plot(alongtrack/1e3, area_14_18);
  grid on;
  title('Vertical Difference between 2018 and 2014 Basal Reflectors');

  h3 = subplot(3,1,3);
  plot(alongtrack/1e3, area_11_18);
  grid on;
  title('Vertical Difference between 2018 and 2011 Basal Reflectors');
 
  h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  %%
%    figure(3);
%   lay_idx = 1;
%   alongtrack = pass(baseline_master_idx).along_track + pass(pass_idx).vel;
%   h_plot(pass_idx) = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%   grid on;
%   leg_str{pass_idx} = num2str(year);
%   hold on;
%   for lay_idx = 2:length(pass(pass_idx).layers)
%     h_plot2 = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%     set(h_plot2,'Color',get(h_plot(pass_idx),'Color'));
%   end
%   xlabel('Along-track (km)');
%   ylabel('Elevation (m,WGS-84)');
  
  fig = figure(9); 
  h1 = subplot(3,1,1);
  plot(alongtrack/1e3, area_11_14);
  grid on;
  title('Vertical Difference between 2014 and 2011 Basal Reflectors');
  
  h2 =  subplot(3,1,2);
  plot(alongtrack/1e3, area_14_18);
  grid on;
  title('Vertical Difference between 2018 and 2014 Basal Reflectors');

  h3 = subplot(3,1,3);
  plot(alongtrack/1e3, area_11_18);
  grid on;
  title('Vertical Difference between 2018 and 2011 Basal Reflectors');
 
  h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  %% Velocity Corrected comparison Petermann Profile 1
  %Example Code from Velocity_coregistration.m
%     figure(3);
%   lay_idx = 1;
%   alongtrack = pass(baseline_master_idx).along_track + pass(pass_idx).vel;
%   h_plot(pass_idx) = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%   grid on;
%   leg_str{pass_idx} = num2str(year);
%   hold on;
%   for lay_idx = 2:length(pass(pass_idx).layers)
%     h_plot2 = plot(alongtrack/1e3, pass(pass_idx).layers(lay_idx).layer_elev);
%     set(h_plot2,'Color',get(h_plot(pass_idx),'Color'));
%   end
%   xlabel('Along-track (km)');
%   ylabel('Elevation (m,WGS-84)');
%   
% Naming each basal velocity pass profile and storing  
alongtrack_vel = pass(baseline_master_idx).along_track + pass(pass_idx).vel;

bV_prof_11 = (pass(1).layers(2).layer_elev); %+ pass(1).vel;
bV_prof_14 = (pass(2).layers(2).layer_elev); %+ pass(2).vel;
bV_prof_18 = (pass(3).layers(2).layer_elev); %+ pass(3).vel;

bv_area_11_14 = bV_prof_11 - bV_prof_14;
annual_11_14 = bv_area_11_14/3;

bv_area_14_18 = bV_prof_14 - bV_prof_18;
annual_14_18 = bv_area_14_18/4;

bv_area_11_18 = bV_prof_11 - bV_prof_18;
annual_11_18 = bv_area_11_18/7;

figure(10);
h1 = plot(alongtrack_vel/1e3, bv_area_11_14);
hold on
h2 = plot(alongtrack_vel/1e3, bv_area_14_18);
% h3 = plot(alongtrack_vel/1e3, bv_area_11_18);
legend('2011-14', '2014-18', 'Location', 'southeast');

h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');

figure(22); 
h_an = plot(alongtrack_vel/1e3, annual_11_18);
hold on
legend('annual melt 2011-18','Location', 'southeast');
title('Petermann Profile 1 2011-2018 annual vertical change');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');

  
figure(23);
h_an = plot(alongtrack_vel/1e3, annual_11_18); %annual height change
hold on
h_an2 = plot(alongtrack_vel/1e3, bv_area_11_18); %differece in total height
h_11 = plot(alongtrack_vel/1e3, Bprof_2011); %2011 bottom profile
h_18 = plot(alongtrack_vel/1e3, Bprof_2018); %2018 bottom profile
legend('annual melt 2011-18','total difference in 2011-18 height', ...
  '2011 botom profile', '2018 bottom profile', 'Location', 'southeast');
title('Petermann Profile 1');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  
%% Petermann Profile 2, 2013, 2014
  %along track distance and velocity correction from pass  
  alongtrack_vel = pass(baseline_master_idx).along_track + pass(pass_idx).vel;

  %Profile storing and naming
  prof2013 = trapz(pass(1).layers(1).layer_elev, pass(1).layers(2).layer_elev);
  prof2014 = trapz(pass(2).layers(1).layer_elev, pass(2).layers(2).layer_elev); 
  
  area_13_14 = prof2013 - prof2014;
  
  Bprof_2013 = (pass(1).layers(2).layer_elev);
  Bprof_2014 = (pass(2).layers(2).layer_elev);
  
  Sprof_2013 = (pass(1).layers(1).layer_elev);
  Sprof_2014 = (pass(2).layers(1).layer_elev);
  
  Bmelt_13_14 = Bprof_2013 - Bprof_2014; 
  Smelt_13_14 = Sprof_2013 - Sprof_2014;

  bV_prof_13 = (pass(1).layers(2).layer_elev); %+ pass(1).vel;
  bV_prof_14 = (pass(2).layers(2).layer_elev); %+ pass(2).vel;
  
  bv_area_13_14 = bV_prof_13 - bV_prof_14;
  annual_13_14 = bv_area_13_14/1;
  
figure(24); 
h_an = plot(alongtrack_vel/1e3, annual_13_14);
hold on
legend('annual melt 2013-14','Location', 'southeast');
title('Petermann Profile 1 2013-2014 annual vertical change');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  
figure(25);
h_an = plot(alongtrack_vel/1e3, annual_13_14); %annual height change
hold on
h_an2 = plot(alongtrack_vel/1e3, bv_area_13_14); %differece in total height
h_11 = plot(alongtrack_vel/1e3, Bprof_2013); %2011 bottom profile
h_18 = plot(alongtrack_vel/1e3, Bprof_2014); %2018 bottom profile
legend('annual melt 2013-14','total difference in 2013-14 height', ...
  '2013 botom profile', '2014 bottom profile', 'Location', 'southeast');
title('Petermann Profile 2');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
%% 79N Profile 1 2014, 2016, 2018
  %along track distance and velocity correction from pass  
  alongtrack_vel = pass(baseline_master_idx).along_track %+ pass(pass_idx).vel;

  %Profile storing and naming
  prof2014 = trapz(pass(1).layers(1).layer_elev, pass(1).layers(2).layer_elev);
  prof2016 = trapz(pass(2).layers(1).layer_elev, pass(2).layers(2).layer_elev);
  prof2018 = trapz(pass(3).layers(1).layer_elev, pass(3).layers(2).layer_elev); 
  
  area_14_16 = prof2014 - prof2016;
  area_16_18 = prof2016 - prof2018;
  area_16_18 = prof2014 - prof2018;
  
  Bprof_2014 = (pass(1).layers(2).layer_elev);
  Bprof_2016 = (pass(2).layers(2).layer_elev);
  Bprof_2018 = (pass(3).layers(2).layer_elev);
  
  Sprof_2014 = (pass(1).layers(1).layer_elev);
  Sprof_2016 = (pass(2).layers(1).layer_elev);
  Sprof_2018 = (pass(3).layers(1).layer_elev);
  
  Bmelt_14_18 = Bprof_2014 - Bprof_2018; 
  Smelt_14_18 = Sprof_2014 - Sprof_2018;

  bV_prof_14 = (pass(1).layers(2).layer_elev); %+ pass(1).vel;
  bV_prof_18 = (pass(3).layers(2).layer_elev); %+ pass(2).vel;
  
  bv_area_14_18 = bV_prof_14 - bV_prof_18;
  annual_14_18 = bv_area_14_18/4;
  
figure(26); 
h_an = plot(alongtrack_vel/1e3, annual_14_18);
hold on
legend('annual melt 2014-18','Location', 'southeast');
title('79N Profile 1 2014-2018 annual vertical change');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  
figure(27);
h_an = plot(alongtrack_vel/1e3, annual_14_18); %annual height change
hold on
h_an2 = plot(alongtrack_vel/1e3, bv_area_14_18); %differece in total height
h_11 = plot(alongtrack_vel/1e3, Bprof_2014); %2011 bottom profile
h_18 = plot(alongtrack_vel/1e3, Bprof_2018); %2014 bottom profile
legend('annual melt 2014-18','total difference in 2014-18 height', ...
  '2014 botom profile', '2018 bottom profile', 'Location', 'southeast');
title('79N Profile 1');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  
    alongtrack_vel = pass(baseline_master_idx).along_track %+ pass(pass_idx).vel;

  %Profile storing and naming
  prof2011 = trapz(pass(1).layers(1).layer_elev, pass(1).layers(2).layer_elev);
  prof2013 = trapz(pass(2).layers(1).layer_elev, pass(2).layers(2).layer_elev);
  prof2014 = trapz(pass(3).layers(1).layer_elev, pass(3).layers(2).layer_elev); 
  
  area_11_13 = prof2011 - prof2013;
  area_13_14 = prof2013 - prof2014;
  area_11_14 = prof2011 - prof2014;
  
  Bprof_2011 = (pass(1).layers(2).layer_elev);
  Bprof_2013 = (pass(2).layers(2).layer_elev);
  Bprof_2014 = (pass(3).layers(2).layer_elev);
  
  Sprof_2011 = (pass(1).layers(1).layer_elev);
  Sprof_2013 = (pass(2).layers(1).layer_elev);
  Sprof_2014 = (pass(3).layers(1).layer_elev);
  
  Bmelt_11_14 = Bprof_2011 - Bprof_2014; 
  Smelt_11_14 = Sprof_2011 - Sprof_2014;

  bV_prof_11 = (pass(1).layers(2).layer_elev); %+ pass(1).vel;
  bV_prof_14 = (pass(3).layers(2).layer_elev); %+ pass(2).vel;
  
  bv_area_11_14 = bV_prof_13 - bV_prof_14;
  annual_11_14 = bv_area_13_14/3;
  
figure(28); 
h_an = plot(alongtrack_vel/1e3, annual_11_14);
hold on
legend('annual melt 2011-14','Location', 'southeast');
title('Petermann Profile 4 2011-2014 annual vertical change');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
  
figure(29);
h_an = plot(alongtrack_vel/1e3, annual_11_14); %annual height change
hold on
h_an2 = plot(alongtrack_vel/1e3, bv_area_11_14); %differece in total height
h_11 = plot(alongtrack_vel/1e3, Bprof_2011); %2011 bottom profile
h_18 = plot(alongtrack_vel/1e3, Bprof_2014); %2014 bottom profile
legend('annual melt 2011-14','total difference in 2011-14 height', ...
  '2011 botom profile', '2014 bottom profile', 'Location', 'southeast');
title('Petermann Profile 4');
h = findobj(gcf, 'type', 'axes');
  set([h.YLabel],'string','\Delta Height(m)');
  set([h.XLabel], 'string','along track distance (km)');
