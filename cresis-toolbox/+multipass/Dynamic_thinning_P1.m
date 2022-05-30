%% DYNAMIC VERSION for autoloading all files in Folder
% Each File must have the same field and same beginning filename for
% processing in loop. If files do not match size and naming scheme this
% loop will fail.
file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\QGIS_processed_XY\');
myfiles = pwd; 
Velocity_CSV_1 = dir(fullfile(myfiles,'*.csv'));
    for i = 1:numel(Velocity_CSV_1)
        import_profile = fullfile(myfiles,Velocity_CSV_1(i).name);
        XY_data(i).data = readmatrix(import_profile);
        XY_data(i).X = XY_data(i).data(:,2).';
        XY_data(i).Y = XY_data(i).data(:,3).';
        XY_data(i).Lons = XY_data(i).data(:,4).';
        XY_data(i).Lats = XY_data(i).data(:,5).';
        XY_data(i).Surf = XY_data(i).data(:,6).';
        XY_data(i).Bed = XY_data(i).data(:,7).';
        XY_data(i).Thickness = XY_data(i).data(:,8).';
        XY_data(i).AT_dist = XY_data(i).data(:,9).';
        XY_data(i).Vel_mag = XY_data(i).data(:,10).';
        XY_data(i).div_12_13_2000m = XY_data(i).data(:,11).';
        XY_data(i).div_14_15_2000m = XY_data(i).data(:,12).';
        XY_data(i).div_15_16_2000m = XY_data(i).data(:,13).';
        XY_data(i).div_16_17_2000m = XY_data(i).data(:,14).';
        XY_data(i).div_17_18_2000m = XY_data(i).data(:,15).';
        XY_data(i).AT_vel_corrected = XY_data(1).AT_dist + XY_data(i).Vel_mag;
        XY_data(i).name_parts = strsplit(string(Velocity_CSV_1(i).name), '_');
        XY_data(i).year = str2double(XY_data(i).name_parts(2));
    end
    

%% 2 - CLIPPING EACH FIELD FOR ALIGNMENT
for i = 1:numel(XY_data) % define as # of passes or use length(fieldnames(XY_data)) on start 
   % Find first intersect point between arrays 
   XY_data(i).Btrack = XY_data(i).AT_vel_corrected(XY_data(i).AT_vel_corrected >= ...
       XY_data(1).AT_vel_corrected(1));
   % Locate Start Point of older array
   XY_data(i).find_AT_val = find(XY_data(i).AT_vel_corrected == ...
       XY_data(i).Btrack(1));
   % Clip Beginning & Ending
   XY_data(i).Beginning_clip = XY_data(i).AT_vel_corrected(XY_data(i).find_AT_val:end);
   XY_data(i).Ending_clip = XY_data(i).Beginning_clip(XY_data(i).Beginning_clip <= 4.698e+04);
   % Clip Beginning off Each Input Field
   XY_data(i).AT_Dist_Beg = XY_data(i).AT_dist(XY_data(i).find_AT_val:end);
   XY_data(i).Surf_Beg = XY_data(i).Surf(XY_data(i).find_AT_val:end); 
   XY_data(i).Bed_Beg = XY_data(i).Bed(XY_data(i).find_AT_val:end);   
   XY_data(i).Thickness_Beg = XY_data(i).Thickness(XY_data(i).find_AT_val:end);
   XY_data(i).div_12_13_2000m_Beg = XY_data(i).div_12_13_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_14_15_2000m_Beg = XY_data(i).div_14_15_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_15_16_2000m_Beg = XY_data(i).div_15_16_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_16_17_2000m_Beg = XY_data(i).div_16_17_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_17_18_2000m_Beg = XY_data(i).div_17_18_2000m(XY_data(i).find_AT_val:end);
   % Clip Endings off Each Input Field
   XY_data(i).AT_Dist_End = XY_data(i).AT_Dist_Beg(1:length(XY_data(i).Ending_clip));
   XY_data(i).Surf_End = XY_data(i).Surf_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Bed_End = XY_data(i).Bed_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Thickness_End = XY_data(i).Thickness_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_12_13_2000m_End = XY_data(i).div_12_13_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_14_15_2000m_End = XY_data(i).div_14_15_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_15_16_2000m_End = XY_data(i).div_15_16_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_16_17_2000m_End = XY_data(i).div_16_17_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_17_18_2000m_End = XY_data(i).div_17_18_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
end

%% 3 - INTERPOLATE SHORTENED DATASETS
% Define Query Array
XY_data(1).query_array = linspace(XY_data(1).Beginning_clip(1),...
    XY_data(1).Ending_clip(end), length(XY_data(1).Ending_clip));

% OR USE for points ever 14.76m
% XY_data(1).query_array = XY_data(1).Beginning_clip(1):14.76:XY_data(1).Ending_clip(end)

% Automated Version
% for i = 1:numel(XY_data)
%      XY_data(i).query_array = linspace(XY_data(i).Beginning_clip(1),...
%         XY_data(i).Ending_clip(end), length(XY_data(i).Ending_clip));
% end

% Inpterpolate Each Data Array
for i = 1:numel(XY_data)
   XY_data(i).Interp_AT_Dist = interp1(XY_data(i).Ending_clip,...
       XY_data(i).AT_Dist_End, XY_data(1).query_array, 'nearest','extrap');
   
   XY_data(i).Interp_Thickness = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Thickness_End, XY_data(1).query_array, 'nearest','extrap');
   
   XY_data(i).Interp_Surf = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Surf_End, XY_data(1).query_array, 'nearest','extrap');
   
   XY_data(i).Interp_Bed = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Bed_End, XY_data(1).query_array, 'nearest','extrap');   
   
   XY_data(i).Interp_Div_12_13 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_12_13_2000m_End, XY_data(1).query_array, 'nearest','extrap');
   
   XY_data(i).Interp_Div_14_15 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_14_15_2000m_End, XY_data(1).query_array, 'nearest','extrap');
   
   XY_data(i).Interp_Div_15_16 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_15_16_2000m_End, XY_data(1).query_array, 'nearest','extrap');
   
   XY_data(i).Interp_Div_16_17 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_16_17_2000m_End, XY_data(1).query_array, 'nearest','extrap'); 

   XY_data(i).Interp_Div_17_18 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_17_18_2000m_End, XY_data(1).query_array, 'nearest','extrap');   
end

%% 4 - Time Fraction Component, Tidal Inputs
% Time component fractions for percentages of years between acquisition
% dates of individual radar echograms (2014-2017). N_years variable denotes
% the total number of years between the 2 observed echogram periods, which
% is the sum of all the partial year components used to correct each
% divergence raster for the observed period. 
time_offset.P12_13 = (256/365);
time_offset.P13_14 = (132/365);
time_offset.P14_15 = (234/365);
time_offset.P15_16 = (365/365);
time_offset.P16_17 = (365/365);
time_offset.P17_18 = (90/365);
time_offset.N_years = time_offset.P14_15 + time_offset.P15_16...
   + time_offset.P16_17 + time_offset.P17_18; 

% Annual Tidal Correction Vertical Shift (m), Includes M2, S2, N2, K2, K1,
% O1, P1 and Q1 tidal constituents to predict tidal offset per image

% Surface Mass Balance Correction

%% 5 ICE DIVERGENCE 
% Save Smoothed Surfaces Into New Arrays
for i = 1:numel(XY_data)
    XY_data(i).Surf_smoothed = movmean(XY_data(i).Interp_Surf,10,'omitnan');
    XY_data(i).Bed_smoothed = movmean(XY_data(i).Interp_Bed,10,'omitnan');
end

% Mean Ice flux calculated from time-fraction corrected divergence tifs
for i = 1:numel(XY_data)
    XY_data(i).Mean_Div_flux = ((XY_data(i).Interp_Div_14_15*time_offset.P14_15)+...
        (XY_data(i).Interp_Div_15_16*time_offset.P15_16) + ...
        (XY_data(i).Interp_Div_16_17*time_offset.P16_17) + ...
        (XY_data(i).Interp_Div_17_18*time_offset.P17_18))/time_offset.N_years;
   
    XY_data(i).Mean_Div_flux_smooth = movmean(((XY_data(i).Interp_Div_14_15*time_offset.P14_15)+...
        (XY_data(i).Interp_Div_15_16*time_offset.P15_16) + ...
        (XY_data(i).Interp_Div_16_17*time_offset.P16_17) + ...
        (XY_data(i).Interp_Div_17_18*time_offset.P17_18))/time_offset.N_years, 10, 'omitnan');
end

% Calcuate Mean Thicknesses from paired years
for i = 1:numel(XY_data)
   XY_data(i).Mean_Thickness = (XY_data(i).Interp_Thickness + XY_data(i+1).Interp_Thickness)/2;
   XY_data(i).Mean_Thickness_smooth = movmean((XY_data(i+1).Interp_Thickness + XY_data(i+1).Interp_Thickness)/2, 10,'omitnan');
   
   XY_data(i).Bed_diff = XY_data(i+1).Interp_Bed - XY_data(i).Interp_Bed;
   XY_data(i).Bed_diff_smooth = movmean((XY_data(i+1).Interp_Bed - XY_data(i).Interp_Bed),10,'omitnan'); 
   
   XY_data(i).Surf_diff = XY_data(i+1).Interp_Surf - XY_data(i).Interp_Surf;
   XY_data(i).Surf_diff_smooth = movmean(XY_data(i+1).Interp_Surf,10,'omitnan') - ...
       movmean(XY_data(i+1).Interp_Surf,10,'omitnan');
   
   XY_data(i).Melt = XY_data(i).Bed_diff - (XY_data(i).Mean_Div_flux.*XY_data(i).Mean_Thickness);
   XY_data(i).Melt_smooth = XY_data(i).Bed_diff_smooth - (XY_data(i).Mean_Div_flux_smooth.*...
       XY_data(i).Mean_Thickness_smooth);
   
   XY_data(i).Dynamic_component = (XY_data(i).Mean_Div_flux.*XY_data(i).Mean_Thickness);
   XY_data(i).Dynamic_component_smooth = (XY_data(i).Mean_Div_flux_smooth.* XY_data(i).Mean_Thickness_smooth);
   if i > i(end)
        break
   end
end

%%
figure(1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt);
xlabel('Along Track Distance (km)');
ylabel('Basal Melt (m)');
title('Petermann Line 2 2014-2017 Basal Melt');

figure(2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Bed_diff, 'color', 'r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Bed_diff, 'color', 'b');
xlabel('Along Track Distance (km)');
ylabel('Bed Elevation difference (m)');
title('Bed Differences');
legend('2013-2014', '2014-2017');

figure(3)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Dynamic_component, 'color', 'b');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(2).Dynamic_component,50,'omitnan'),'color','r'); 
xlabel('Along Track Distance (km)');
ylabel('Bed Elevation difference (m)');
title('Dynamic thinning Component');
legend('2014-2017', '50 element moving mean');

figure(4)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt, 'color', 'r'); 
hold on 
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Bed_diff, 'color', 'b'); 
xlabel('Along Track Distance (km)');
ylabel('Melt/Bed Elevation Difference (m)');
title('Bed Elevation difference and Dynamic thinning corrected Melt Signal 2014-2017');
legend('Dynamic Thinning corrected', 'Delta Bed Elevation Diff');

figure(5)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf,'color','r');
hold on 
plot(XY_data(2).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf,'color','g');
plot(XY_data(3).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf,'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed,'color','r');
plot(XY_data(2).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed,'color','g');
plot(XY_data(3).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed,'color','b');
xlabel('Along Track Distance (km)');
ylabel('Basal Melt (m)');
title('Interpolated Ice Drafts Petermann line 2')
legend('2013','2014','2017');

figure(6)
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(2).Melt,10,'omitnan'),'color','r'); 
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(2).Melt,50,'omitnan'),'color','g'); 
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(2).Melt,100,'omitnan'),'color','b'); 
plot(XY_data(1).Interp_AT_Dist/1e3, sgolayfilt(XY_data(2).Melt,6,101),'color','k'); 
xlabel('Along Track Distance (km)');
ylabel('Melt Rate (m)');
legend('10 element mean','50 element mean', '100 element mean');
title('Moving Mean of 2014-2017 Basal Melt Signal');
%%
figure(10)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Dynamic_component,'color','r'); 
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Dynamic_component_smooth,'color','b'); 
legend('original','smoothed');
xlabel('along track distance (km)');
ylabel('Dynamic thinning component');

figure(11)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Bed_diff,'color','r'); 
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Bed_diff_smooth,'color','b'); 
legend('original','smoothed');
xlabel('along track distance (km)');
ylabel('Dynamic thinning component');
%% Smoothed Data
fig = figure(7); 
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_smooth/time_offset.N_years);
patch([5.2801 6.1048 6.1048 5.2801],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([6.525 7.29 7.29 6.525],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([9.0151 9.5101 9.5101 9.0151],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([11.25 12.21 12.21 11.25],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([17.1 17.895 17.895 17.1],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([24.255 28.14 28.14 24.255],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([29.475 30.225 30.225 29.475],[-500 -500 200 200],'r','FaceAlpha',0.25);
h1 = line([5.9399 5.9399],[-500 200]);
h2 = line([6.9151 6.9151],[-500 200]);
h3 = line([9.2851 9.2851],[-500 200]);
h4 = line([11.985 11.985],[-500 200]);
h5 = line([17.505 17.505],[-500 200]);
h6 = line([25.89 25.89],[-500 200]);
h7 = line([26.865 26.865],[-500 200]);
h8 = line([29.85 29.85],[-500 200]);
set([h1 h2 h3 h4 h5 h6 h7 h8], 'Color','k','LineStyle','--');
ylabel('Annual Basal Melt (m/yr)');
ylim([-10 60]);
legend('2014-2017 mean annual melt');

subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(2).Interp_Surf, 10, 'omitnan'),'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(3).Interp_Surf, 10, 'omitnan'),'color','k');
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(2).Interp_Bed, 10, 'omitnan'),'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, movmean(XY_data(3).Interp_Bed, 10, 'omitnan'),'color','k');
patch([5.2801 6.1048 6.1048 5.2801],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([6.525 7.29 7.29 6.525],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([9.0151 9.5101 9.5101 9.0151],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([11.25 12.21 12.21 11.25],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([17.1 17.895 17.895 17.1],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([24.255 28.14 28.14 24.255],[-500 -500 200 200],'r','FaceAlpha',0.25);
patch([29.475 30.225 30.225 29.475],[-500 -500 200 200],'r','FaceAlpha',0.25);
h1 = line([5.9399 5.9399],[-500 200]);
h2 = line([6.9151 6.9151],[-500 200]);
h3 = line([9.2851 9.2851],[-500 200]);
h4 = line([11.985 11.985],[-500 200]);
h5 = line([17.505 17.505],[-500 200]);
h6 = line([25.89 25.89],[-500 200]);
h7 = line([26.865 26.865],[-500 200]);
h8 = line([29.85 29.85],[-500 200]);
set([h1 h2 h3 h4 h5 h6 h7 h8], 'Color','k','LineStyle','--');
ylabel('WGS84 Elevation (m)');
ylim([-500 200]);
title('Interpolated Ice Drafts Petermann line 2')
legend('2014','2017','Location','Best');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Along Track Distance (km)');
title(han,'Smoothed Basal Melt & Along Track Profile');

%% Original Data
fig = figure(8); 
subplot(2,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt/time_offset.N_years);
xlabel('Along Track Distance (km)');
ylabel('Annual Basal Melt (m/yr)');
ylim([-10 60]);

subplot(2,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf,'color','g');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf,'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed,'color','g');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed,'color','b');
ylabel('WGS84 Elevation (m)');
title('Interpolated Ice Drafts Petermann line 2')
legend('2014','2017','Location','Best');


han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
xlabel(han,'Along Track Distance (km)');
title(han,'Original Basal Melt & Along Track Profile');
