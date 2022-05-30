% Dynamic_Thinning_Analysis
% DESCRIPTION: Computes both the Langragian full-shelf dynamic thinning 
% rate and the Eulierian local dynamic thinning rate of each individual
% flight line or crevasse. Langragian full-shelf thinning rates are
% computed using the sum of DU/dx and DV/dy for 2 separate years, that are
% then averaged. Each velocity/distance difference is computer first for
% each year, then they are averaged after finding each respective years
% DU/dx and DV/dy components. 
% 
% INPUTS: pass.along_track, XY velocity components (from QGIS), profile 
% thickness, time difference (year or days), crevasse position & heights
%
% OUTPUT: Excel Spreadsheet with new column for Dynamic thinning, plots of
% dyanmic thinning rates, and plots for local crevasse thinning
%
% AUTHOR: Cody Barnett
%
%% 1) Load in new XY velocity appended data derived from QGIS raster 
% sampling process. Prior data was loaded from Export_multipass_struct
% script where data was uploaded into QGIS and rasters were sampled
    
% Change file path to location of XY data
file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\QGIS_processed_XY\');
myfiles = pwd; 
    
% change string to load desired years
Velocity_CSV = dir(fullfile(myfiles,'P2_divergence_y17_resamples*')); 
for i = 1:numel(Velocity_CSV)
    Import_Profile = fullfile(myfiles,Velocity_CSV(i).name);
    XY_data(i).data = readmatrix(Import_Profile);
    XY_data(i).X = XY_data(i).data(:,2).';
    XY_data(i).Y = XY_data(i).data(:,3).';
    XY_data(i).Lons = XY_data(i).data(:,4).';
    XY_data(i).Lats = XY_data(i).data(:,5).';
    XY_data(i).Surf = XY_data(i).data(:,6).';
    XY_data(i).Bed = XY_data(i).data(:,7).';
    XY_data(i).Thickness = XY_data(i).data(:,8).';
    XY_data(i).AT_dist = XY_data(i).data(:,9).';
    XY_data(i).Vel_mag = XY_data(i).data(:,10).';
    XY_data(i).div_200m = XY_data(i).data(:,11).';
    XY_data(i).div_500m = XY_data(i).data(:,12).';
    XY_data(i).div_1000m = XY_data(i).data(:,13).';
    XY_data(i).div_1500m = XY_data(i).data(:,14).';
    XY_data(i).div_2000m = XY_data(i).data(:,15).';
%     Old Format     
%     XY_data(i).XVel_14_15 = XY_data(i).data(:,11).';
%     XY_data(i).XVel_15_16 = XY_data(i).data(:,12).';
%     XY_data(i).XVel_16_17 = XY_data(i).data(:,13).';
%     XY_data(i).YVel_14_15 = XY_data(i).data(:,14).';
%     XY_data(i).YVel_15_16 = XY_data(i).data(:,15).';
%     XY_data(i).YVel_16_17 = XY_data(i).data(:,16).';        
end
%%
for i = 1:numel(XY_data)
   XY_data(i).dynamic_thinning_200m = XY_data(i).div_200m.*XY_data(i).Thickness;
   XY_data(i).dynamic_thinning_500m = XY_data(i).div_500m.*XY_data(i).Thickness;
   XY_data(i).dynamic_thinning_1000m = XY_data(i).div_1000m.*XY_data(i).Thickness;
   XY_data(i).dynamic_thinning_1500m = XY_data(i).div_1500m.*XY_data(i).Thickness;
   XY_data(i).dynamic_thinning_2000m = XY_data(i).div_2000m.*XY_data(i).Thickness;
end
figure(5)
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_200m);
hold on
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_500m);
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_1000m);
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_1500m);
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_2000m);
xlabel('along track distance (km)');
ylabel('dynamic thinning (m)');
title('Dynamic thinning Signals 2017/2018');
legend('200m','500m','1000m','1500m','2000m');

fig = figure(6);
subplot(5,1,1);
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_200m);
legend('200m','Location','southeast');

subplot(5,1,2)
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_500m);
legend('500m','Location','southeast');

subplot(5,1,3)
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_1000m);
legend('1000m','Location','southeast');

subplot(5,1,4)
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_1500m);
legend('1500m','Location','southeast');
ylim([-50 10]);

subplot(5,1,5)
plot(XY_data.AT_dist/1e3, XY_data.dynamic_thinning_2000m);
legend('2000m','Location','southeast');
ylim([-30 10]);

han = axes(fig, 'visible','off');
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Dynamic thinning');
xlabel(han,'along track distance (km)');
title(han,'Individual dynamic Thinning signals - 2017/2018');
%% Compute Ice Flux and Dynamic thinning
% Flux multipled by thickness
for i = 1:numel(XY_data)
%    XY_data(i).dynamic_thinning_14 = XY_data(i).div_14_15.*XY_data(i).Thickness;
%    XY_data(i).dynamic_thinning_15 = XY_data(i).div_15_16.*XY_data(i).Thickness;
   XY_data(i).dynamic_thinning_16 = XY_data(i).div_16_17.*XY_data(i).Thickness;
end
mean_dynamic = mean(XY_data(i).dynamic_thinning_16);
smoothed_rate = smooth(XY_data(1).div_16_17,0.1,'rloess');
smoothed_rate_t = smoothed_rate.';
smooth_thickness = smoothed_rate_t.*XY_data(1).Thickness;
%%
figure(4)
plot(XY_data(1).AT_dist/1000, XY_data(1).dynamic_thinning_16);
hold on 
plot(XY_data(1).AT_dist/1000, zeros(3450));
plot(XY_data(1).AT_dist/1000, smooth_thickness,'Color','k');
xlabel('along track distance (km)');
ylabel('dynamic thinning (m/yr)');
title('500m Dynamic thinning signal');


figure(5)
plot(XY_data(1).AT_dist/1000, XY_data(1).div_16_17);
hold on 
plot(XY_data(1).AT_dist/1000, zeros(3450));
plot(XY_data(1).AT_dist/1000, smoothed_rate,'Color','k');
xlabel('along track distance (km)');
ylabel('ice divergence');
legend('500m Signal', 'zero line', 'smoothed 500m');
title('500m Divergence and Smoothed signal');
%%
% Change file path to location of XY data
file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\QGIS_processed_XY\');
myfiles = pwd; 
    
% change string to load desired years
Velocity_CSV = dir(fullfile(myfiles,'N79*')); 
for i = 1:numel(Velocity_CSV)
    Import_Profile = fullfile(myfiles,Velocity_CSV(i).name);
    XY_data(i).data = readmatrix(Import_Profile);
    XY_data(i).Lat = XY_data(i).data(:,1).';
    XY_data(i).Lon = XY_data(i).data(:,2).';
    XY_data(i).time = XY_data(i).data(:,3).';
    XY_data(i).thickness = XY_data(i).data(:,4).';
    XY_data(i).elev = XY_data(i).data(:,5).';
    XY_data(i).frame = XY_data(i).data(:,6).';
    XY_data(i).surf = XY_data(i).data(:,7).';
    XY_data(i).bed = XY_data(i).data(:,8).';
    XY_data(i).quality = XY_data(i).data(:,9).';
    XY_data(i).div_16_17 = XY_data(i).data(:,10).';
%     Old Format     
%     XY_data(i).XVel_14_15 = XY_data(i).data(:,11).';
%     XY_data(i).XVel_15_16 = XY_data(i).data(:,12).';
%     XY_data(i).XVel_16_17 = XY_data(i).data(:,13).';
%     XY_data(i).YVel_14_15 = XY_data(i).data(:,14).';
%     XY_data(i).YVel_15_16 = XY_data(i).data(:,15).';
%     XY_data(i).YVel_16_17 = XY_data(i).data(:,16).';    
end
figure(10)
plot3(XY_data.Lon, XY_data.Lat, XY_data.div_16_17)%*XY_data.thickness);
hold on 
% plot3(XY_data.Lon, XY_data.Lat, XY_data.surf);
% plot3(XY_data.Lon, XY_data.Lat, XY_data.bed)
xlabel('longitude');
ylabel('Latitude');
zlabel('divergence');
%%
% Mean of 2014 and 2017
mean_thinning = (XY_data(1).dynamic_thinning_14 + XY_data(2).dynamic_thinning_16)/2;
array_size = 3450;
xfit = linspace(min(XY_data(1).AT_dist/1000), max(XY_data(1).AT_dist/1000), array_size);
yfit = polyval(best_fit, xfit);

best_fit = polyfit(XY_data(1).AT_dist/1000, mean_thinning, 2);
best_fit_vals = polyval(best_fit, XY_data(1).AT_dist/1000);

figure(1)
plot(XY_data(1).AT_dist/1000, mean_thinning);
hold on
plot(xfit, yfit);
xlabel('Along Track distance (km)');
ylabel('Thinning Rate (m)');
legend show;
title('Dynamic Thinning Petermann Line 2, 2014-2017 Average');

figure(2)
plot(XY_data(1).AT_dist/1000, XY_data(1).div_14_15);
hold on
plot(XY_data(2).AT_dist/1000, XY_data(2).div_16_17);
xlabel('Along Track distance (km)');
ylabel('del U (m)');
legend('2014','2017');
title('Dynamic Thinning Petermann Line 2, 2014-2017 Average');

%% Compute compounded component velocity values for difference between years (OLD FORMAT)
for i = 1:numel(XY_data)
    XY_data(i).X_compound = XY_data(i).XVel_14_15 + XY_data(i).XVel_15_16 + XY_data(i).XVel_16_17;
    XY_data(i).Y_compound = XY_data(i).YVel_14_15 + XY_data(i).YVel_15_16 + XY_data(i).YVel_16_17;
end

% Find the difference between sequential points X, Y, X_vel, Y_vel fields
for i = 1:numel(XY_data)
    % Along Track and thickness Distance difference for plotting
    XY_data(i).AT_dif = XY_data(i).AT_dist(2:end) - XY_data(i).AT_dist(1:end-1);
    XY_data(i).Thickness_dif = XY_data(i).Thickness(2:end) - XY_data(i).Thickness(1:end-1);
    % X & Y position differences
    XY_data(i).X_dist = XY_data(i).X(2:end) - XY_data(i).X(1:end-1);
    XY_data(i).Y_dist = XY_data(i).Y(2:end) - XY_data(i).Y(1:end-1);
    % X-Vel annual difference
    XY_data(i).X_vel_dif_14 = XY_data(i).XVel_14_15(2:end) - XY_data(i).XVel_14_15(1:end-1);
    XY_data(i).X_vel_dif_15 = XY_data(i).XVel_15_16(2:end) - XY_data(i).XVel_15_16(1:end-1);
    XY_data(i).X_vel_dif_16 = XY_data(i).XVel_16_17(2:end) - XY_data(i).XVel_16_17(1:end-1);
    % Y-Vel annual difference
    XY_data(i).Y_vel_dif_14 = XY_data(i).YVel_14_15(2:end) - XY_data(i).YVel_14_15(1:end-1);
    XY_data(i).Y_vel_dif_15 = XY_data(i).YVel_15_16(2:end) - XY_data(i).YVel_15_16(1:end-1);
    XY_data(i).Y_vel_dif_16 = XY_data(i).YVel_16_17(2:end) - XY_data(i).YVel_16_17(1:end-1);
end

% Compute ice flux for 2014
for i = 1:numel(XY_data) 
   % Compute DU_DX and DV_DY 
   XY_data(i).DU_DX_14 = XY_data(i).X_vel_dif_14 ./ XY_data(i).X_dist;
   XY_data(i).DV_DY_14 = XY_data(i).Y_vel_dif_14 ./ XY_data(i).Y_dist; 
   % Add flux components to get total divergence
   XY_data(i).Divergence_14 = XY_data(i).DU_DX_14 + XY_data(i).DV_DY_14;
   % Multiply by thickness to get total ice flux
   XY_data(i).ice_flux_14 = XY_data(i).Divergence_14.*XY_data(i).Thickness_dif;
end

% Compute ice flux for 2017
for i = 1:numel(XY_data) 
   % Compute DU_DX and DV_DY 
   XY_data(i).DU_DX_16 = XY_data(i).X_vel_dif_16./XY_data(i).X_dist;
   XY_data(i).DV_DY_16 = XY_data(i).Y_vel_dif_16./XY_data(i).Y_dist; 
   % Add flux components to get total divergence
   XY_data(i).Divergence_16 = XY_data(i).DU_DX_16 + XY_data(i).DV_DY_16;
   % Multiply by thickness to get total ice flux
   XY_data(i).ice_flux_16 = XY_data(i).Divergence_16.*XY_data(i).Thickness_dif;
end

% mean dynamic thinning
mean_dynamic_thinning = (XY_data(1).ice_flux_14 + XY_data(2).ice_flux_16)/2;

% Plot dynamic thinning component
figure(1)
plot(XY_data(1).AT_dist(1:end-1)/1000, mean_dynamic_thinning)
xlabel('along track distance (km)') ;
ylabel('dyanmic thinning offset (m)');
title('Dynamic thinning 2014-2017 Petermann Line 2'); 
legend show;

%% Unqiue X coordinate ID approach (FIX LATER)
for i = 1:numel(XY_data)
    % Find unique x values, their index and index within original array
    [XY_data(i).unique_x, XY_data(i).x_orig_idx, XY_data(i).x_new_idx] = unique(XY_data(i).X, 'first', 'legacy');
    % Flip arrays so they are in correct order
    XY_data(i).x_orig_idx = flip(XY_data(i).x_orig_idx);
    XY_data(i).x_new_idx = flip(XY_data(i).x_new_idx);
    % mask out y-array so that it has values only based on x index
    XY_data(i).unique_y_mask = (XY_data(i).Y(XY_data(i).x_new_idx));
    %XY_data(i).unique_y = 
end


%% 2) Load in along_track distance, lidar surface, bed, and velocity 
% magnitude and save to individual arrays within the pass structure

file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\CSV_export_files\');
myfiles = pwd;

% Calculate midpoint between each array
for j = 1:numel(pass)
    % pull along track and subtract second element from first
    along_track = pass(j).along_track;
    pass(j).point_dist = along_track(2:end) - along_track(1:end-1);
    
    % pull velocity magnitude and subtract second element from first (km)
    vel_magnitude = pass(j).vel;
    pass(j).vel_mag_difference = (vel_magnitude(2:end)/1000) - (vel_magnitude(1:end-1)/1000);
    
    
end