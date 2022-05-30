%% DYNAMIC VERSION for autoloading all files in Folder - Petermann 1
% Each File must have the same field and same beginning filename for
% processing in loop. If files do not match size and naming scheme this
% loop will fail.
file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Extended_Passes\');
myfiles = pwd; 
Velocity_CSV_1 = dir(fullfile(myfiles,'P1M_div*.csv'));
    for i = 1:numel(Velocity_CSV_1)
        import_profile = fullfile(myfiles,Velocity_CSV_1(i).name);
        XY_data(i).data = readmatrix(import_profile);
        XY_data(i).fid = XY_data(i).data(1:4436,1).';
        XY_data(i).X = XY_data(i).data(1:4436,2).';
        XY_data(i).Y = XY_data(i).data(1:4436,3).';
        XY_data(i).Lons = XY_data(i).data(1:4436,4).';
        XY_data(i).Lats = XY_data(i).data(1:4436,5).';
        XY_data(i).Surf = XY_data(i).data(1:4436,6).';
        XY_data(i).Bed = XY_data(i).data(1:4436,7).';
        XY_data(i).Thickness = XY_data(i).data(1:4436,8).';
        XY_data(i).AT_dist = XY_data(i).data(1:4436,9).';
        XY_data(i).Vel_mag = XY_data(i).data(1:4436,10).';
        XY_data(i).div_09_10_2000m = XY_data(i).data(1:4436,11).';
        XY_data(i).div_12_13_2000m = XY_data(i).data(1:4436,12).';
        XY_data(i).div_14_15_2000m = XY_data(i).data(1:4436,13).';
        XY_data(i).div_15_16_2000m = XY_data(i).data(1:4436,14).';
        XY_data(i).div_16_17_2000m = XY_data(i).data(1:4436,15).';
        XY_data(i).div_17_18_2000m = XY_data(i).data(1:4436,16).';
        XY_data(i).SMB_10 = XY_data(i).data(1:4436,17).';
        XY_data(i).SMB_11 = XY_data(i).data(1:4436,18).';
        XY_data(i).SMB_12 = XY_data(i).data(1:4436,19).';
        XY_data(i).SMB_13 = XY_data(i).data(1:4436,20).';
        XY_data(i).SMB_14 = XY_data(i).data(1:4436,21).';
        XY_data(i).SMB_15 = XY_data(i).data(1:4436,22).';        
        XY_data(i).SMB_16 = XY_data(i).data(1:4436,23).';
        XY_data(i).SMB_17 = XY_data(i).data(1:4436,24).';
        XY_data(i).SMB_18 = XY_data(i).data(1:4436,25).';
        XY_data(i).SMB_19 = XY_data(i).data(1:4436,26).';       
        XY_data(i).AT_vel_corrected = XY_data(1).AT_dist + XY_data(i).Vel_mag;
        XY_data(i).name_parts = strsplit(string(Velocity_CSV_1(i).name), '_');
        XY_data(i).name_parts2 = strsplit(string(XY_data(i).name_parts(4)),'.');
        XY_data(i).year = str2double(XY_data(i).name_parts2(1));
    end
XY_data = rmfield(XY_data,'name_parts2');
%%
figure(12)
plot(XY_data(1).AT_dist/1e3, XY_data(1).Surf)
hold on
plot(XY_data(1).AT_dist/1e3, XY_data(1).Bed)
figure(13)
plot(XY_data(2).AT_dist/1e3, XY_data(2).Surf)
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).Bed)
figure(14)
plot(XY_data(3).AT_dist/1e3, XY_data(3).Surf)
hold on
plot(XY_data(3).AT_dist/1e3, XY_data(3).Bed)
figure(15)
plot(XY_data(4).AT_dist/1e3, XY_data(4).Surf)
hold on
plot(XY_data(4).AT_dist/1e3, XY_data(4).Bed)
figure(16)
plot(XY_data(5).AT_dist/1e3, XY_data(5).Surf)
hold on
plot(XY_data(5).AT_dist/1e3, XY_data(5).Bed)
figure(17)
plot(XY_data(6).AT_dist/1e3, XY_data(6).Surf)
hold on
plot(XY_data(6).AT_dist/1e3, XY_data(6).Bed)
%% TIDAL CORRECTIONS
% define mean tides over observation period. If values are negative
% subtract from first value instead of adding. If mean tidal value is
% positive, the value is subtracted from total surface elevation. If tidal
% value is negative, the value is added in order to set tidal applitude to
% 0 between all echograms.

% tidal Amplitudes
XY_data(1).mean_tidal_c = (0.3302 + 0.3390 + 0.3452)/3;
XY_data(2).mean_tidal_c = (-0.1859 + (-0.5356))/2;
XY_data(3).mean_tidal_c = (1.7635 + 1.7762)/2;
XY_data(4).mean_tidal_c = (0.4308 + 0.4128)/2;
XY_data(5).mean_tidal_c = (0.1733 + 0.1709)/2;
XY_data(6).mean_tidal_c = (1.887 + 1.9808)/2;

% Extract Ocean Surface with tidal correction, set to 0 tidal level
XY_data(1).ocean_surface = XY_data(1).Surf(3891:4225) - XY_data(1).mean_tidal_c;
XY_data(2).ocean_surface = XY_data(2).Surf(3696:4030) + XY_data(2).mean_tidal_c;
XY_data(3).ocean_surface = XY_data(3).Surf(3770:4104) - XY_data(3).mean_tidal_c;
XY_data(4).ocean_surface = XY_data(4).Surf(3930:4264) - XY_data(4).mean_tidal_c;
XY_data(5).ocean_surface = XY_data(5).Surf(4045:4339) - XY_data(5).mean_tidal_c;
XY_data(6).ocean_surface = XY_data(6).Surf(4102:4436) - XY_data(6).mean_tidal_c;

% Extract AT array for each ocean surface, set to 0 tidal level
XY_data(1).AT_ocean = XY_data(1).AT_dist(3891:4225); 
XY_data(2).AT_ocean = XY_data(2).AT_dist(3696:4030);
XY_data(3).AT_ocean = XY_data(3).AT_dist(3770:4104);
XY_data(4).AT_ocean = XY_data(4).AT_dist(3930:4264);
XY_data(5).AT_ocean = XY_data(5).AT_dist(4045:4339);
XY_data(6).AT_ocean = XY_data(6).AT_dist(4102:4436);

% 1) Calculate distance between terminus and edge of echogram
% 2) set start distance to 0 for each along track array
% 3) Obtain mask of all along track values less than 5km 
% 4) Obtain values of along track that is less than the 5km
% 5) Clip ocean surface to the same size as the along track clipped array
% 6) Calculate mean elevation offset of ocean surface
for i = 1:numel(XY_data)
    XY_data(i).terminus_dist = (XY_data(i).AT_ocean(end) - XY_data(i).AT_ocean(1))/1000;
    XY_data(i).AT_ocean_0 = (XY_data(i).AT_ocean - XY_data(i).AT_ocean(1))/1000;
    XY_data(i).AT_ocean_5km_mask = XY_data(i).AT_ocean_0 < 5;
    XY_data(i).AT_ocean_5km = XY_data(i).AT_ocean_0(XY_data(i).AT_ocean_5km_mask == 1);
    XY_data(i).ocean_surface_c = XY_data(i).ocean_surface(1:length(XY_data(i).AT_ocean_5km));
    XY_data(i).mean_ocean_surf_5km = nanmean(XY_data(i).ocean_surface_c);
end

% Remove mean sea surface offset and tidal offset from overall surface/bed 
% signals, setting mean elevation offset value to 0 per profile
XY_data(1).surf_tidal_elevation_align = XY_data(1).Surf - XY_data(1).mean_ocean_surf_5km/2 + XY_data(1).mean_tidal_c;
XY_data(2).surf_tidal_elevation_align = XY_data(2).Surf - XY_data(2).mean_ocean_surf_5km - XY_data(2).mean_tidal_c;
XY_data(3).surf_tidal_elevation_align = XY_data(3).Surf - XY_data(3).mean_ocean_surf_5km + XY_data(3).mean_tidal_c;
XY_data(4).surf_tidal_elevation_align = XY_data(4).Surf - XY_data(4).mean_ocean_surf_5km - XY_data(4).mean_tidal_c;
XY_data(5).surf_tidal_elevation_align = XY_data(5).Surf - XY_data(5).mean_ocean_surf_5km - XY_data(5).mean_tidal_c;
XY_data(6).surf_tidal_elevation_align = XY_data(6).Surf - XY_data(6).mean_ocean_surf_5km - XY_data(6).mean_tidal_c;

XY_data(1).bed_tidal_elevation_align = XY_data(1).Bed - XY_data(1).mean_ocean_surf_5km/2 + XY_data(1).mean_tidal_c;
XY_data(2).bed_tidal_elevation_align = XY_data(2).Bed - XY_data(2).mean_ocean_surf_5km - XY_data(2).mean_tidal_c;
XY_data(3).bed_tidal_elevation_align = XY_data(3).Bed - XY_data(3).mean_ocean_surf_5km + XY_data(3).mean_tidal_c;
XY_data(4).bed_tidal_elevation_align = XY_data(4).Bed - XY_data(4).mean_ocean_surf_5km - XY_data(4).mean_tidal_c;
XY_data(5).bed_tidal_elevation_align = XY_data(5).Bed - XY_data(5).mean_ocean_surf_5km - XY_data(5).mean_tidal_c;
XY_data(6).bed_tidal_elevation_align = XY_data(6).Bed - XY_data(6).mean_ocean_surf_5km - XY_data(6).mean_tidal_c;

% surface
fig = figure(113);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).surf_tidal_elevation_align);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).surf_tidal_elevation_align);
plot(XY_data(3).AT_dist/1e3, XY_data(3).surf_tidal_elevation_align);
plot(XY_data(4).AT_dist/1e3, XY_data(4).surf_tidal_elevation_align);
plot(XY_data(5).AT_dist/1e3, XY_data(5).surf_tidal_elevation_align);
plot(XY_data(6).AT_dist/1e3, XY_data(6).surf_tidal_elevation_align);
title('2011,2014,2015,2017,2018,2019 Unaligned surfaces');
ax1.YLim = [-10 110];

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).surf_tidal_elevation_align);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).surf_tidal_elevation_align);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).surf_tidal_elevation_align);
plot(XY_data(4).AT_vel_corrected/1e3, XY_data(4).surf_tidal_elevation_align);
plot(XY_data(5).AT_vel_corrected/1e3, XY_data(5).surf_tidal_elevation_align);
plot(XY_data(6).AT_vel_corrected/1e3, XY_data(6).surf_tidal_elevation_align);
title('2011,2014,2015,2017,2018,2019 Aligned surfaces');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Elevation WGS (m)');
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'xy')

% bed
fig = figure(114);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).bed_tidal_elevation_align);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).bed_tidal_elevation_align);
plot(XY_data(3).AT_dist/1e3, XY_data(3).bed_tidal_elevation_align);
plot(XY_data(4).AT_dist/1e3, XY_data(4).bed_tidal_elevation_align);
plot(XY_data(5).AT_dist/1e3, XY_data(5).bed_tidal_elevation_align);
plot(XY_data(6).AT_dist/1e3, XY_data(6).bed_tidal_elevation_align);
title('2011,2014,2015,2017,2018,2019 Unaligned surfaces');
%ax1.YLim = [-10 110];

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).bed_tidal_elevation_align);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).bed_tidal_elevation_align);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).bed_tidal_elevation_align);
plot(XY_data(4).AT_vel_corrected/1e3, XY_data(4).bed_tidal_elevation_align);
plot(XY_data(5).AT_vel_corrected/1e3, XY_data(5).bed_tidal_elevation_align);
plot(XY_data(6).AT_vel_corrected/1e3, XY_data(6).bed_tidal_elevation_align);
title('2011,2014,2015,2017,2018,2019 Aligned surfaces');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Elevation WGS (m)');
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'xy')

%% Polynomial fitted curves for surface and bed
XY_data(1).surf_tidal_align_polyfit = polyfit(XY_data(1).AT_dist, XY_data(1).surf_tidal_elevation_align, 50);
XY_data(2).surf_tidal_align_polyfit = polyfit(XY_data(2).AT_dist, XY_data(2).surf_tidal_elevation_align, 50);
XY_data(3).surf_tidal_align_polyfit = polyfit(XY_data(3).AT_dist, XY_data(3).surf_tidal_elevation_align, 50);
XY_data(4).surf_tidal_align_polyfit = polyfit(XY_data(4).AT_dist, XY_data(4).surf_tidal_elevation_align, 50);
XY_data(5).surf_tidal_align_polyfit = polyfit(XY_data(5).AT_dist, XY_data(5).surf_tidal_elevation_align, 50);
XY_data(6).surf_tidal_align_polyfit = polyfit(XY_data(6).AT_dist, XY_data(6).surf_tidal_elevation_align, 50);

XY_data(1).surf_tidal_align_polyval = polyval(XY_data(1).surf_tidal_align_polyfit, XY_data(1).AT_dist);
XY_data(2).surf_tidal_align_polyval = polyval(XY_data(2).surf_tidal_align_polyfit, XY_data(2).AT_dist);
XY_data(3).surf_tidal_align_polyval = polyval(XY_data(3).surf_tidal_align_polyfit, XY_data(3).AT_dist);
XY_data(4).surf_tidal_align_polyval = polyval(XY_data(4).surf_tidal_align_polyfit, XY_data(4).AT_dist);
XY_data(5).surf_tidal_align_polyval = polyval(XY_data(5).surf_tidal_align_polyfit, XY_data(5).AT_dist);
XY_data(6).surf_tidal_align_polyval = polyval(XY_data(6).surf_tidal_align_polyfit, XY_data(6).AT_dist);

XY_data(1).bed_tidal_align_polyfit = polyfit(XY_data(1).AT_dist, XY_data(1).bed_tidal_elevation_align, 50);
XY_data(2).bed_tidal_align_polyfit = polyfit(XY_data(2).AT_dist, XY_data(2).bed_tidal_elevation_align, 50);
XY_data(3).bed_tidal_align_polyfit = polyfit(XY_data(3).AT_dist, XY_data(3).bed_tidal_elevation_align, 50);
XY_data(4).bed_tidal_align_polyfit = polyfit(XY_data(4).AT_dist, XY_data(4).bed_tidal_elevation_align, 50);
XY_data(5).bed_tidal_align_polyfit = polyfit(XY_data(5).AT_dist, XY_data(5).bed_tidal_elevation_align, 50);
XY_data(6).bed_tidal_align_polyfit = polyfit(XY_data(6).AT_dist, XY_data(6).bed_tidal_elevation_align, 50);

XY_data(1).bed_tidal_align_polyval = polyval(XY_data(1).bed_tidal_align_polyfit, XY_data(1).AT_dist);
XY_data(2).bed_tidal_align_polyval = polyval(XY_data(2).bed_tidal_align_polyfit, XY_data(2).AT_dist);
XY_data(3).bed_tidal_align_polyval = polyval(XY_data(3).bed_tidal_align_polyfit, XY_data(3).AT_dist);
XY_data(4).bed_tidal_align_polyval = polyval(XY_data(4).bed_tidal_align_polyfit, XY_data(4).AT_dist);
XY_data(5).bed_tidal_align_polyval = polyval(XY_data(5).bed_tidal_align_polyfit, XY_data(5).AT_dist);
XY_data(6).bed_tidal_align_polyval = polyval(XY_data(6).bed_tidal_align_polyfit, XY_data(6).AT_dist);

% polyfit
fig = figure(115);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).surf_tidal_align_polyval);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).surf_tidal_align_polyval);
plot(XY_data(3).AT_dist/1e3, XY_data(3).surf_tidal_align_polyval);
plot(XY_data(4).AT_dist/1e3, XY_data(4).surf_tidal_align_polyval);
plot(XY_data(5).AT_dist/1e3, XY_data(5).surf_tidal_align_polyval);
plot(XY_data(6).AT_dist/1e3, XY_data(6).surf_tidal_align_polyval);
title('2011,2014,2015,2017,2018,2019 Aligned Polyfit Surface');
ax1.YLim = [-10 110];

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).bed_tidal_align_polyval);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).bed_tidal_align_polyval);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).bed_tidal_align_polyval);
plot(XY_data(4).AT_vel_corrected/1e3, XY_data(4).bed_tidal_align_polyval);
plot(XY_data(5).AT_vel_corrected/1e3, XY_data(5).bed_tidal_align_polyval);
plot(XY_data(6).AT_vel_corrected/1e3, XY_data(6).bed_tidal_align_polyval);
title('2011,2014,2015,2017,2018,2019 Aligned Polyfit Bed');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Elevation WGS (m)');
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'xy')

%% Smoothed 200
XY_data(1).surf_tidal_align_smooth_200 = smooth(XY_data(1).surf_tidal_elevation_align, 200);
XY_data(2).surf_tidal_align_smooth_200 = smooth(XY_data(2).surf_tidal_elevation_align, 200);
XY_data(3).surf_tidal_align_smooth_200 = smooth(XY_data(3).surf_tidal_elevation_align, 200);
XY_data(4).surf_tidal_align_smooth_200 = smooth(XY_data(4).surf_tidal_elevation_align, 200);
XY_data(5).surf_tidal_align_smooth_200 = smooth(XY_data(5).surf_tidal_elevation_align, 200);
XY_data(6).surf_tidal_align_smooth_200 = smooth(XY_data(6).surf_tidal_elevation_align, 200);

XY_data(1).bed_tidal_align_smooth_200 = smooth(XY_data(1).bed_tidal_elevation_align, 200);
XY_data(2).bed_tidal_align_smooth_200 = smooth(XY_data(2).bed_tidal_elevation_align, 200);
XY_data(3).bed_tidal_align_smooth_200 = smooth(XY_data(3).bed_tidal_elevation_align, 200);
XY_data(4).bed_tidal_align_smooth_200 = smooth(XY_data(4).bed_tidal_elevation_align, 200);
XY_data(5).bed_tidal_align_smooth_200 = smooth(XY_data(5).bed_tidal_elevation_align, 200);
XY_data(6).bed_tidal_align_smooth_200 = smooth(XY_data(6).bed_tidal_elevation_align, 200);

% smooth 150
fig = figure(117);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).surf_tidal_align_smooth_200);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).surf_tidal_align_smooth_200);
plot(XY_data(3).AT_dist/1e3, XY_data(3).surf_tidal_align_smooth_200);
plot(XY_data(4).AT_dist/1e3, XY_data(4).surf_tidal_align_smooth_200);
plot(XY_data(5).AT_dist/1e3, XY_data(5).surf_tidal_align_smooth_200);
plot(XY_data(6).AT_dist/1e3, XY_data(6).surf_tidal_align_smooth_200);
title('2011,2014,2015,2017,2018,2019 Aligned 200pt Smooth Surface');
ax1.YLim = [-10 110];
ylabel('Elevation WGS (m)');

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).bed_tidal_align_smooth_200);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).bed_tidal_align_smooth_200);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).bed_tidal_align_smooth_200);
plot(XY_data(4).AT_vel_corrected/1e3, XY_data(4).bed_tidal_align_smooth_200);
plot(XY_data(5).AT_vel_corrected/1e3, XY_data(5).bed_tidal_align_smooth_200);
plot(XY_data(6).AT_vel_corrected/1e3, XY_data(6).bed_tidal_align_smooth_200);
title('2011,2014,2015,2017,2018,2019 Aligned 200pt Smooth Bed');
ylabel('Elevation WGS (m)');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'x')

%% Smoothed 150
XY_data(1).surf_tidal_align_smooth_150 = smooth(XY_data(1).surf_tidal_elevation_align, 150);
XY_data(2).surf_tidal_align_smooth_150 = smooth(XY_data(2).surf_tidal_elevation_align, 150);
XY_data(3).surf_tidal_align_smooth_150 = smooth(XY_data(3).surf_tidal_elevation_align, 150);
XY_data(4).surf_tidal_align_smooth_150 = smooth(XY_data(4).surf_tidal_elevation_align, 150);
XY_data(5).surf_tidal_align_smooth_150 = smooth(XY_data(5).surf_tidal_elevation_align, 150);
XY_data(6).surf_tidal_align_smooth_150 = smooth(XY_data(6).surf_tidal_elevation_align, 150);

XY_data(1).bed_tidal_align_smooth_150 = smooth(XY_data(1).bed_tidal_elevation_align, 150);
XY_data(2).bed_tidal_align_smooth_150 = smooth(XY_data(2).bed_tidal_elevation_align, 150);
XY_data(3).bed_tidal_align_smooth_150 = smooth(XY_data(3).bed_tidal_elevation_align, 150);
XY_data(4).bed_tidal_align_smooth_150 = smooth(XY_data(4).bed_tidal_elevation_align, 150);
XY_data(5).bed_tidal_align_smooth_150 = smooth(XY_data(5).bed_tidal_elevation_align, 150);
XY_data(6).bed_tidal_align_smooth_150 = smooth(XY_data(6).bed_tidal_elevation_align, 150);

% smooth 150
fig = figure(116);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).surf_tidal_align_smooth_150);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).surf_tidal_align_smooth_150);
plot(XY_data(3).AT_dist/1e3, XY_data(3).surf_tidal_align_smooth_150);
plot(XY_data(4).AT_dist/1e3, XY_data(4).surf_tidal_align_smooth_150);
plot(XY_data(5).AT_dist/1e3, XY_data(5).surf_tidal_align_smooth_150);
plot(XY_data(6).AT_dist/1e3, XY_data(6).surf_tidal_align_smooth_150);
title('2011,2014,2015,2017,2018,2019 Aligned 150pt Smooth Surface');
ax1.YLim = [-10 110];
ylabel('Elevation WGS (m)');

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).bed_tidal_align_smooth_150);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).bed_tidal_align_smooth_150);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).bed_tidal_align_smooth_150);
plot(XY_data(4).AT_vel_corrected/1e3, XY_data(4).bed_tidal_align_smooth_150);
plot(XY_data(5).AT_vel_corrected/1e3, XY_data(5).bed_tidal_align_smooth_150);
plot(XY_data(6).AT_vel_corrected/1e3, XY_data(6).bed_tidal_align_smooth_150);
title('2011,2014,2015,2017,2018,2019 Aligned 150pt Smooth Bed');
ylabel('Elevation WGS (m)');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'x')

%% Smooth 100
XY_data(1).surf_tidal_align_smooth_100 = smooth(XY_data(1).surf_tidal_elevation_align, 100);
XY_data(2).surf_tidal_align_smooth_100 = smooth(XY_data(2).surf_tidal_elevation_align, 100);
XY_data(3).surf_tidal_align_smooth_100 = smooth(XY_data(3).surf_tidal_elevation_align, 100);
XY_data(4).surf_tidal_align_smooth_100 = smooth(XY_data(4).surf_tidal_elevation_align, 100);
XY_data(5).surf_tidal_align_smooth_100 = smooth(XY_data(5).surf_tidal_elevation_align, 100);
XY_data(6).surf_tidal_align_smooth_100 = smooth(XY_data(6).surf_tidal_elevation_align, 100);

XY_data(1).bed_tidal_align_smooth_100 = smooth(XY_data(1).bed_tidal_elevation_align, 100);
XY_data(2).bed_tidal_align_smooth_100 = smooth(XY_data(2).bed_tidal_elevation_align, 100);
XY_data(3).bed_tidal_align_smooth_100 = smooth(XY_data(3).bed_tidal_elevation_align, 100);
XY_data(4).bed_tidal_align_smooth_100 = smooth(XY_data(4).bed_tidal_elevation_align, 100);
XY_data(5).bed_tidal_align_smooth_100 = smooth(XY_data(5).bed_tidal_elevation_align, 100);
XY_data(6).bed_tidal_align_smooth_100 = smooth(XY_data(6).bed_tidal_elevation_align, 100);

% smooth 100
fig = figure(117);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).surf_tidal_align_smooth_100);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).surf_tidal_align_smooth_100);
plot(XY_data(3).AT_dist/1e3, XY_data(3).surf_tidal_align_smooth_100);
plot(XY_data(4).AT_dist/1e3, XY_data(4).surf_tidal_align_smooth_100);
plot(XY_data(5).AT_dist/1e3, XY_data(5).surf_tidal_align_smooth_100);
plot(XY_data(6).AT_dist/1e3, XY_data(6).surf_tidal_align_smooth_100);
title('2011,2014,2015,2017,2018,2019 Aligned 100pt Smooth Surface');
ax1.YLim = [-10 110];

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).bed_tidal_align_smooth_100);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).bed_tidal_align_smooth_100);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).bed_tidal_align_smooth_100);
plot(XY_data(4).AT_vel_corrected/1e3, XY_data(4).bed_tidal_align_smooth_100);
plot(XY_data(5).AT_vel_corrected/1e3, XY_data(5).bed_tidal_align_smooth_100);
plot(XY_data(6).AT_vel_corrected/1e3, XY_data(6).bed_tidal_align_smooth_100);
title('2011,2014,2015,2017,2018,2019 Aligned 100pt Smooth Bed');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Elevation WGS (m)');
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'xy')
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
   XY_data(i).Ending_clip = XY_data(i).Beginning_clip(XY_data(i).Beginning_clip <= 6.80e4);
   XY_data(i).Ending_clip = XY_data(i).Ending_clip(XY_data(i).Ending_clip >= -100);
   XY_data(i).Ending_clip = unique(XY_data(i).Ending_clip);
   % Clip Beginning off Each Input Field
   XY_data(i).AT_Dist_Beg = XY_data(i).AT_dist(XY_data(i).find_AT_val:end);
   XY_data(i).Surf_Beg = XY_data(i).Surf(XY_data(i).find_AT_val:end);
   XY_data(i).Bed_Beg = XY_data(i).Bed(XY_data(i).find_AT_val:end);   
   XY_data(i).Thickness_Beg = XY_data(i).Thickness(XY_data(i).find_AT_val:end);
   XY_data(i).Vel_Mag_Beg = XY_data(i).Vel_mag(XY_data(i).find_AT_val:end);
   XY_data(i).Lats_Beg = XY_data(i).Lats(XY_data(i).find_AT_val:end);
   XY_data(i).Lons_Beg = XY_data(i).Lons(XY_data(i).find_AT_val:end);
   XY_data(i).X_Beg = XY_data(i).X(XY_data(i).find_AT_val:end);
   XY_data(i).Y_Beg = XY_data(i).Y(XY_data(i).find_AT_val:end);
   XY_data(i).div_09_10_2000m_Beg = XY_data(i).div_09_10_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_12_13_2000m_Beg = XY_data(i).div_12_13_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_14_15_2000m_Beg = XY_data(i).div_14_15_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_15_16_2000m_Beg = XY_data(i).div_15_16_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_16_17_2000m_Beg = XY_data(i).div_16_17_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).div_17_18_2000m_Beg = XY_data(i).div_17_18_2000m(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_10_Beg = XY_data(i).SMB_10(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_11_Beg = XY_data(i).SMB_11(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_12_Beg = XY_data(i).SMB_12(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_13_Beg = XY_data(i).SMB_13(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_14_Beg = XY_data(i).SMB_14(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_15_Beg = XY_data(i).SMB_15(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_16_Beg = XY_data(i).SMB_16(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_17_Beg = XY_data(i).SMB_17(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_18_Beg = XY_data(i).SMB_18(XY_data(i).find_AT_val:end);
   XY_data(i).SMB_19_Beg = XY_data(i).SMB_19(XY_data(i).find_AT_val:end);   
   XY_data(i).Surf_tidal_elevation_align_Beg = XY_data(i).surf_tidal_elevation_align(XY_data(i).find_AT_val:end);
   XY_data(i).Surf_poly_Beg = XY_data(1).surf_tidal_align_polyval(XY_data(i).find_AT_val:end);
   XY_data(i).Surf_smooth_150_Beg = XY_data(1).surf_tidal_align_smooth_150(XY_data(i).find_AT_val:end);
   XY_data(i).Surf_smooth_100_Beg = XY_data(1).surf_tidal_align_smooth_100(XY_data(i).find_AT_val:end);
   XY_data(i).Bed_tidal_elevation_align_Beg = XY_data(i).bed_tidal_elevation_align(XY_data(i).find_AT_val:end);
   XY_data(i).Bed_poly_Beg = XY_data(1).surf_tidal_align_polyval(XY_data(i).find_AT_val:end);
   XY_data(i).Bed_smooth_150_Beg = XY_data(1).bed_tidal_align_smooth_150(XY_data(i).find_AT_val:end);
   XY_data(i).Bed_smooth_100_Beg = XY_data(1).bed_tidal_align_smooth_100(XY_data(i).find_AT_val:end);
   % Clip Endings off Each Input Field
   XY_data(i).AT_Dist_End = XY_data(i).AT_Dist_Beg(1:length(XY_data(i).Ending_clip));
   XY_data(i).Surf_End = XY_data(i).Surf_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Bed_End = XY_data(i).Bed_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Thickness_End = XY_data(i).Thickness_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Vel_Mag_End = XY_data(i).Vel_Mag_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Lats_End = XY_data(i).Lats_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Lons_End = XY_data(i).Lons_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).X_end = XY_data(i).X_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Y_end = XY_data(i).Y_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_09_10_2000m_End = XY_data(i).div_09_10_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_12_13_2000m_End = XY_data(i).div_12_13_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_14_15_2000m_End = XY_data(i).div_14_15_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_15_16_2000m_End = XY_data(i).div_15_16_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_16_17_2000m_End = XY_data(i).div_16_17_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).div_17_18_2000m_End = XY_data(i).div_17_18_2000m_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_10_End = XY_data(i).SMB_10_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_11_End = XY_data(i).SMB_11_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_12_End = XY_data(i).SMB_12_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_13_End = XY_data(i).SMB_13_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_14_End = XY_data(i).SMB_14_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_15_End = XY_data(i).SMB_15_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_16_End = XY_data(i).SMB_16_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_17_End = XY_data(i).SMB_17_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_17_End = XY_data(i).SMB_17_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_18_End = XY_data(i).SMB_18_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).SMB_19_End = XY_data(i).SMB_19_Beg(1:length(XY_data(i).AT_Dist_End));      
   XY_data(i).Surf_tidal_elevation_align_End = XY_data(i).Surf_tidal_elevation_align_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Surf_poly_End = XY_data(i).Surf_poly_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Surf_smooth_150_End = XY_data(1).Surf_smooth_150_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Surf_smooth_100_End = XY_data(1).Surf_smooth_100_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Bed_tidal_elevation_align_End = XY_data(i).Bed_tidal_elevation_align_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Bed_poly_End = XY_data(i).Bed_poly_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Bed_smooth_150_End = XY_data(1).Bed_smooth_150_Beg(1:length(XY_data(i).AT_Dist_End));
   XY_data(i).Bed_smooth_100_End = XY_data(1).Bed_smooth_100_Beg(1:length(XY_data(i).AT_Dist_End));
   
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
       XY_data(i).AT_Dist_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Surf = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Surf_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Bed = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Bed_End, XY_data(1).query_array, 'makima','extrap');   

   XY_data(i).Interp_Thickness = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Thickness_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Vel_mag = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Vel_Mag_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Lats = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Lats_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Lons = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Lons_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_X = interp1(XY_data(i).Ending_clip,...
       XY_data(i).X_end, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Y = interp1(XY_data(i).Ending_clip,...
       XY_data(i).Y_end, XY_data(1).query_array, 'makima','extrap'); 
   
   XY_data(i).Interp_Div_09_10 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_09_10_2000m_End, XY_data(1).query_array, 'makima','extrap');

   XY_data(i).Interp_Div_12_13 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_12_13_2000m_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Div_14_15 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_14_15_2000m_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Div_15_16 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_15_16_2000m_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_Div_16_17 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_16_17_2000m_End, XY_data(1).query_array, 'makima','extrap'); 

   XY_data(i).Interp_Div_17_18 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).div_17_18_2000m_End, XY_data(1).query_array, 'makima','extrap');   
   
   XY_data(i).Interp_SMB_10 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_10_End, XY_data(1).query_array, 'makima','extrap');   
   
   XY_data(i).Interp_SMB_11 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_11_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_SMB_12 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_12_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_SMB_13 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_13_End, XY_data(1).query_array, 'makima','extrap');   
   
   XY_data(i).Interp_SMB_14 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_14_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_SMB_15 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_15_End, XY_data(1).query_array, 'nearest','extrap');
   
   XY_data(i).Interp_SMB_16 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_16_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_SMB_17 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_17_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_SMB_18 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_18_End, XY_data(1).query_array, 'makima','extrap');
   
   XY_data(i).Interp_SMB_19 = interp1(XY_data(i).Ending_clip,...
       XY_data(i).SMB_19_End, XY_data(1).query_array, 'makima','extrap');   
   
   XY_data(i).Interp_Surf_Corrected = interp1(XY_data(i).Ending_clip,...
        XY_data(i).Surf_tidal_elevation_align_End, XY_data(1).query_array, 'makima','extrap');
    
   XY_data(i).Interp_Surf_Smooth_150 = interp1(XY_data(i).Ending_clip,...
        XY_data(i).Surf_smooth_150_End, XY_data(1).query_array, 'makima','extrap');
    
   XY_data(i).Interp_Surf_Smooth_100 = interp1(XY_data(i).Ending_clip,...
        XY_data(i).Surf_smooth_100_End, XY_data(1).query_array, 'makima','extrap');
    
   XY_data(i).Interp_Bed_Corrected = interp1(XY_data(i).Ending_clip,...
        XY_data(i).Bed_tidal_elevation_align_End, XY_data(1).query_array, 'makima','extrap');
    
   XY_data(i).Interp_Bed_Smooth_150 = interp1(XY_data(i).Ending_clip,...
        XY_data(i).Bed_smooth_150_End, XY_data(1).query_array, 'makima','extrap');
    
   XY_data(i).Interp_Bed_Smooth_100 = interp1(XY_data(i).Ending_clip,...
        XY_data(i).Bed_smooth_100_End, XY_data(1).query_array, 'makima','extrap');  
   
   XY_data(i).Interp_Thickness_Smooth_150 = interp1(XY_data(i).Ending_clip,...
       (XY_data(i).Surf_smooth_150_End - XY_data(i).Bed_smooth_150_End), XY_data(1).query_array, 'makima','extrap');   
   
   XY_data(i).Interp_Thickness_Smooth_100 = interp1(XY_data(i).Ending_clip,...
       (XY_data(i).Surf_smooth_100_End - XY_data(i).Bed_smooth_100_End), XY_data(1).query_array, 'makima','extrap');
end

%% 4 - Time Fraction Component, Tidal Inputs
% Time component fractions for percentages of years between acquisition
% dates of individual radar echograms (2014-2017). N_years variable denotes
% the total number of years between the 2 observed echogram periods, which
% is the sum of all the partial year components used to correct each
% divergence raster for the observed period. 
clear time_offset 

time_offset.P11 = (239/365); % 2011
time_offset.P12 = (366/366); % 2012
time_offset.P13 = (365/365); % 2013
time_offset.P14 = (124/365); % 2014
time_offset.P14_2 = (241/365); % 2014 post echogram
time_offset.P15 = (124/365); % 2015
time_offset.P15_2 = (241/365); % 2015 post echogram
time_offset.P16 = (366/366); % 2016
time_offset.P17 = (89/365); % 2017
time_offset.P17_2 = (276/365); % 2017 post echogram
time_offset.P18 = (94/365); % 2018
time_offset.P18_2 = (271/365); % 2018 post echogram
time_offset.P19 = (107/365); % 2019

% total time offset
time_offset.N_years = (time_offset.P11 + time_offset.P12 + time_offset.P13...
    + time_offset.P14 + time_offset.P14_2 + time_offset.P15 + time_offset.P15_2...
    + time_offset.P16 + time_offset.P17 + time_offset.P17_2 + time_offset.P18...
    + time_offset.P18_2 + time_offset.P19); 

%% SMB 
% Surface Mass Balance Correction
for i = 1:numel(XY_data)
    % 2011 - 2014
    XY_data(i).SMB_11t = (XY_data(i).Interp_SMB_11*time_offset.P11);
    XY_data(i).SMB_12t = (XY_data(i).Interp_SMB_12*time_offset.P12);
    XY_data(i).SMB_13t = (XY_data(i).Interp_SMB_13*time_offset.P13);
    XY_data(i).SMB_14t = (XY_data(i).Interp_SMB_14*time_offset.P14);
    % 2014 - 2015
    XY_data(i).SMB_14_2t = (XY_data(i).Interp_SMB_14*time_offset.P14_2);
    XY_data(i).SMB_15t = (XY_data(i).Interp_SMB_15*time_offset.P15);
    % 2015 - 2017
    XY_data(i).SMB_15_2t = (XY_data(i).Interp_SMB_15*time_offset.P15_2);
    XY_data(i).SMB_16t = (XY_data(i).Interp_SMB_16*time_offset.P16);
    XY_data(i).SMB_17t = (XY_data(i).Interp_SMB_17*time_offset.P17);
    % 2017 - 2018
    XY_data(i).SMB_17_2t = (XY_data(i).Interp_SMB_17*time_offset.P17_2);
    XY_data(i).SMB_18t = (XY_data(i).Interp_SMB_18*time_offset.P18);
    % 2018 - 2019
    XY_data(i).SMB_18_2t = (XY_data(i).Interp_SMB_18*time_offset.P18_2);
    XY_data(i).SMB_19t = (XY_data(i).Interp_SMB_19*time_offset.P19);
    %2011 - 2019 mean
    XY_data(i).SMB_11_19_mean = (XY_data(i).SMB_11t + XY_data(i).SMB_12t + ...
      XY_data(i).SMB_13t + XY_data(i).SMB_14t + XY_data(i).SMB_14_2t + ...
      XY_data(i).SMB_15t + XY_data(i).SMB_15_2t + XY_data(i).SMB_16t + ...
      XY_data(i).SMB_17t + XY_data(i).SMB_17_2t + XY_data(i).SMB_18t + ...
      XY_data(i).SMB_18_2t + XY_data(i).SMB_19t);

    % firn correction (in meters)
    XY_data(i).Firn_thickness = 0;  
end
%% 5 ICE DIVERGENCE 
% Note: For years that are missing (e.g.10_11, 11_12, 13_14, 18_19, 19_20) 
% alternative years are used and multipled by time correction varible. 
for i = 1:numel(XY_data)
    % 2011 - 2014
    XY_data(i).div_11t = (XY_data(i).Interp_Div_09_10*time_offset.P11);
    XY_data(i).div_12t = (XY_data(i).Interp_Div_12_13*time_offset.P12);
    XY_data(i).div_13t = (XY_data(i).Interp_Div_12_13*time_offset.P13);
    XY_data(i).div_14t = (XY_data(i).Interp_Div_14_15*time_offset.P14);
    % 2014 - 2015
    XY_data(i).div_14_2t = (XY_data(i).Interp_Div_14_15*time_offset.P14_2);
    XY_data(i).div_15t = (XY_data(i).Interp_Div_15_16*time_offset.P15);
    % 2015 - 2017
    XY_data(i).div_15_2t = (XY_data(i).Interp_Div_15_16*time_offset.P15_2);
    XY_data(i).div_16t = (XY_data(i).Interp_Div_16_17*time_offset.P16);
    XY_data(i).div_17t = (XY_data(i).Interp_Div_17_18*time_offset.P17);    
    % 2017 - 2018
    XY_data(i).div_17_2t = (XY_data(i).Interp_Div_17_18*time_offset.P17_2);
    XY_data(i).div_18t = (XY_data(i).Interp_Div_17_18*time_offset.P18);
    % 2018 - 2019
    XY_data(i).div_18_2t = (XY_data(i).Interp_Div_17_18*time_offset.P18_2);
    XY_data(i).div_19t = (XY_data(i).Interp_Div_17_18*time_offset.P19);    
    % 2011 - 2019
    XY_data(i).div_11_19_mean = (XY_data(i).div_11t + XY_data(i).div_12t + ...
      XY_data(i).div_13t + XY_data(i).div_14t + XY_data(i).div_14_2t + ...
      XY_data(i).div_15t + XY_data(i).div_15_2t + XY_data(i).div_16t + ...
      XY_data(i).div_17t + XY_data(i).div_17_2t + XY_data(i).div_18t + ...
      XY_data(i).div_18_2t + XY_data(i).div_19t);
end

% Calcuate Mean Thicknesses from paired years
for i = 1:numel(XY_data)
   XY_data(i).Mean_Thickness = ((XY_data(i).Interp_Surf_Corrected - XY_data(i).Interp_Bed_Corrected) + ...
       (XY_data(i+1).Interp_Surf_Corrected - XY_data(i+1).Interp_Bed_Corrected))/2;
   %XY_data(i).Mean_Thickness_Poly = (XY_data(i).Interp_Thickness_Poly + XY_data(i+1).Interp_Thickness_Poly)/2;
   XY_data(i).Mean_Thickness_Smooth150 = (XY_data(i).Interp_Thickness_Smooth_150 + XY_data(i+1).Interp_Thickness_Smooth_150)/2;
   XY_data(i).Mean_Thickness_Smooth100 = (XY_data(i).Interp_Thickness_Smooth_100 + XY_data(i+1).Interp_Thickness_Smooth_100)/2;
   
   XY_data(i).Mean_Thickness_smooth = movmean((XY_data(i+1).Interp_Thickness + XY_data(i+1).Interp_Thickness)/2, 10,'omitnan');
   
   XY_data(i).Bed_diff = XY_data(i+1).Interp_Bed - XY_data(i).Interp_Bed;
   XY_data(i).Bed_diff_smooth = movmean((XY_data(i+1).Interp_Bed - XY_data(i).Interp_Bed),10,'omitnan'); 
   
   XY_data(i).Surf_diff = XY_data(i+1).Interp_Surf - XY_data(i).Interp_Surf;
   XY_data(i).Surf_diff_smooth = movmean(XY_data(i+1).Interp_Surf,10,'omitnan') - ...
       movmean(XY_data(i+1).Interp_Surf,10,'omitnan');
   
   XY_data(i).Melt_no_SMB = XY_data(i).Bed_diff - (XY_data(i).div_11_19_mean.*XY_data(i).Mean_Thickness);   
   XY_data(i).Dynamic_component = (XY_data(i).div_11_19_mean.*XY_data(i).Mean_Thickness);
   if i > i(end)
        break
   end
end

%% ANNUAL MELT RATES 
% 2011 - 2014: Melt, Dynamic Thinning, Mean SMB
% ask leigh if I should use total time difference or just divide by number
% of input components
% Removed time_offset, mean is just the sum of the years
% original XY_data(1).Mean_SMB = ((XY_data(1).SMB_10t + XY_data(1).SMB_11t)/(time_offset.P10 + time_offset.P11));
% Time Offset
XY_data(1).total_time_offset = time_offset.P11 + time_offset.P12 + time_offset.P13 + time_offset.P14;

% SMB
XY_data(1).Mean_SMB = (XY_data(1).SMB_11t + XY_data(1).SMB_12t + XY_data(1).SMB_13t + XY_data(1).SMB_14t);

% Dynamic Thinning
XY_data(1).Dynamic_component = ((XY_data(1).div_11t + XY_data(1).div_12t + XY_data(1).div_13t + XY_data(1).div_14t)...
    .*XY_data(1).Mean_Thickness);
%XY_data(1).Dynamic_component_Poly = ((XY_data(1).div_11t + XY_data(1).div_12t + XY_data(1).div_13t + XY_data(1).div_14t)...
%    .*XY_data(1).Mean_Thickness_Poly);
XY_data(1).Dynamic_component_S150 = ((XY_data(1).div_11t + XY_data(1).div_12t + XY_data(1).div_13t + XY_data(1).div_14t)...
    .*XY_data(1).Mean_Thickness_Smooth150);
XY_data(1).Dynamic_component_S100 = ((XY_data(1).div_11t + XY_data(1).div_12t + XY_data(1).div_13t + XY_data(1).div_14t)...
    .*XY_data(1).Mean_Thickness_Smooth100);

% Total Melt
XY_data(1).Melt_SMB = (XY_data(1).Interp_Thickness - XY_data(2).Interp_Thickness) - XY_data(1).Dynamic_component - XY_data(1).Mean_SMB;

% Annual Melt 
XY_data(1).Melt_SMB_yr = ((XY_data(1).Interp_Thickness - XY_data(2).Interp_Thickness) - XY_data(1).Dynamic_component - XY_data(1).Mean_SMB)...
    /XY_data(1).total_time_offset;
%XY_data(1).Melt_SMB_poly_yr = ((XY_data(1).Interp_Thickness_Poly - XY_data(2).Interp_Thickness_Poly) - XY_data(1).Dynamic_component_Poly - XY_data(1).Mean_SMB)...
%    /XY_data(1).total_time_offset;
XY_data(1).Melt_SMB_smooth150_yr = ((XY_data(1).Interp_Thickness_Smooth_150 - XY_data(2).Interp_Thickness_Smooth_150) - XY_data(1).Dynamic_component_S150 - XY_data(1).Mean_SMB)...
    /XY_data(1).total_time_offset;
XY_data(1).Melt_SMB_smooth100_yr = ((XY_data(1).Interp_Thickness_Smooth_100 - XY_data(2).Interp_Thickness_Smooth_100) - XY_data(1).Dynamic_component_S100 - XY_data(1).Mean_SMB)...
    /XY_data(1).total_time_offset;


fig1 = figure(101);
subplot(5,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(1).Mean_SMB);
%ylim([-2.5 -1]);
ylabel('SMB (m/yr)');
title('Mean SMB')
legend('Time offset','Location','eastoutside');
subplot(5,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(1).Dynamic_component);
%ylim([-20 10]);
ylabel('Thinning Rate (m/yr)');
title('Mean Dynamic Thinning')
legend('Time offset','Location','eastoutside');
subplot(5,1,3)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB);
ylim([-100 200]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,4)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr);
ylim([-50 100]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,5)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 100]);
suptitle('Petermann Line 1 2011-2014')
legend('2011', '2014', 'Location','southeast');
han = axes(fig1,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

%%
figure(122)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','g');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','c');
title('Elevation Profiles');
ylabel('elevation (m)');
suptitle('Petermann Line 4 2010-2011')
legend('2010', '2011', 'Location','southeast');

figure(123)
plot(XY_data(1).AT_dist/1e3, XY_data(1).Surf, 'color','r');
hold on
plot(XY_data(1).AT_dist/1e3, XY_data(2).Surf, 'color','b');
plot(XY_data(1).AT_dist/1e3, XY_data(1).Bed, 'color','r');
plot(XY_data(1).AT_dist/1e3, XY_data(2).Bed, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
suptitle('Petermann Line 4 2010-2011')
legend('2010', '2011', 'Location','southeast');

%%
fig1 = figure(1001);
subplot(3,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_yr, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_poly_yr, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_smooth50_yr, 'color','g');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Melt_SMB_smooth100_yr, 'color','k');
xlim([-0 50]);
ylim([-50 100]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Original', 'Poly fit', 'Smooth150', 'Smooth100', 'Location','eastoutside');
subplot(3,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Poly, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Smooth_50, 'color','g');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Smooth_100, 'color','k');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Poly, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Smooth_50, 'color','g');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Smooth_100, 'color','k');
xlim([-0 50]);
ylim([-500 150]);
title('Elevation Profile 1');
legend('Original', 'Poly fit', 'Smooth150', 'Smooth100', 'Location','eastoutside');
subplot(3,1,3)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Poly, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Smooth_50, 'color','g');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Smooth_100, 'color','k');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Poly, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Smooth_50, 'color','g');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Smooth_100, 'color','k');
title('Elevation Profile 2');
ylabel('elevation (m)');
xlim([-0 50]);
ylim([-500 150]);
suptitle('Petermann Line 1 2011-2014')
legend('2011', '2014', 'Location','southeast');
han = axes(fig1,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

%% 2014 - 2015: Melt, Dynamic Thinning, Mean SMB
% original XY_data(2).Mean_SMB = ((XY_data(2).SMB_11_2t + XY_data(2).SMB_12t + XY_data(2).SMB_13t)/(time_offset.P11_2 + time_offset.P12 + time_offset.P13));
% Total Time offset
XY_data(2).total_time_offset = time_offset.P14_2 + time_offset.P15;

% SMB
XY_data(2).Mean_SMB = (XY_data(2).SMB_14_2t + XY_data(2).SMB_15t);

% Dynamic Thinning
XY_data(2).Dynamic_component = ((XY_data(2).div_14_2t + XY_data(2).div_15t)...
    .*XY_data(2).Mean_Thickness);
%XY_data(2).Dynamic_component_Poly = ((XY_data(2).div_14_2t + XY_data(2).div_15t)...
%    .*XY_data(2).Mean_Thickness_Poly);
XY_data(2).Dynamic_component_S150 = ((XY_data(2).div_14_2t + XY_data(2).div_15t)...
    .*XY_data(2).Mean_Thickness_Smooth150);
XY_data(2).Dynamic_component_S100 = ((XY_data(2).div_14_2t + XY_data(2).div_15t)...
    .*XY_data(2).Mean_Thickness_Smooth100);

% Total Melt
XY_data(2).Melt_SMB = (XY_data(2).Interp_Thickness - XY_data(3).Interp_Thickness) - XY_data(2).Dynamic_component - XY_data(2).Mean_SMB;

% Annual Melt 
XY_data(2).Melt_SMB_yr = ((XY_data(2).Interp_Thickness - XY_data(3).Interp_Thickness) - XY_data(2).Dynamic_component - XY_data(2).Mean_SMB)...
    /XY_data(2).total_time_offset;
%XY_data(2).Melt_SMB_poly_yr = ((XY_data(2).Interp_Thickness_Poly - XY_data(3).Interp_Thickness_Poly) - XY_data(2).Dynamic_component_Poly - XY_data(2).Mean_SMB)...
%    /XY_data(2).total_time_offset;
XY_data(2).Melt_SMB_smooth150_yr = ((XY_data(2).Interp_Thickness_Smooth_150 - XY_data(3).Interp_Thickness_Smooth_150) - XY_data(2).Dynamic_component_S150 - XY_data(2).Mean_SMB)...
    /XY_data(2).total_time_offset;
XY_data(2).Melt_SMB_smooth100_yr = ((XY_data(2).Interp_Thickness_Smooth_100 - XY_data(3).Interp_Thickness_Smooth_100) - XY_data(2).Dynamic_component_S100 - XY_data(2).Mean_SMB)...
    /XY_data(2).total_time_offset;


fig2 = figure(102);
subplot(5,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(2).Mean_SMB);
%ylim([-2.5 -1]);
ylabel('SMB (m/yr)');
title('Mean SMB')
legend('Time offset','Location','eastoutside');
subplot(5,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(2).Dynamic_component);
%ylim([-20 10]);
ylabel('Thinning Rate (m/yr)');
title('Mean Dynamic Thinning')
legend('Time offset','Location','eastoutside');
subplot(5,1,3)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB);
%ylim([-10 100]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,4)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt_SMB_yr);
%ylim([-10 100]);
ylabel('Basal Melt (m)');
title('Annual Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,5)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
suptitle('Petermann Line 1 2014-2015')
legend('2014', '2015', 'Location','southeast');
han = axes(fig2,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

%% 2015 - 2017: Melt, Dynamic Thinning, Mean SMB
% Time offset
XY_data(3).total_time_offset = time_offset.P15_2 + time_offset.P16 + time_offset.P17;

% SMB
XY_data(3).Mean_SMB = (XY_data(3).SMB_15_2t + XY_data(1).SMB_16t + XY_data(1).SMB_17t);

% Dynamic Thinning
XY_data(3).Dynamic_component = ((XY_data(3).div_15_2t + XY_data(3).div_16t + XY_data(3).div_17t)...
    .*XY_data(3).Mean_Thickness);
XY_data(3).Dynamic_component_S150 = ((XY_data(3).div_15_2t + XY_data(3).div_16t + XY_data(3).div_17t)...
    .*XY_data(3).Mean_Thickness_Smooth150);
XY_data(3).Dynamic_component_S100 = ((XY_data(3).div_15_2t + XY_data(3).div_16t + XY_data(3).div_17t)...
    .*XY_data(3).Mean_Thickness_Smooth100);

% Total Melt
XY_data(3).Melt_SMB = (XY_data(3).Interp_Thickness - XY_data(4).Interp_Thickness) - XY_data(3).Dynamic_component - XY_data(3).Mean_SMB;

% Annual Melt 
XY_data(3).Melt_SMB_yr = ((XY_data(3).Interp_Thickness - XY_data(4).Interp_Thickness) - XY_data(3).Dynamic_component - XY_data(3).Mean_SMB)...
    /XY_data(3).total_time_offset;
XY_data(3).Melt_SMB_smooth150_yr = ((XY_data(3).Interp_Thickness_Smooth_150 - XY_data(4).Interp_Thickness_Smooth_150) - XY_data(3).Dynamic_component_S150 - XY_data(3).Mean_SMB)...
    /XY_data(3).total_time_offset;
XY_data(3).Melt_SMB_smooth100_yr = ((XY_data(3).Interp_Thickness_Smooth_100 - XY_data(4).Interp_Thickness_Smooth_100) - XY_data(3).Dynamic_component_S100 - XY_data(3).Mean_SMB)...
    /XY_data(3).total_time_offset;

fig3 = figure(103);
subplot(5,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(3).Mean_SMB);
%ylim([-2.2 -1]);
ylabel('SMB (m/yr)');
title('Mean SMB')
legend('Time offset','Location','eastoutside');
subplot(5,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(3).Dynamic_component);
%ylim([-20 10]);
ylabel('Thinning Rate (m/yr)');
title('Mean Dynamic Thinning')
legend('Time offset','Location','eastoutside');
subplot(5,1,3)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Melt_SMB);
%ylim([-20 200]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,4)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Melt_SMB_yr);
ylim([-100 200]);
ylabel('Basal Melt (m)');
title('Annual Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,5)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(3).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 200]);
suptitle('Petermann Line 1 2015-2017')
legend('2015', '2017', 'Location','southeast');
han = axes(fig3,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');

%% 2017 - 2018: Melt, Dynamic Thinning, Mean SMB
% Time Offset
XY_data(4).total_time_offset = time_offset.P17_2 + time_offset.P18;

% SMB
XY_data(4).Mean_SMB = (XY_data(4).SMB_14_2t + XY_data(4).SMB_15t + XY_data(4).SMB_16t + XY_data(4).SMB_17t);

% Dynamic Thinning
XY_data(4).Dynamic_component = ((XY_data(4).div_14_2t + XY_data(4).div_15t + XY_data(4).div_16t + XY_data(4).div_17t)...
    .*XY_data(4).Mean_Thickness);
XY_data(4).Dynamic_component_S150 = ((XY_data(4).div_14_2t + XY_data(4).div_15t + XY_data(4).div_16t + XY_data(4).div_17t)...
    .*XY_data(4).Mean_Thickness_Smooth150);
XY_data(4).Dynamic_component_S100 = ((XY_data(4).div_14_2t + XY_data(4).div_15t + XY_data(4).div_16t + XY_data(4).div_17t)...
    .*XY_data(4).Mean_Thickness_Smooth100);

% Total Melt
XY_data(4).Melt_SMB = ((XY_data(4).Interp_Thickness - XY_data(5).Interp_Thickness) - XY_data(4).Dynamic_component - XY_data(4).Mean_SMB);

% Annual Melt 
XY_data(4).Melt_SMB_yr = ((XY_data(4).Interp_Thickness - XY_data(5).Interp_Thickness) - XY_data(4).Dynamic_component - XY_data(4).Mean_SMB)...
    /XY_data(4).total_time_offset;
XY_data(4).Melt_SMB_smooth150_yr = ((XY_data(4).Interp_Thickness_Smooth_150 - XY_data(5).Interp_Thickness_Smooth_150) - XY_data(4).Dynamic_component_S150 - XY_data(4).Mean_SMB)...
    /XY_data(4).total_time_offset;
XY_data(4).Melt_SMB_smooth100_yr = ((XY_data(4).Interp_Thickness_Smooth_100 - XY_data(5).Interp_Thickness_Smooth_100) - XY_data(4).Dynamic_component_S100 - XY_data(4).Mean_SMB)...
    /XY_data(4).total_time_offset;

fig4 = figure(104);
%tiledlayout(4,1)
subplot(5,1,1)
%ax1 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(4).Mean_SMB);
%ylim([-2.2 -1]);
ylabel('SMB (m/yr)');
title('Mean SMB')
legend('Time offset','Location','eastoutside');
subplot(5,1,2)
%ax2 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(4).Dynamic_component);
%ylim([-50 20]);
ylabel('Thinning Rate (m/yr)');
title('Mean Dynamic Thinning')
legend('Time offset','Location','eastoutside');
subplot(5,1,3)
%ax3 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB);
%ylim([-50 250]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,4)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB_yr);
%ylim([-50 250]);
ylabel('Basal Melt (m)');
title('Annual Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,5)
%ax4 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 250]);
suptitle('Petermann Line 1 2017-2018')
legend('2017', '2018', 'Location','southeast');
han = axes(fig4,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2 ax3 ax4],'x');
ax1.XLim = [0 50];

%% 2018 - 2019: Melt, Dynamic Thinning, Mean SMB
% Time Offset
XY_data(5).total_time_offset = time_offset.P18_2 + time_offset.P19;

% SMB
XY_data(5).Mean_SMB = (XY_data(5).SMB_18_2t + XY_data(5).SMB_19t);

% Dynamic Thinning
XY_data(5).Dynamic_component = ((XY_data(5).div_18_2t + XY_data(5).div_19t).*XY_data(5).Mean_Thickness);
XY_data(5).Dynamic_component_S150 = ((XY_data(5).div_18_2t + XY_data(5).div_19t).*XY_data(5).Mean_Thickness_Smooth150);
XY_data(5).Dynamic_component_S100 = ((XY_data(5).div_18_2t + XY_data(5).div_19t).*XY_data(5).Mean_Thickness_Smooth100);

% Total Melt
XY_data(5).Melt_SMB = ((XY_data(5).Interp_Thickness - XY_data(6).Interp_Thickness) - XY_data(5).Dynamic_component - XY_data(5).Mean_SMB);

% Annual Melt 
XY_data(5).Melt_SMB_yr = ((XY_data(5).Interp_Thickness - XY_data(6).Interp_Thickness) - XY_data(5).Dynamic_component - XY_data(5).Mean_SMB)...
    /XY_data(5).total_time_offset;
XY_data(5).Melt_SMB_smooth150_yr = ((XY_data(5).Interp_Thickness_Smooth_150 - XY_data(6).Interp_Thickness_Smooth_150) - XY_data(5).Dynamic_component_S150 - XY_data(5).Mean_SMB)...
    /XY_data(5).total_time_offset;
XY_data(5).Melt_SMB_smooth100_yr = ((XY_data(5).Interp_Thickness_Smooth_100 - XY_data(6).Interp_Thickness_Smooth_100) - XY_data(5).Dynamic_component_S100 - XY_data(5).Mean_SMB)...
    /XY_data(5).total_time_offset;

fig5 = figure(105);
%tiledlayout(4,1)
subplot(5,1,1)
%ax1 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(4).Mean_SMB);
%ylim([-2.2 -1]);
ylabel('SMB (m/yr)');
title('Mean SMB')
legend('Time offset','Location','eastoutside');
subplot(5,1,2)
%ax2 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(4).Dynamic_component);
%ylim([-50 20]);
ylabel('Thinning Rate (m/yr)');
title('Mean Dynamic Thinning')
legend('Time offset','Location','eastoutside');
subplot(5,1,3)
%ax3 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB);
%ylim([-50 250]);
ylabel('Basal Melt (m)');
title('Total Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,4)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Melt_SMB_yr);
%ylim([-50 250]);
ylabel('Basal Melt (m)');
title('Annual Basal Melt');
legend('Time offset thickness','Location','eastoutside');
subplot(5,1,5)
%ax4 = nexttile;
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(4).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Bed_Corrected, 'color','b');
title('Elevation Profiles');
ylabel('elevation (m)');
ylim([-500 250])
suptitle('Petermann Line 1 2018-2019')
legend('2018', '2019', 'Location','southeast');
han = axes(fig4,'visible','off');
han.XLabel.Visible ='on';
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2 ax3 ax4],'x');
ax1.XLim = [0 50];

%% Define export file structure with fields of interest
for i = 1:5
export_melt_rates(i).data_frame = cat(2, XY_data(i).Interp_AT_Dist.', XY_data(i).Interp_Lons.', XY_data(i).Interp_Lats.',...
    XY_data(i).Interp_X.', XY_data(i).Interp_Y.', XY_data(i).Melt_SMB.', XY_data(i).Melt_SMB_yr.', XY_data(i).Melt_SMB_smooth150_yr.', ...
    XY_data(i).Melt_SMB_smooth100_yr.', XY_data(i).Interp_Surf_Corrected.', XY_data(i).Interp_Bed_Corrected.', XY_data(i).Interp_Thickness.', ...
    XY_data(i).Mean_SMB.', XY_data(i).Dynamic_component.');

export_melt_rates(i).data_frame2 = cat(2, XY_data(i).Interp_AT_Dist.', XY_data(i).Interp_Lons.', XY_data(i).Interp_Lats.',...
    XY_data(i).Interp_X.', XY_data(i).Interp_Y.', XY_data(i).Interp_Vel_mag.', XY_data(i).Interp_SMB_10.', XY_data(i).Interp_SMB_11.',...
    XY_data(i).Interp_SMB_12.', XY_data(i).Interp_SMB_13.', XY_data(i).Interp_SMB_14.', XY_data(i).Interp_SMB_15.', ... 
    XY_data(i).Interp_SMB_16.', XY_data(i).Interp_SMB_17.', XY_data(i).Interp_SMB_18.', XY_data(i).Interp_SMB_19.',...
    XY_data(i).Melt_SMB.', XY_data(i).Melt_SMB_yr.', XY_data(i).Melt_SMB_smooth150_yr.', ...
    XY_data(i).Melt_SMB_smooth100_yr.', XY_data(i).Interp_Surf_Corrected.', XY_data(i).Interp_Bed_Corrected.', XY_data(i).Interp_Thickness.', ...
    XY_data(i).Mean_SMB.', XY_data(i).Dynamic_component.');
end

%% Define Header Array of strings and vertically concatenate to data 
cheader = {'AT Distance','Lons', 'Lats', 'X', 'Y', 'total melt', 'annual melt', 'annual melt S150', ...
    'annual melt S100','surface','Bed', 'Thickness', 'SMB', 'Dynamic thinning'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

% change folder
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Melt_CSV\'

% Change file name per file prior to exporting, comment out excess file
% names if there are fewer years than 5 being exported to CSV

for j = 1:numel(export_melt_rates)
    name(j).file = append('P1M_melt_E_', string(XY_data(j).year),'_',string(XY_data(j+1).year),'.csv');
    fid = fopen(name(j).file,'w');
    fprintf(fid,'%s\n',textHeader);
    fclose(fid);
    dlmwrite(name(j).file, export_melt_rates(j).data_frame, '-append');
end
%% Define Header Array of strings and vertically concatenate to data 2
cheader1 = {'AT Distance','Lons', 'Lats', 'X', 'Y', 'Vel Mag', 'SMB_10', 'SMB_11', 'SMB_12','SMB_13',...
    'SMB_14','SMB_15','SMB_16','SMB_17','SMB_18','SMB_19','total melt', 'annual melt', 'annual melt S150', ...
    'annual melt S100','surface','Bed', 'Thickness', 'SMB', 'Dynamic thinning'}; % header
commaHeader1 = [cheader1;repmat({','},1,numel(cheader1))];
commaHeader1 = commaHeader1(:)';
textHeader1 = cell2mat(commaHeader1);

% change folder
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Melt_CSV\'

% Change file name per file prior to exporting, comment out excess file
% names if there are fewer years than 5 being exported to CSV

for j = 1:numel(export_melt_rates)
    name(j).file = append('P1M_melt_E_new_', string(XY_data(j).year),'_',string(XY_data(j+1).year),'.csv');
    fid = fopen(name(j).file,'w');
    fprintf(fid,'%s\n',textHeader1);
    fclose(fid);
    dlmwrite(name(j).file, export_melt_rates(j).data_frame2, '-append');
end




%% Test figure ti check surface alignment
fig = figure(113);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).Surf_tidal_elevation_align);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).Surf_tidal_elevation_align);
plot(XY_data(3).AT_dist/1e3, XY_data(3).Surf_tidal_elevation_align);
title('2013,2013,2017 Unaligned surfaces');
ax1.YLim = [-10 110];

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).Surf_tidal_elevation_align);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).Surf_tidal_elevation_align);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).Surf_tidal_elevation_align);
title('2013,2013,2017 Aligned surfaces');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Elevation WGS (m)');
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'xy')

fig = figure(114);
tiledlayout(2,1)
ax1 = nexttile;
plot(XY_data(1).AT_dist/1e3, XY_data(1).Bed_tidal_elevation_align);
hold on
plot(XY_data(2).AT_dist/1e3, XY_data(2).Bed_tidal_elevation_align);
plot(XY_data(3).AT_dist/1e3, XY_data(3).Bed_tidal_elevation_align);
title('2013,2013,2017 Unaligned surfaces');
%ax1.YLim = [-10 110];

ax2 = nexttile;
plot(XY_data(1).AT_vel_corrected/1e3, XY_data(1).Bed_tidal_elevation_align);
hold on
plot(XY_data(2).AT_vel_corrected/1e3, XY_data(2).Bed_tidal_elevation_align);
plot(XY_data(3).AT_vel_corrected/1e3, XY_data(3).Bed_tidal_elevation_align);
title('2013,2013,2017 Aligned surfaces');

han=axes(fig,'visible','off'); 
han.Title.Visible='on';
han.XLabel.Visible='on';
han.YLabel.Visible='on';
ylabel(han,'Elevation WGS (m)');
xlabel(han,'Along Track Distance (km)');
linkaxes([ax1 ax2],'xy')

%% TEST PLOTS FOR CLIPPING, IDENTIFY CLIPPING POINTS HERE
% 2010
figure(109)
plot(XY_data(1).AT_dist/1e3, smooth(XY_data(1).Surf,'rloess'));
hold on 
plot(XY_data(1).AT_dist()/1e3, XY_data(1).Bed());
plot(XY_data(1).AT_dist()/1e3, XY_data(1).Surf());
xlabel('AT DIST KM');
ylabel('Surface Elev m');
title('2010');
% 2011
figure(110)
plot(XY_data(2).AT_dist/1e3, smooth(XY_data(2).Surf,'rloess'));
hold on 
plot(XY_data(2).AT_dist()/1e3, XY_data(2).Bed());
plot(XY_data(2).AT_dist()/1e3, XY_data(2).Surf());
xlabel('AT DIST KM');
ylabel('Surface Elev m');
title('2011');
%% 2013
figure(111)
plot(XY_data(3).AT_dist/1e3, smooth(XY_data(3).Surf,'rloess'));
hold on 
plot(XY_data(3).AT_dist(1:3487)/1e3, XY_data(3).Bed(1:3487));
plot(XY_data(3).AT_dist(3493:end)/1e3, XY_data(3).Surf(3493:end));
xlabel('AT DIST KM');
ylabel('Surface Elev m');
title('2013');
% 2014
figure(112)
plot(XY_data(4).AT_dist/1e3, smooth(XY_data(4).Surf,'rloess'));
hold on 
plot(XY_data(4).AT_dist(1:3570)/1e3, XY_data(4).Bed(1:3570));
plot(XY_data(4).AT_dist(3569:end)/1e3, XY_data(4).Surf(3569:end));
xlabel('AT DIST KM');
ylabel('Surface Elev m');
title('2014');
% 2017
figure(113)
plot(XY_data(5).AT_dist/1e3, smooth(XY_data(5).Surf,'rloess'));
hold on 
plot(XY_data(5).AT_dist()/1e3, XY_data(5).Bed());
plot(XY_data(5).AT_dist()/1e3, XY_data(5).Surf());
xlabel('AT DIST KM');
ylabel('Surface Elev m');
title('2017');
%%
figure(114)
plot((XY_data(1).AT_dist + XY_data(1).Vel_mag)/1e3, XY_data(1).Surf);
hold on 
plot((XY_data(1).AT_dist + XY_data(1).Vel_mag)/1e3, XY_data(1).Bed);
plot((XY_data(2).AT_dist + XY_data(2).Vel_mag)/1e3, XY_data(2).Surf);
plot((XY_data(2).AT_dist + XY_data(2).Vel_mag)/1e3, XY_data(2).Bed);
plot((XY_data(3).AT_dist + XY_data(3).Vel_mag)/1e3, XY_data(3).Surf);
plot((XY_data(3).AT_dist + XY_data(3).Vel_mag)/1e3, XY_data(3).Bed);
plot((XY_data(4).AT_dist + XY_data(4).Vel_mag)/1e3, XY_data(4).Surf);
plot((XY_data(4).AT_dist + XY_data(4).Vel_mag)/1e3, XY_data(4).Bed);
plot((XY_data(5).AT_dist + XY_data(5).Vel_mag)/1e3, XY_data(5).Surf);
plot((XY_data(5).AT_dist + XY_data(5).Vel_mag)/1e3, XY_data(5).Bed);
legend show 

figure(115)
plot(XY_data(1).AT_dist/1e3, XY_data(1).Surf);
hold on 
plot(XY_data(1).AT_dist/1e3, XY_data(1).Bed);
plot(XY_data(1).AT_dist/1e3, XY_data(2).Surf);
plot(XY_data(1).AT_dist/1e3, XY_data(2).Bed);
plot(XY_data(1).AT_dist/1e3, XY_data(3).Surf);
plot(XY_data(1).AT_dist/1e3, XY_data(3).Bed);
plot(XY_data(1).AT_dist/1e3, XY_data(4).Surf);
plot(XY_data(1).AT_dist/1e3, XY_data(4).Bed);
plot(XY_data(1).AT_dist/1e3, XY_data(5).Surf);
plot(XY_data(1).AT_dist/1e3, XY_data(5).Bed);
legend show

%% PLOTS 2010 - 2017 melt
XY_data(6).Mean_Thickness = (XY_data(1).Interp_Thickness + XY_data(5).Interp_Thickness)/2;
XY_data(6).Bed_diff = XY_data(5).Interp_Bed - XY_data(1).Interp_Bed;
XY_data(6).Surf_diff = XY_data(5).Interp_Surf - XY_data(1).Interp_Surf;
XY_data(6).Melt_SMB = XY_data(6).Bed_diff - (XY_data(5).div_10_17_mean.*XY_data(6).Mean_Thickness) - XY_data(1).SMB_10_17_mean;
XY_data(6).Dynamic_component = (XY_data(5).div_10_17_mean.*XY_data(6).Mean_Thickness);

figure(1)
subplot(4,1,1)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(1).SMB_10_17_mean/time_offset.N_years);
ylabel('mean SMB (m/yr)');
title('mean annual SMB')

subplot(4,1,2)
plot(XY_data(1).Interp_AT_Dist/1e3,  XY_data(6).Dynamic_component/time_offset.N_years);
ylabel('mean Dyn. Thin. (m/yr)');
title('mean annual dynamic thinning')

subplot(4,1,3)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(6).Melt_SMB/time_offset.N_years);
ylabel('Basal Melt (m)');

subplot(4,1,4)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Surf_Corrected, 'color','r');
hold on
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Surf_Corrected, 'color','b');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(1).Interp_Bed_Corrected, 'color','r');
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(5).Interp_Bed_Corrected, 'color','b');
xlabel('Along Track Distance (km)');
title('Petermann Line 4 2010-2017 Basal Melt');
%% Test Figure Melt and dynamic component 2010-2017
figure(1)
plot(XY_data(1).Interp_AT_Dist/1e3, XY_data(2).Melt);
xlabel('Along Track Distance (km)');
ylabel('Basal Melt (m)');
title('Petermann Line 2 2014-2017 Basal Melt');

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
