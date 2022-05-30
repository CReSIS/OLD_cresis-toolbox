%% Basal melt rate calculation script developed using methods described in 
% Shean et al, 2017 section 3 equation 5. Basal melt rates are calculated
% using six separate inputs, namely 1) change in surface elevation over
% time, 2) annual surface elevation (lidar), 3) along-flow ice velocity, 4)
% density of water, 5) density of ice, and 6) RACMO-derived surface mass 
% balance value arrays.

%% Define Constants
rho_ice = 917; %kg/m^3
rho_water = 1000; %kg/m^3
firn_air_correction = 0; 
water_ice_density_ratio = rho_water/(rho_water - rho_ice);

%% load in arrays SMB & SEC data arrays


%% Obtain lidar & radar surfaces for all years, correct for firn air 
% correction (FAC) and then produce radar and lidar derived bed and
% thickness arrays for given years. Additionally all velocity data saved
for i = 1:numel(pass)
lidar_surface(i).surf = pass(i).layers(3).layer_elev(1:end); 
lidar_surface_FAC(i).surf = lidar_surface(i).surf - firn_air_correction; 
radar_surface(i).surf = pass(i).layers(1).layer_elev(1:end);
radar_surface_FAC(i).surf = radar_surface(i).surf - firn_air_correction;
radar_bed(i).bed = pass(i).layers(2).layer_elev(1:end);
thickness_lidar_FAC(i).layer = lidar_surface_FAC(i).surf - radar_bed(i).bed;
thickness_radar_FAC(i).layer = radar_surface_FAC(i).surf - radar_bed(i).bed;
pass_velocity(i).vel = pass(i).vel(1:end);
end

% Load/write in SMB structure
SMB = zeros(size(thickness_lidar_FAC(1).layer));

for i = 1:numel(pass)
    lidar_SEC(i).year = lidar_surface_FAC(i).surf - lidar_surface_FAC(i+1).surf;
    radar_SEC(i).year = radar_surface_FAC(i).surf - radar_surface_FAC(i+1).surf;
    if i == i(end)
        break
    end
end

% Calculate ice divergence

% Calculate melt using Shean et al., 2019 Method (includes ice divergence)
%basal_melt_rate_13_14 = (((lidar_SEC(1).year/1) + lidar_surface_FAC(1).surf*(divergence(1).year))*water_ice_density_ratio) + SMB;

%% 
figure(1)
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, lidar_surface_FAC(1).surf(1:end));
hold on
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, radar_surface_FAC(1).surf(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, lidar_surface_FAC(2).surf(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, radar_surface_FAC(2).surf(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, lidar_surface_FAC(3).surf(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, radar_surface_FAC(3).surf(1:end));
xlabel('along track distance (km)');
ylabel('elevation above WGS84 0-surface (m)');
legend show

figure(2)
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(1).layers(1).layer_elev(1:end));
hold on
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(1).layers(3).layer_elev(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(2).layers(1).layer_elev(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(2).layers(3).layer_elev(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(3).layers(1).layer_elev(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(3).layers(3).layer_elev(1:end));
xlabel('along track distance (km)');
ylabel('elevation above WGS84 0-surface (m)');
legend show

figure(4)
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(1).layers(2).layer_elev(1:end));
hold on
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(2).layers(2).layer_elev(1:end));
plot(pass(baseline_master_idx).along_track/1e3 + pass(1).vel/1e3, pass(3).layers(2).layer_elev(1:end));
xlabel('along track distance (km)');
ylabel('elevation above WGS84 0-surface (m)');
legend show

%% Lebrocq et al., 2013 method of calculating predicted basal melt
% [(thickness*velocity) - (thickness*velocity@GL)]/distanceGL - SMB = B
SMB = zeros(size(thickness_lidar_FAC(1).layer));
clipped_thickness_lidar_FAC = thickness_lidar_FAC(1).layer(404:end);
clipped_velocity = pass_velocity(1).vel(404:end);
clipped_distance = pass(1).along_track(404:end);

for j = 1:numel(thickness_lidar_FAC)
    clipped_H_lidar_FAC(j).layer = thickness_lidar_FAC(j).layer(404:end);
    clipped_vel(j).vel = pass_velocity(j).vel(404:end);
    for i = 1:length(clipped_vel(1).vel)
        AT_point(j).layer(i) = clipped_H_lidar_FAC(j).layer(i)*clipped_vel(j).vel(i);
        GL_point(j).layer(i) = clipped_H_lidar_FAC(j).layer(1)*clipped_vel(j).vel(1);
        distance(i) = clipped_distance(i);
        clipped_melt(j).layer(i) = (AT_point(j).layer(i) - GL_point(j).layer(i))/distance(i) - SMB(i);
    end
end


%% clipped loop
for i = 1:length(clipped_thickness_lidar_FAC)
    AT_point(i)= clipped_thickness_lidar_FAC(i)*clipped_velocity(i);
    GL_point = clipped_thickness_lidar_FAC(1)*clipped_velocity(1);
    distance(i) = clipped_distance(i);
    clipped_basal_melt(i) = (AT_point(i) - GL_point)/distance(i) - SMB(i);
end
clipped_mean_melt = mean(clipped_basal_melt);

figure(9)
plot(clipped_distance/1e3 + clipped_velocity/1e3, clipped_basal_melt);
hold on
%plot(pass(baseline_master_idx).along_track(404:end)/1e3 + pass(1).vel(404:end)/1e3, pass(1).layers(2).layer_elev(404:end));
xlabel('along track distance (km)');
ylabel('basal melt (m/yr)');
legend show
%% original loop
for i = 1:length(thickness_lidar_FAC(1).layer)
    AT_point(i) = (thickness_lidar_FAC(1).layer(i)*pass_velocity(1).vel(i));
    GL_point = (thickness_lidar_FAC(1).layer(404)*pass_velocity(1).vel(404));
    distance(i) = pass(1).along_track(i);
    basal_melt(i) = (AT_point(i) - GL_point)/distance(i) - SMB(i);
end

figure(10)
plot(pass(baseline_master_idx).along_track(404:end)/1e3 + pass(1).vel(404:end)/1e3, basal_melt(404:end));
hold on
plot(pass(baseline_master_idx).along_track(404:end)/1e3 + pass(1).vel(404:end)/1e3, pass(1).layers(2).layer_elev(404:end));
xlabel('along track distance (km)');
ylabel('basal melt (m/yr)');
legend show

%% Load spreadsheet of original data and x/y velocity components
    file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\XY_velocity_spreadsheets\');
    myfiles = pwd;
    % Define years 
    
    Velocity_CSV = dir(fullfile(myfiles,'P2*')); 
    for i = 1:numel(Velocity_CSV)
        Import_Profile = fullfile(myfiles,Velocity_CSV(i).name);
        XY_data(i).data = readmatrix(Import_Profile);
        XY_data(i).Lons = XY_data(i).data(:,2).';
        XY_data(i).Lats = XY_data(i).data(:,3).';
        XY_data(i).Thickness = XY_data(i).data(:,4).';
        XY_data(i).Surf = XY_data(i).data(:,5).';
        XY_data(i).Lidar = XY_data(i).data(:,6).';
        XY_data(i).Depth = XY_data(i).data(:,7).';
        XY_data(i).XVel_14_15 = XY_data(i).data(:,8).';
        XY_data(i).XVel_15_16 = XY_data(i).data(:,9).';
        XY_data(i).XVel_16_17 = XY_data(i).data(:,10).';
        XY_data(i).YVel_14_15 = XY_data(i).data(:,11).';
        XY_data(i).YVel_15_16 = XY_data(i).data(:,12).';
        XY_data(i).YVel_16_17 = XY_data(i).data(:,13).';        
    end

%     DH_dt = (AT_data.interp_data.P14_thickness - AT_data.interp_data.P17_thickness)/3 ;
%     H_mean = (AT_data.interp_data.P14_thickness + AT_data.interp_data.P17_thickness)/2 ;
%     XY_data(1).XVel_mean = (XY_data(1).XVel_14_15(1:end) + XY_data(1).XVel_15_16(1:end) + XY_data(1).XVel_16_17(1:end))/3 ;
%     mean_xvel = mean(XY_data(1).XVel_mean(1:end));
%     mean_yvel = mean(XY_data(1).YVel_mean(1:end));
%     XY_data(1).YVel_mean = (XY_data(1).YVel_14_15(1:end) + XY_data(1).YVel_15_16(1:end) + XY_data(1).YVel_16_17(1:end))/3 ;
%     Ice_Div = divergence(XY_data(1).XVel_mean, XY_data(1).YVel_mean);
%     b = (DH_dt) - H_mean*(div*u_mean) + SMB_mean ;
%% Obtain Arrays of Position and Velocity, obtain difference between points (Langrangian approach)
% obtain mid point for both position x, velocity and thickness
% save arrays per year and then take average for du/dx + dv/dy
% input into dH/Dt equation

for j = 1:numel(pass)
    % pull along track and subtract second element from first
    along_track = pass(j).along_track;
    pass(j).point_dist = along_track(2:end) - along_track(1:end-1);
    
    % pull velocity magnitude and subtract second element from first (km)
    vel_magnitude = pass(j).vel;
    pass(j).vel_mag_difference = (vel_magnitude(2:end)/1000) - (vel_magnitude(1:end-1)/1000);
    
    % pull x-velocity magnitude and subtract second element from first (km)
    vel_magnitude_x = pass(j).vel;
    pass(j).vel_mag_difference_x = (vel_magnitude_x(2:end)/1000) - (vel_magnitude_x(1:end-1)/1000);
    
    % pull y-velocity magnitude and subtract second element from first (km)
    vel_magnitude_y = pass(j).vel;
    pass(j).vel_mag_difference_y = (vel_magnitude_y(2:end)/1000) - (vel_magnitude_y(1:end-1)/1000);
end

for j = 1:numel(pass)    
    point_dist(j).x = pass(j).along_track;
    vel_dif(j).vel = pass(j).vel;
    for i = 1:length(pass(1).along_track)
        point_distx(j).pt(i) = point_dist(j).x(i+1) - point_dist(j).x(i); %- pass(1).along_track(i);
        vel_difv(j).vel(i) = vel_dif(j).vel(i+1) - vel_dif(j).vel(i); %- pass(j).vel(i);
        if i == i(end)
            break
        end
    end
end

% describe to group results and workflow 
% draw out approach and 
% describe and write out why i do/dont need partial/full derivative 
% scematic for work and workflow

% works
for j = 1:numel(thickness_lidar_FAC)
    clipped_H_lidar_FAC(j).layer = thickness_lidar_FAC(j).layer(404:end);
    clipped_vel(j).vel = pass_velocity(j).vel(404:end);
    for i = 1:length(clipped_vel(1).vel)
        AT_point(j).layer(i) = clipped_H_lidar_FAC(j).layer(i)*clipped_vel(j).vel(i);
        GL_point(j).layer(i) = clipped_H_lidar_FAC(j).layer(1)*clipped_vel(j).vel(1);
        distance(i) = clipped_distance(i);
        clipped_melt(j).layer(i) = (AT_point(j).layer(i) - GL_point(j).layer(i))/distance(i) - SMB(i);
    end
end

for i = 1:length(pass(1).along_track)
    point_dist(i).pt = pass(i).along_track - pass(i+1).along_track;
    vel_dif(i).vel = pass(i).vel - pass(i+1).vel;
    if i == i(end)
        break  
    end
end 

