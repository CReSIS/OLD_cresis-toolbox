if isfield(pass,'velx')
    for i = 1:numel(pass)
        plot_data(i).velx = pass(i).velx;
        if i == 5
            plot_data(5).velx = zeros(size(pass(4).velx));
        end
    end
end

if isfield(pass,'vely')
    for i = 1:numel(pass)
        plot_data(i).vely = pass(i).vely;
        if i == 5
            plot_data(5).vely = zeros(size(pass(4).vely));
        end
    end
end

%%
for i = 1:numel(pass)
    plot_data(i).surf = pass(i).layers(1).layer_elev.';
    plot_data(i).bed =  pass(i).layers(2).layer_elev.';
    plot_data(i).lon = interp1(pass(i).lon, pass(i).lon, pass(5).lon);
    plot_data(i).lat = interp1(pass(i).lat, pass(i).lat, pass(5).lat);
end

for i = 1:numel(plot_data)
    [plot_data(i).x , plot_data(i).y] = projfwd(vel_proj, plot_data(i).lat, plot_data(i).lon);
    plot_data(i).x_cor = plot_data(i).x + plot_data(i).velx ; 
    plot_data(i).y_cor = plot_data(i).y + plot_data(i).vely ; 
end 
%% Run for x
plot_data(1).velx = pass(1).velx;
plot_data(2).velx = pass(2).velx;
plot_data(3).velx = pass(3).velx;
plot_data(4).velx = pass(4).velx;
plot_data(5).velx = pass(5).velx;
%% run for y
plot_data(1).vely = pass(1).vely;
plot_data(2).vely = pass(2).vely;
plot_data(3).vely = pass(3).vely;
plot_data(4).vely = pass(4).vely;
plot_data(5).vely = pass(5).vely;

%% Y Velcotity plot
figure(111)
plot(pass(5).along_track/1e3, plot_data(1).vely);
hold on
plot(pass(5).along_track/1e3, plot_data(2).vely);
plot(pass(5).along_track/1e3, plot_data(3).vely);
plot(pass(5).along_track/1e3, plot_data(4).vely);
%plot(pass(5).along_track/1e3, plot_data(5).vely);
legend on
ylim([-14e9 2e9]);
xlabel('along track distance (km)');
ylabel('velocity correction (m/year');
title('Y-Velocity Component');


figure(112)
plot(pass(5).along_track/1e3, plot_data(1).velx);
hold on
plot(pass(5).along_track/1e3, plot_data(2).velx);
plot(pass(5).along_track/1e3, plot_data(3).velx);
plot(pass(5).along_track/1e3, plot_data(4).velx);
%plot(pass(5).along_track/1e3, plot_data(5).velx);
legend on
ylim([-14e9 2e9]);
xlabel('along track distance (km)');
ylabel('velocity correction (m/year');
title('X-Velocity Component');
%%
for i = 1:numel(plot_data)
    [plot_data(i).vel_lon, plot_data(i).vel_lat] = projinv(vel_proj,plot_data(i).x_cor, plot_data(i).y_cor);
end

figure(22)
h1 = plot3(plot_data(1).vel_lon, plot_data(1).vel_lat, plot_data(1).bed);
hold on
h2 = plot3(plot_data(2).vel_lon, plot_data(2).vel_lat, plot_data(2).bed);
h3 = plot3(plot_data(3).vel_lon, plot_data(3).vel_lat, plot_data(3).bed);
h4 = plot3(plot_data(4).vel_lon, plot_data(4).vel_lat, plot_data(4).bed);
h5 = plot3(plot_data(5).vel_lon, plot_data(5).vel_lat, plot_data(5).bed);
zlabel('elevation relative to WGS84 SSE')
xlabel('lon')
ylabel('lat')
xlim([80 81]);
ylim([-120 -117]);
%%
figure(21)
h1 = plot3(plot_data(1).x, plot_data(1).y, plot_data(1).bed);
hold on
h2 = plot3(plot_data(2).x, plot_data(2).y, plot_data(2).bed);
h3 = plot3(plot_data(3).x, plot_data(3).y, plot_data(3).bed);
h4 = plot3(plot_data(4).x, plot_data(4).y, plot_data(4).bed);
h5 = plot3(plot_data(5).x, plot_data(5).y, plot_data(5).bed);
zlabel('elevation relative to WGS84 SSE')
xlabel('x projected (m)')
ylabel('y projected (m)')
%% Section 1 - 3D plotting and loading of echogram profiles
    % Change working directory to the CSV_export_files path
    file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\CSV_export_files\');
    myfiles = pwd;

    % Read in all files in folder 
    % split radar data file into fields Lon, Lat, Surf, Bed and take transpose
    Radar_data = dir(fullfile(myfiles,'*.csv')); 
    for i = 1:numel(Radar_data)
        Import_Profile = fullfile(myfiles,Radar_data(i).name);
        Radar_data(i).data = readmatrix(Import_Profile);
        Radar_data(i).Lon = Radar_data(i).data(:,1).';
        Radar_data(i).Lat = Radar_data(i).data(:,2).';
        Radar_data(i).Surf = Radar_data(i).data(:,3).';
        Radar_data(i).Bed = Radar_data(i).data(:,4).';
    end

    % plotting 
    figure(20)
    h1 = plot3(Radar_data(1).Lon,Radar_data(1).Lat,Radar_data(1).Bed, 'color','#0072BD');
    hold on
    h2 = plot3(Radar_data(2).Lon,Radar_data(2).Lat,Radar_data(2).Bed, 'color','#D95319');
    h3 = plot3(Radar_data(3).Lon,Radar_data(3).Lat,Radar_data(3).Bed, 'color','#EDB120');
    h4 = plot3(Radar_data(4).Lon,Radar_data(4).Lat,Radar_data(4).Bed, 'color','#77AC30');
    h5 = plot3(Radar_data(5).Lon,Radar_data(5).Lat,Radar_data(5).Bed, 'color','r');
    h6 = plot3(Radar_data(1).Lon,Radar_data(1).Lat,Radar_data(1).Surf, 'color','#0072BD');
    h7 = plot3(Radar_data(2).Lon,Radar_data(2).Lat,Radar_data(2).Surf, 'color','#D95319');
    h8 = plot3(Radar_data(3).Lon,Radar_data(3).Lat,Radar_data(3).Surf, 'color','#EDB120');
    h9 = plot3(Radar_data(4).Lon,Radar_data(4).Lat,Radar_data(4).Surf, 'color','#77AC30');
    h10 = plot3(Radar_data(5).Lon,Radar_data(5).Lat,Radar_data(5).Surf, 'color','r');
    
%% Sectoin 2 - Basal Crevasse height Change Plotting
    % Change working directory to the CSV_export_files path
    file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Crevasse_points_files\');
    myfiles = pwd;
    % Define years 
    date_array = [10 11 13 14];
    Radar_data = dir(fullfile(myfiles,'P4*.csv')); 
    for i = 1:numel(Radar_data)
        Import_Profile = fullfile(myfiles,Radar_data(i).name);
        Radar_data(i).data = readmatrix(Import_Profile);
        Radar_data(i).Lon = Radar_data(i).data(:,1).';
        Radar_data(i).Lat = Radar_data(i).data(:,2).';
        Radar_data(i).Depth = Radar_data(i).data(:,3).';
        % Second loop takes date_array and makes field that is the length 
        % of the other fields and adds the column of single date values
        for k = 1:numel(date_array)
            Radar_data(k).date = zeros(length(Radar_data(i).Depth),1) + date_array(k);
        end
    end
    
    % loop for making DT_plot field in each column
    for i = 1:numel(Radar_data) 
        Radar_data(i).DT_plot = [2010.5 2012 2013.5];
    end 
    % loop for obtaining height difference between years
    for i = 1:numel(Radar_data)
        for k = 1:length(Radar_data(i).Depth)
        Radar_data(i).DH(k) = (Radar_data(i).Depth(k) - Radar_data(i+1).Depth(k));
        Radar_data(i).DT(k) = (Radar_data(i).date(k) - Radar_data(i+1).date(k)); 
        end 
        if i == i(end)
            break
        end
    end
    
    %% Crevasse height change loop, produces array for plotting 
    for k = 1:length(Radar_data(1).Depth)
        Crevasse_id(k).Profile = -1*[Radar_data(1).DH(k) Radar_data(2).DH(k) Radar_data(3).DH(k)];
    end
    
   %% Plotting
   figure(100)
   hold on
   for k = 1:length(Crevasse_id)
       text_1 = ['Crevasse ', num2str(k)];
       plot(Radar_data(1).DT_plot, Crevasse_id(k).Profile, 'DisplayName', text_1);
       xlim([2010 2014]);
       ylim([-50 100]);
       title('DH/DT crevasse apex change');
       xlabel('Date (years)');
       ylabel('DH (m)');
       legend show
       grid on
   end
