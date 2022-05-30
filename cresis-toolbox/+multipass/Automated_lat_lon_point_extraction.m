% Points array for each year
% NOTE: WHEN PICKING POINTS PICK 1 LESS THAN THE LAST DIGIT OFF GRAPH
% TO ACCOUNT FOR ROUNDING ERRORS WHEN PLOTTING
CP.Points.Y11 = [14.9947 15.7007 17.5761 18.8462 21.4078 22.7891 23.1544...
    24.2345 28.651 36.1104 39.9339 41.3364 42.6618 44.0113 49.5285...
    51.4564 52.8096 54.487 56.5378];

CP.Points.Y14 = [14.9355 15.7716 17.7242 18.9134 21.4272 22.2746 22.4415...
    24.7764 28.2662 36.4379 39.5901 41.0930 42.4679 43.9417 48.8925 ...
    51.3161 52.9193 54.5392 56.5654];

CP.Points.Y15 = [15.0301 16.0365 17.84 18.9260 20.9326 22.1983 22.6077...
    25.104 28.4144 36.4310 39.6212 40.6842 42.4915 43.8703 48.7739 ...
    51.3455 52.854 54.5158 56.537];

CP.Points.Y17 = [15.0599 15.9298 17.6901 18.8128 21.3071 21.8764 22.4195...
    25.4552 28.2531 36.2321 39.6269 40.9183 42.412 43.9808 48.4615 ...
    51.1957 53.0641 54.374 56.4656];

CP.Points.Y18 = [15.1048 16.005 18.165 19.1554 21.3748 21.8764 22.8298 ...
    25.5747 28.3647 36.2403 39.9448 40.9347 42.5396 43.9648 48.405 ...
    51.1054 52.9946 54.4203 56.5493];

%% SHORTEST POINTS 
% Left and right position denoted in binary format by 0 = left, 1 = right
CP.Points_short.Y11 = [15.2292 15.7006 17.3571 19.0493 21.2853 22.9286 22.9286...
    24.3236 28.3220 35.3155 40.0205 41.259 42.9738 44.8427 49.7146 51.6482...
    52.6856 54.4338 56.8320];
% No crevasse for 2011 at 15.7006 
%    R  R  L  R  L  R  L  R  L  L  R  L  R  R  R  R  L  L  R
CP.Points_short.Y14 = [14.7748 15.5424 17.4603 19.0747 21.613 22.3331 22.558...
    24.9065 28.1535 36.5880 39.6752 41.6103 42.4325 44.5938 49.6617 51.378...
    52.8877 54.9795 56.9382];
% 28.1536 might have error in bed picking
%    L  L  L  R  R  R  R  R  L  R  R  R  L  R  R  R  L  R  R
CP.Points_short.Y15 = [15.1958 15.9257 17.7718 18.374 21.2501 22.3828 22.7438...
    25.1890 28.3033 36.3390 39.7190 40.6060 42.882 44.737 49.6458 51.360 ...
    52.8392 54.68 57.055];
%    R  L  L  L  R  R  R  R  L  L  R  L  R  R  R  L  L  R  R
CP.Points_short.Y17 = [14.8970 16.0810 17.5313 19.1612 21.3756 21.6305 22.5428...
    25.3304 27.5677 36.6504 39.7477 40.8207 42.8106 45.1926 49.9680 51.1503...
    53.0036 54.9297 57.1111];
%    L  R  L  R  R  L  R  L  L  R  R  L  R  R  R  L  L  R  R
CP.Points_short.Y18 = [15.3448 16.56 18.04 19.2002 21.4349 21.9301 22.8749 ...
    25.4697 28.1699 36.7349 39.7348 40.6049 43.0945 45.1051 48.8996 51.2104 ...
    52.9198 55.1256 57.0292];
%    R  R  L  R  R  R  R  L  L  R  L  L  R  R  R  R  L  R  R

%% LONGEST POINTS
CP.Points_long.Y11 = [14.7127 15.7006 17.756 18.5687 21.819 22.6972 23.3446...
    23.9600 28.9587 36.6103 39.7914 41.7665 42.4811 43.88 49.281 51.0613...
    53.1404 54.6170 56.2511];
% 36.6103 uses entire width of cavity, 51.0613 uses base of vertical wall
% 56.2511 uses the base of the wall at the end of the shelf
%    L  L  R  L  R  L  R  L  R  R  L  R  L  L  L  L  R  R  L
CP.Points_long.Y14 = [15.1675 15.9272 17.8147 18.617 20.8972 22.0783 22.3331...
    24.6385 28.5971 36.3140 39.5196 40.9722 43.0034 43.7475 48.8013 51.2532...
    52.961 54.1823 56.2644];
% 56.2644 was chosen at base of the left wall, unsure where to pick from
%    R  R  R  L  L  L  L  L  R  L  L  L  R  L  L  L  R  L  L
CP.Points_long.Y15 = [14.979 16.0861 17.9042 19.3435 20.8516 22.0471 22.5074...
    24.9110 28.7385 36.4868 39.4414 40.9977 42.356 43.8178 48.6243 51.454... 
    52.8985 54.1711 56.3573];
%   L  R  R  R  L  L  L  L  R  R  L  R  L  L  L  R  R  L  L
CP.Points_long.Y17 = [15.2378 15.7724 17.9884 18.6706 21.2653 22.0223 22.3928...
    25.5523 28.5826 36.0972 39.5061 41.4445 42.1100 43.8187 48.2181 51.4057...
    53.2588 53.9761 56.1011];
%   R   L  R  L  L  R  L  R  R  L  L  R  L  L  L  R  R  R  R
CP.Points_long.Y18 = [14.8500 15.7348 18.6152 19.0803 21.2696 21.8548 22.7398 ...
    25.6347 28.5446 35.9700 40.1698 41.3097 41.9395 43.8000 48.1948 51.0753...
    53.0396 54.0448 55.9795];
%   L  L  R  L  L  L  L  R  R  L  R  R  L   L  L  L  R  L  L

%% WIDTHS
CP.Left.Y11 = [14.7127 15.7006 17.3571 18.5687 21.2853 22.6972 22.9286 ...
    23.9600 28.3220 35.3155 39.7914 41.259 42.4811 43.88 49.281 51.0613 ...
    52.6856 54.4338 56.2511];
CP.Right.Y11 = [15.2292 15.7006 17.756 19.0493 21.819 22.9286 23.3446 ...
    24.3236 28.9587 36.6103 40.0205 41.7665 42.9738 44.8427 49.7146 ...
    51.6482 53.1404 54.6170 56.8320];

CP.Left.Y14 = [14.7748 15.5424 17.4603 18.617 20.8972 22.0783 22.3331...
    24.6385 28.1535 36.3140 39.5196 40.9722 42.4325 43.7475 48.8013 51.2532...
    52.8877 54.1823 56.2644];
CP.Right.Y14 = [15.1675 15.9272 17.8147 19.0747 21.613 22.3331 22.558...
    24.9065 28.5971 36.5880 39.6752 41.6103 43.0034 44.5938 49.6617 51.378... 
    52.961 54.9795 56.9382];

CP.Left.Y15 = [14.979 15.9257 17.7718 18.374 20.8516 22.0471 22.5074...
    24.9110 28.3033 36.3390 39.4414 40.6060 42.356 43.8178 48.6243...
    51.360 52.8392 54.1711 56.3573];
CP.Right.Y15 = [15.1958 16.0861 17.9042 19.3435 21.2501 22.3828 22.7438...
    25.1890 28.7385 36.4868 39.7190 40.9977 42.882 44.737 49.6458...
    51.454 52.8985 54.68 57.055];

CP.Left.Y17 = [14.8970 15.7724 17.5313 18.6706 21.2653 21.6305 22.3928...
    25.3304 27.5677 36.0972 39.5061 40.8207 42.1100 43.8187 48.2181...
    51.1503 53.0036 53.9761 56.1011];
CP.Right.Y17 = [15.2378 16.0810 17.9884 19.1612 21.3756 22.0223 22.5428...
    25.5523 28.5826 36.6504 39.7477 41.4445 42.8106 45.1926 49.9680 ...
    51.4057 53.2588 54.9297 57.1111];

CP.Left.Y18 = [14.8500 15.7348 18.04 19.0803 21.2696 21.8548 22.7398...
    25.4697 28.1699 35.9700 39.7348 40.6049 41.9395 43.8000 48.1948 ...
    51.0753 52.9198 54.0448 55.9795];
CP.Right.Y18 = [15.3448 16.56 18.6152 19.2002 21.4349 21.9301 22.8749...
    25.6347 28.5446 36.7349 40.1698 41.3097 43.0945 45.1051 48.8996 ...
    51.2104 53.0396 55.1256 57.0292];

for i = 1:length(CP.Points.Y11)
    CP.width.Y11(i) = CP.Right.Y11(i) - CP.Left.Y11(i);
    CP.width.Y14(i) = CP.Right.Y14(i) - CP.Left.Y14(i);
    CP.width.Y15(i) = CP.Right.Y15(i) - CP.Left.Y15(i);
    CP.width.Y17(i) = CP.Right.Y17(i) - CP.Left.Y17(i);
    CP.width.Y18(i) = CP.Right.Y18(i) - CP.Left.Y18(i);
end

%% SHORT AND LONG COMPONENT ELEVATION LOOPS
% SHORT
for i = 1:length(CP.Points.Y11)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P11/1e3 >= CP.Points_short.Y11(i));
    CP.Z_short.Y11(i) = AT_data.elev_End_Clip.P2011(CP.AT_start_short(1));
    CP.AT_short.Y11(i) = AT_data.elev_End_Clip_PASS.P2011(CP.AT_start_short(1))/1e3;
end
for i = 1:length(CP.Points.Y14)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P14/1e3 >= CP.Points_short.Y14(i));
    CP.Z_short.Y14(i) = AT_data.elev_End_Clip.P2014(CP.AT_start_short(1));
    CP.AT_short.Y14(i) = AT_data.elev_End_Clip_PASS.P2014(CP.AT_start_short(1))/1e3;
end
for i = 1:length(CP.Points.Y15)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P15/1e3 >= CP.Points_short.Y15(i));
    CP.Z_short.Y15(i) = AT_data.elev_End_Clip.P2015(CP.AT_start_short(1));
    CP.AT_short.Y15(i) = AT_data.elev_End_Clip_PASS.P2015(CP.AT_start_short(1))/1e3;
end
for i = 1:length(CP.Points.Y17)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P17/1e3 >= CP.Points_short.Y17(i));
    CP.Z_short.Y17(i) = AT_data.elev_End_Clip.P2017(CP.AT_start_short(1));
    CP.AT_short.Y17(i) = AT_data.elev_End_Clip_PASS.P2017(CP.AT_start_short(1))/1e3;
end
for i = 1:length(CP.Points.Y18)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P18/1e3 >= CP.Points_short.Y18(i));
    CP.Z_short.Y18(i) = AT_data.elev_End_Clip.P2018(CP.AT_start_short(1));
    CP.AT_short.Y18(i) = AT_data.elev_End_Clip_PASS.P2018(CP.AT_start_short(1))/1e3;
end
% LONG
for i = 1:length(CP.Points.Y11)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P11/1e3 >= CP.Points_long.Y11(i));
    CP.Z_long.Y11(i) = AT_data.elev_End_Clip.P2011(CP.AT_start_long(1));
    CP.AT_long.Y11(i) = AT_data.elev_End_Clip_PASS.P2011(CP.AT_start_long(1))/1e3;
end
for i = 1:length(CP.Points.Y14)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P14/1e3 >= CP.Points_long.Y14(i));
    CP.Z_long.Y14(i) = AT_data.elev_End_Clip.P2014(CP.AT_start_long(1));
    CP.AT_long.Y14(i) = AT_data.elev_End_Clip_PASS.P2014(CP.AT_start_long(1))/1e3;
end
for i = 1:length(CP.Points.Y15)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P15/1e3 >= CP.Points_long.Y15(i));
    CP.Z_long.Y15(i) = AT_data.elev_End_Clip.P2015(CP.AT_start_long(1));
    CP.AT_long.Y15(i) = AT_data.elev_End_Clip_PASS.P2015(CP.AT_start_long(1))/1e3;
end
for i = 1:length(CP.Points.Y17)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P17/1e3 >= CP.Points_long.Y17(i));
    CP.Z_long.Y17(i) = AT_data.elev_End_Clip.P2017(CP.AT_start_long(1));
    CP.AT_long.Y17(i) = AT_data.elev_End_Clip_PASS.P2017(CP.AT_start_long(1))/1e3;
end
for i = 1:length(CP.Points.Y18)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P18/1e3 >= CP.Points_long.Y18(i));
    CP.Z_long.Y18(i) = AT_data.elev_End_Clip.P2018(CP.AT_start_long(1));
    CP.AT_long.Y18(i) = AT_data.elev_End_Clip_PASS.P2018(CP.AT_start_long(1))/1e3;
end
% MEAN (SJORT+LONG/2)
for i = 1:length(CP.Points.Y18)
    CP.Z_mean.Y11(i) = (CP.Z_short.Y11(i) + CP.Z_long.Y11(i))/2 ;
    CP.Z_mean.Y14(i) = (CP.Z_short.Y14(i) + CP.Z_long.Y14(i))/2 ;
    CP.Z_mean.Y15(i) = (CP.Z_short.Y15(i) + CP.Z_long.Y15(i))/2 ;
    CP.Z_mean.Y17(i) = (CP.Z_short.Y17(i) + CP.Z_long.Y17(i))/2 ;
    CP.Z_mean.Y18(i) = (CP.Z_short.Y18(i) + CP.Z_long.Y18(i))/2 ;
    CP.AT_mean.Y11(i) = (CP.AT_short.Y11(i) + CP.AT_long.Y11(i))/2 ;
    CP.AT_mean.Y14(i) = (CP.AT_short.Y14(i) + CP.AT_long.Y14(i))/2 ;
    CP.AT_mean.Y15(i) = (CP.AT_short.Y15(i) + CP.AT_long.Y15(i))/2 ;
    CP.AT_mean.Y17(i) = (CP.AT_short.Y17(i) + CP.AT_long.Y17(i))/2 ;
    CP.AT_mean.Y18(i) = (CP.AT_short.Y18(i) + CP.AT_long.Y18(i))/2 ;
end

%% Loops for getting each data component (lat,lon, z, AT)
for i = 1:length(CP.Points.Y11)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P11/1e3 >= CP.Points.Y11(i));
    CP.lat.Y11(i) = AT_data.elev_End_Clip_LAT.P2011(CP.AT_start(1)); 
    CP.lon.Y11(i) = AT_data.elev_End_Clip_LON.P2011(CP.AT_start(1));
    CP.Z.Y11(i) = AT_data.elev_End_Clip.P2011(CP.AT_start(1));
    CP.AT.Y11(i) = AT_data.elev_End_Clip_PASS.P2011(CP.AT_start(1))/1e3;
end

for i = 1:length(CP.Points.Y14)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P14/1e3 >= CP.Points.Y14(i));
    CP.lat.Y14(i) = AT_data.elev_End_Clip_LAT.P2014(CP.AT_start(1)); 
    CP.lon.Y14(i) = AT_data.elev_End_Clip_LON.P2014(CP.AT_start(1));
    CP.Z.Y14(i) = AT_data.elev_End_Clip.P2014(CP.AT_start(1));
    CP.AT.Y14(i) = AT_data.elev_End_Clip_PASS.P2014(CP.AT_start(1))/1e3;
end

for i = 1:length(CP.Points.Y15)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P15/1e3 >= CP.Points.Y15(i));
    CP.lat.Y15(i) = AT_data.elev_End_Clip_LAT.P2015(CP.AT_start(1)); 
    CP.lon.Y15(i) = AT_data.elev_End_Clip_LON.P2015(CP.AT_start(1));
    CP.Z.Y15(i) = AT_data.elev_End_Clip.P2015(CP.AT_start(1));
    CP.AT.Y15(i) = AT_data.elev_End_Clip_PASS.P2015(CP.AT_start(1))/1e3;
end

for i = 1:length(CP.Points.Y17)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P17/1e3 >= CP.Points.Y17(i));
    CP.lat.Y17(i) = AT_data.elev_End_Clip_LAT.P2018(CP.AT_start(1)); 
    CP.lon.Y17(i) = AT_data.elev_End_Clip_LON.P2018(CP.AT_start(1));
    CP.Z.Y17(i) = AT_data.elev_End_Clip.P2017(CP.AT_start(1));
    CP.AT.Y17(i) = AT_data.elev_End_Clip_PASS.P2017(CP.AT_start(1))/1e3;
end

for i = 1:length(CP.Points.Y18)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P18/1e3 >= CP.Points.Y18(i));
    CP.lat.Y18(i) = AT_data.elev_End_Clip_LAT.P2018(CP.AT_start(1)); 
    CP.lon.Y18(i) = AT_data.elev_End_Clip_LON.P2018(CP.AT_start(1));
    CP.Z.Y18(i) = AT_data.elev_End_Clip.P2018(CP.AT_start(1));
    CP.AT.Y18(i) = AT_data.elev_End_Clip_PASS.P2018(CP.AT_start(1))/1e3;
end
%% Take transpose of each field
CP.lat = structfun(@transpose,CP.lat,'UniformOutput',false);
CP.lon = structfun(@transpose,CP.lon,'UniformOutput',false);
CP.Z = structfun(@transpose,CP.Z,'UniformOutput',false);
CP.Z_short = structfun(@transpose,CP.Z_short,'UniformOutput',false);
CP.Z_long = structfun(@transpose,CP.Z_long,'UniformOutput',false);
CP.Z_mean = structfun(@transpose,CP.Z_mean,'UniformOutput',false);
CP.AT = structfun(@transpose,CP.AT,'UniformOutput',false);
CP.width = structfun(@transpose,CP.width,'UniformOutput',false);

%% Horizontal concatenation
CP.export.P11 = cat(2, CP.lon.Y11, CP.lat.Y11, CP.Z.Y11, CP.Z_short.Y11,  CP.Z_long.Y11, CP.Z_mean.Y11, CP.AT.Y11, CP.width.Y11); 
CP.export.P14 = cat(2, CP.lon.Y14, CP.lat.Y14, CP.Z.Y14, CP.Z_short.Y14,  CP.Z_long.Y14, CP.Z_mean.Y14, CP.AT.Y14, CP.width.Y14); 
CP.export.P15 = cat(2, CP.lon.Y15, CP.lat.Y15, CP.Z.Y15, CP.Z_short.Y15,  CP.Z_long.Y15, CP.Z_mean.Y15, CP.AT.Y15, CP.width.Y15); 
CP.export.P17 = cat(2, CP.lon.Y17, CP.lat.Y17, CP.Z.Y17, CP.Z_short.Y17,  CP.Z_long.Y17, CP.Z_mean.Y17, CP.AT.Y17, CP.width.Y17); 
CP.export.P18 = cat(2, CP.lon.Y18, CP.lat.Y18, CP.Z.Y18, CP.Z_short.Y18,  CP.Z_long.Y18, CP.Z_mean.Y18, CP.AT.Y18, CP.width.Y18); 

% Change directory to save new transposed files to
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Crevasse_points_files';

%% NEW HEADER METHOD: Header Array of strings and vertically concatenate to data 
cheader = {'Lons', 'Lats', 'Depth','min_depth', 'Max_depth','Mean_depth', 'AT Position', 'width'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

%write header to file 2011
fid = fopen('P1_Crevasse_points_11.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P1_Crevasse_points_11.csv', CP.export.P11, '-append');

%write header to file 2014
fid = fopen('P1_Crevasse_points_14.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P1_Crevasse_points_14.csv', CP.export.P14, '-append');

%write header to file 2015
fid = fopen('P1_Crevasse_points_15.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P1_Crevasse_points_15.csv', CP.export.P15, '-append');

%write header to file 2017
fid = fopen('P1_Crevasse_points_17.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P1_Crevasse_points_17.csv', CP.export.P17, '-append');

%write header to file 2018
fid = fopen('P1_Crevasse_points_18.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P1_Crevasse_points_18.csv', CP.export.P18, '-append');

%% Sectoin 2 - Basal Crevasse height Change Plotting
    % Change working directory to the CSV_export_files path
    file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Crevasse_points_files\');
    myfiles = pwd;
    % Define years 
    
    Radar_data = dir(fullfile(myfiles,'P1*.csv')); 
    for i = 1:numel(Radar_data)
        Import_Profile = fullfile(myfiles,Radar_data(i).name);
        Radar_data(i).data = readmatrix(Import_Profile);
        Radar_data(i).Lon = Radar_data(i).data(:,1).';
        Radar_data(i).Lat = Radar_data(i).data(:,2).';
        Radar_data(i).Depth = Radar_data(i).data(:,3).';
        Radar_data(i).Depth_min = Radar_data(i).data(:,4).';
        Radar_data(i).Depth_max = Radar_data(i).data(:,5).';
        Radar_data(i).Depth_mean = Radar_data(i).data(:,6).';
        Radar_data(i).Pass = Radar_data(i).data(:,7).';
        Radar_data(i).width = Radar_data(i).data(:,8).';
    end

    % loop for obtaining height difference between years
    for i = 1:numel(Radar_data)
        for k = 1:length(Radar_data(i).Depth)
        Radar_data(i).DH = (Radar_data(i).Depth - Radar_data(i+1).Depth);
        end 
        if i == i(end)
            break
        end
    end
    
    % Crevasse height change loop, produces array for plotting 
%     for k = 1:length(Radar_data(1).Depth)
%         Crevasse_id(k).Profile = -1*[Radar_data(1).DH(k) Radar_data(2).DH(k) Radar_data(3).DH(k)];
%     end

%% 
for i = 1:numel(Radar_data)
    Radar_data(i).elev_min = Radar_data(i).Depth - Radar_data(i).Depth_min;
    Radar_data(i).elev_max = Radar_data(i).Depth - Radar_data(i).Depth_max;
    Radar_data(i).elev_mean = (Radar_data(i).elev_max + Radar_data(i).elev_min)/2;
end
%for k = 1:numel(Radar_data)
    for i = 1:length(Radar_data(i).Depth)
        DH_points(i).points = [Radar_data(1).Depth(i) Radar_data(2).Depth(i) Radar_data(3).Depth(i) Radar_data(4).Depth(i) Radar_data(5).Depth(i)];
        DH_pass(i).points = [Radar_data(1).Pass(i) Radar_data(2).Pass(i) Radar_data(3).Pass(i) Radar_data(4).Pass(i) Radar_data(5).Pass(i)];
        DH_Min_points(i).points = [Radar_data(1).elev_min(i) Radar_data(2).elev_min(i) Radar_data(3).elev_min(i) Radar_data(4).elev_min(i) Radar_data(5).elev_min(i)];
        DH_Max_points(i).points = [Radar_data(1).elev_max(i) Radar_data(2).elev_max(i) Radar_data(3).elev_max(i) Radar_data(4).elev_max(i) Radar_data(5).elev_max(i)];
        DH_Mean_points(i).points = [Radar_data(1).elev_mean(i) Radar_data(2).elev_mean(i) Radar_data(3).elev_mean(i) Radar_data(4).elev_mean(i) Radar_data(5).elev_mean(i)];
        DH_width_points(i).points = [Radar_data(1).width(i) Radar_data(2).width(i) Radar_data(3).width(i) Radar_data(4).width(i) Radar_data(5).width(i)];
    end 
        % New min/max/mean crevasse heights
%         DH_Min_points(i).points = [Radar_data(1).elev_min(i) Radar_data(2).elev_min(i) Radar_data(3).elev_min(i) Radar_data(4).elev_min(i) Radar_data(5).elev_min(i)];
%         DH_Max_points(i).points = [Radar_data(1).elev_max(i) Radar_data(2).elev_max(i) Radar_data(3).elev_max(i) Radar_data(4).elev_max(i) Radar_data(5).elev_max(i)];

%         Min_points = [(Radar_data(1).Depth(i) - Radar_data(1).Depth_min(i))  (Radar_data(2).Depth(i) - Radar_data(2).Depth_min(i)) ...
%             (Radar_data(3).Depth(i) - Radar_data(3).Depth_min(i)) (Radar_data(4).Depth(i) - Radar_data(4).Depth_min(i)) ...
%             (Radar_data(5).Depth(i) - Radar_data(5).Depth_min(i))];
%         Max_points = [(Radar_data(1).Depth(i) - Radar_data(1).Depth_max(i))  (Radar_data(2).Depth(i) - Radar_data(2).Depth_max(i)) ...
%             (Radar_data(3).Depth(i) - Radar_data(3).Depth_max(i)) (Radar_data(4).Depth(i) - Radar_data(4).Depth_max(i)) ...
%             (Radar_data(5).Depth(i) - Radar_data(5).Depth_max(i))];
%         Mean_points = [(Radar_data(1).Depth_min(i) + Radar_data(1).Depth_max(i))/2  (Radar_data(2).Depth_min(i) + Radar_data(2).Depth_max(i))/2 ...
%             (Radar_data(3).Depth_min(i) + Radar_data(3).Depth_max(i))/2 (Radar_data(4).Depth_min(i) + Radar_data(4).Depth_max(i))/2 ...
%             (Radar_data(5).Depth_min(i) + Radar_data(5).Depth_max(i))/2];

   %% Plotting ORIGINAL
   % Shows height change from base elevation
   figure(14)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       plot(DH_pass(k).points, DH_points(k).points, '-o', 'DisplayName', text_1);
   end
   title('Crevasse Apex Elevation change with Respect to Distance from Grounding Line (DH/DX)');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Apex Elevation (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   
   figure(15)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       plot(DH_pass(k).points, DH_Min_points(k).points, '-o', 'DisplayName', text_1);
   end
   title('Minimum Crevasse Height with Respect to Distance from Grounding Line');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Minimum Crevasse Height (apex - min) (m)');
   xlim([0 60]);
   ylim([0 140]);
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   
   figure(16)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       plot(DH_pass(k).points, DH_Max_points(k).points, '-o', 'DisplayName', text_1);
   end
   title('Maximum Crevasse Height with Respect to Distance from Grounding Line');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Maximum Crevasse height (apex - max) (m)');
   xlim([0 60]);
   ylim([0 140]);
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
  
   figure(17)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       plot(DH_pass(k).points, DH_Mean_points(k).points, '-o', 'DisplayName', text_1);
   end
   title('Mean Crevasse Height with Respect to Distance from Grounding Line');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Mean Crevasse height (max-min)/2 (m)');
   xlim([0 60]);
   ylim([0 140]);
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
  
   figure(18)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       plot(DH_pass(k).points, DH_width_points(k).points, '-o', 'DisplayName', text_1);
   end
   title('Crevasse widths with Respect to Distance from Grounding Line');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Crevasse Width (km)');
   xlim([0 60]);
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   %% Depicts crevasse apex elevation as function of distance and crevasse 
   % width (circle size). The crevasse widths were normalized using the
   % equation z = (x - u)/sigma, where x is a given element, u is the mean
   % and sigma is the standard deveiation of the dataset. All points were
   % then 1) multipled by 100 and 2) had their absolute values taken for
   % plotting. 
   colors = colormap(flipud(turbo(length(DH_Max_points))));
   total_std = std([DH_width_points.points]);
   total_mean = mean([DH_width_points.points]);
   figure(20)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       %sizes = (DH_width_points(k).points + 1)*100;
       sizes = abs(((DH_width_points(k).points - total_mean)/total_std)*100);
       plot(DH_pass(k).points, DH_points(k).points, '-o', 'MarkerFaceColor',colors(k,:), 'Color',colors(k,:),  'DisplayName', text_1);
       h = plot(DH_pass(k).points, DH_points(k).points, 'o', 'Color','k');
       h.Annotation.LegendInformation.IconDisplayStyle = 'off';
       h2 = scatter(DH_pass(k).points, DH_points(k).points, sizes, 'MarkerFaceColor', colors(k,:), 'MarkerEdgeColor', 'k', 'LineWidth', 1);
       h2.Annotation.LegendInformation.IconDisplayStyle = 'off';
   end
   title('Crevasse Apex Elevation change with Respect to WGS84 0-surface as a function of Distance from Grounding Line and Crevasse Width');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Apex Elevation (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   
   %%
   scaled
   figure(21)
   hold on
   for k = 1:numel(DH_pass)
       text_1 = ['Crevasse ', num2str(k)];
       scatter(DH_pass(k).points, DH_points(k).points, DH_width_points(k).points, 'DisplayName', text_1);
   end
   title('Crevasse Apex Elevation change with Respect to Distance from Grounding Line (DH/DX)');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Apex Elevation (m)');
   lgd = legend;
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   grid on;
   legend show;
   
   %%
   for i= 1:numel(DH_pass)
       sizes2(i).points = (DH_width_points(i).points + 1);
   end
   %%
   figure(23)
   hold on
   colors = colormap(flipud(turbo(length(DH_Max_points))));
   total_std = std([DH_width_points.points]);
   total_mean = mean([DH_width_points.points]);
   
   for i= 1:numel(DH_pass)
       sizes = (DH_width_points(i).points + 1)*100;
       sizes = abs(((DH_width_points(i).points - total_mean)/total_std)*100);
       text_2 = ['Crevasse ', num2str(i)];
       scatter(DH_pass(i).points, DH_points(i).points, sizes2, 'MarkerFaceColor', colors(i,:), 'MarkerEdgeColor', 'k', 'DisplayName', text_2);
   end
   %%
   sizes.points = (DH_width_points(1).points + 1)*100;
   figure(22)
   text = ['Crevasse ', num2str(k)];
   scatter(DH_pass(1).points, DH_points(1).points, sizes, 'DisplayName', text);
   
   
     %% Plotting
   figure(14)
   for k = 1:numel(DH_pass)
       text = ['Crevasse ', num2str(k)];
       plot(DH_pass(k).points, DH_points(k).points, '-o', 'DisplayName', text);
       hold on
   end
   title('Crevasse Apex Elevation change with Respect to Distance from Grounding Line (DH/DX)');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Apex Elevation (m)');
   grid on;
   legend('show', 'Location','southeast')
   %% Plotting
   figure(14)
   hold on
   Legend = cell(numel(DH_pass),1);
   for j = 1:numel(DH_pass)
       Legend(j) = strcat('Crevasse ', num2str(j));
   end
   ldg = legend(Legend);
   lgd.NumColumns = 2;
   lgd.Location = 'best';
   ldg show;
   %%
   for k = 1:numel(DH_pass)
       plot(DH_pass(k).points, DH_points(k).points, '-o');
   end
   hold off
   title('Crevasse Apex Elevation change with Respect to Distance from Grounding Line (DH/DX)');
   xlabel('Distance from Grounding Line (Km)');
   ylabel('Apex Elevation (m)');
   grid on;
  
   
%%
figure(19)
h1 = plot(AT_data.Btrack_End_Clip.P11/1e3, AT_data.elev_End_Clip.P2011);
hold on
h2 = plot(AT_data.Btrack_End_Clip.P14/1e3, AT_data.elev_End_Clip.P2014);
h3 = plot(AT_data.Btrack_End_Clip.P15/1e3, AT_data.elev_End_Clip.P2015);
h4 = plot(AT_data.Btrack_End_Clip.P17/1e3, AT_data.elev_End_Clip.P2017);
h5 = plot(AT_data.Btrack_End_Clip.P18/1e3, AT_data.elev_End_Clip.P2018);
h6 = plot(CP.Points.Y11, CP.Z.Y11, 'o', 'color', 'k', 'MarkerSize', 5);
h7 = plot(CP.Points.Y14, CP.Z.Y14, 'o', 'color', 'k', 'MarkerSize', 5);
h8 = plot(CP.Points.Y15, CP.Z.Y15, 'o', 'color', 'k', 'MarkerSize', 5);
h9 = plot(CP.Points.Y17, CP.Z.Y17, 'o', 'color', 'k', 'MarkerSize', 5);
h10 = plot(CP.Points.Y18, CP.Z.Y18, 'o', 'color', 'k', 'MarkerSize', 5);
%h6 = plot(AT_data.Btrack_End_Clip.P19/1e3, AT_data.elev_End_Clip.P2019);
title('Test Plot 1 - Alignment of Profiles is Correct');
xlabel('Along Track (km)');
ylabel('Elevation (m)');
legend('original 2011','original 2014', 'original 2015', 'original 2017',...
    'original 2018','Location', 'southeast');

 %% ORIGNAL POINT SELECTION, DO NOT CHANGE? DELETE 
% CP.Points.Y11 = [14.9948 15.7007 17.5923 18.8632 21.4078 22.7891 23.1544...
%     24.2345 28.6626 36.0908 39.9339 41.3364 42.6618 44.0113 49.5285...
%     51.4564 52.8096 54.487 56.5378];
% 
% CP.Points.Y14 = [14.9355 15.7716 17.7243 18.9134 21.4272 22.2747 22.4416...
%     24.7925 28.2662 36.4379 39.5901 41.0931 42.4679 43.9417 48.8925 ...
%     51.3318 52.9194 54.5392 56.5654];
% 
% CP.Points.Y15 = [15.0301 15.8398 17.84 18.9261 20.9326 22.1983 22.6077...
%     25.1319 28.4144 36.4310 39.6212 40.6842 42.4915 43.8704 48.7586 ...
%     51.3456 52.869 54.5158 56.537];
% 
% CP.Points.Y17 = [15.0599 16.005 17.6901 18.8128 21.3071 21.8922 22.4195...
%     25.4553 28.2531 36.2321 39.6269 40.9183 42.4276 43.9808 48.4615 ...
%     51.1957 53.0641 54.374 56.4656];
% 
% CP.Points.Y18 = [15.1048 15.8656 18.165 19.1402 21.3749 21.8701 22.6948 ...
%     25.5598 28.3647 36.2403 39.9449 40.9347 42.5396 43.9649 48.405 ...
%     51.1054 52.9947 54.4203 56.5493];