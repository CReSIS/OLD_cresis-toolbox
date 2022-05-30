% Points array for each year
% NOTE: WHEN PICKING POINTS PICK 1 LESS THAN THE LAST DIGIT OFF GRAPH
% TO ACCOUNT FOR ROUNDING ERRORS WHEN PLOTTING
CP.Points.Y13 = [11.2149, 11.7464 13.5059 13.9425 15.0454 16.6193 18.6578];

CP.Points.Y14 = [11.169 11.6007 13.5451 14.008 15.1224 16.7912 18.7925];

CP.Points.Y17 = [11.1301 11.4901 13.5751 13.965 15.1201 16.7701 18.7501];

%% SHORTEST POINTS 
% Left and right position denoted in binary format by 0 = left, 1 = right
CP.Points_short.Y13 = [11.2702 11.9471 13.5898 14.105 14.9274 16.9093 18.2647];
%    R  R  R  R  L  R  L
CP.Points_short.Y14 = [11.2384 12.0153 13.7574 14.2189 15.0273 16.2501 18.8958];
%    R  R  R  R  L  L  R
CP.Points_short.Y17 = [11.1601 11.2501 13.5301 13.6051 15.2101 17.2052 18.1651];
%    R  L  L  L  R  R  L
% uneven crevasse 15.2101, 17.2052
%% LONGEST POINTS
CP.Points_long.Y13 = [11.0745 11.46 13.3666 13.8265 15.2151 16.2818 18.7037];
%    L  L  L  L  R  L  R
% uneven 18.7037 
CP.Points_long.Y14 = [11.0651 11.3397 13.3874 13.8474 15.3364 17.0254 18.3962];
%    L  L  L  L  R  R  L
CP.Points_long.Y17 = [11.0851 12.2251 13.6051 14.4902 15.0451 16.125 18.87];
%    L  R  R  R  L  L  R
%% WIDTHS
CP.Left.Y13 = [11.0745 11.46 13.3666 13.8265 14.9274 16.2818 18.2647];
CP.Right.Y13 = [11.2702 11.9471 13.5898 14.105 15.2151 16.9093 18.7037];

CP.Left.Y14 = [11.0651 11.3397 13.3874 13.8474 15.0273 16.2501 18.3962];
CP.Right.Y14 = [11.2384 12.0153 13.7574 14.2189 15.3364 17.0254 18.8958];

CP.Left.Y17 = [11.0851 11.2501 13.5301 13.6051 15.0451 16.125 18.1651];
CP.Right.Y17 = [11.1601 12.2251 13.6051 14.4902 15.2101 17.2052 18.87];

for i = 1:length(CP.Points.Y13)
    CP.width.Y13(i) = CP.Right.Y13(i) - CP.Left.Y13(i);
    CP.width.Y14(i) = CP.Right.Y14(i) - CP.Left.Y14(i);
    CP.width.Y17(i) = CP.Right.Y17(i) - CP.Left.Y17(i);
end

%% SHORT AND LONG COMPONENT ELEVATION LOOPS
% SHORT
for i = 1:length(CP.Points.Y13)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P13/1e3 >= CP.Points_short.Y13(i));
    CP.Z_short.Y13(i) = AT_data.elev_End_Clip.P2013(CP.AT_start_short(1));
    CP.AT_short.Y13(i) = AT_data.elev_End_Clip_PASS.P2013(CP.AT_start_short(1))/1e3;
end
for i = 1:length(CP.Points.Y14)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P14/1e3 >= CP.Points_short.Y14(i));
    CP.Z_short.Y14(i) = AT_data.elev_End_Clip.P2014(CP.AT_start_short(1));
    CP.AT_short.Y14(i) = AT_data.elev_End_Clip_PASS.P2014(CP.AT_start_short(1))/1e3;
end
for i = 1:length(CP.Points.Y17)
    CP.AT_start_short = find(AT_data.Btrack_End_Clip.P17/1e3 >= CP.Points_short.Y17(i));
    CP.Z_short.Y17(i) = AT_data.elev_End_Clip.P2017(CP.AT_start_short(1));
    CP.AT_short.Y17(i) = AT_data.elev_End_Clip_PASS.P2017(CP.AT_start_short(1))/1e3;
end
% LONG
for i = 1:length(CP.Points.Y13)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P13/1e3 >= CP.Points_long.Y13(i));
    CP.Z_long.Y13(i) = AT_data.elev_End_Clip.P2013(CP.AT_start_long(1));
    CP.AT_long.Y13(i) = AT_data.elev_End_Clip_PASS.P2013(CP.AT_start_long(1))/1e3;
end
for i = 1:length(CP.Points.Y14)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P14/1e3 >= CP.Points_long.Y14(i));
    CP.Z_long.Y14(i) = AT_data.elev_End_Clip.P2014(CP.AT_start_long(1));
    CP.AT_long.Y14(i) = AT_data.elev_End_Clip_PASS.P2014(CP.AT_start_long(1))/1e3;
end
for i = 1:length(CP.Points.Y17)
    CP.AT_start_long = find(AT_data.Btrack_End_Clip.P17/1e3 >= CP.Points_long.Y17(i));
    CP.Z_long.Y17(i) = AT_data.elev_End_Clip.P2017(CP.AT_start_long(1));
    CP.AT_long.Y17(i) = AT_data.elev_End_Clip_PASS.P2017(CP.AT_start_long(1))/1e3;
end
% MEAN (SJORT+LONG/2)
for i = 1:length(CP.Points.Y13)
    CP.Z_mean.Y13(i) = (CP.Z_short.Y13(i) + CP.Z_long.Y13(i))/2 ;
    CP.Z_mean.Y14(i) = (CP.Z_short.Y14(i) + CP.Z_long.Y14(i))/2 ;
    CP.Z_mean.Y17(i) = (CP.Z_short.Y17(i) + CP.Z_long.Y17(i))/2 ;
    CP.AT_mean.Y13(i) = (CP.AT_short.Y13(i) + CP.AT_long.Y13(i))/2 ;
    CP.AT_mean.Y14(i) = (CP.AT_short.Y14(i) + CP.AT_long.Y14(i))/2 ;
    CP.AT_mean.Y17(i) = (CP.AT_short.Y17(i) + CP.AT_long.Y17(i))/2 ;
end

%% Loops for getting each data component (lat,lon, z, AT)
for i = 1:length(CP.Points.Y13)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P13/1e3 >= CP.Points.Y13(i));
    CP.lat.Y13(i) = AT_data.elev_End_Clip_LAT.P2013(CP.AT_start(1)); 
    CP.lon.Y13(i) = AT_data.elev_End_Clip_LON.P2013(CP.AT_start(1));
    CP.Z.Y13(i) = AT_data.elev_End_Clip.P2013(CP.AT_start(1));
    CP.AT.Y13(i) = AT_data.elev_End_Clip_PASS.P2013(CP.AT_start(1))/1e3;
end

for i = 1:length(CP.Points.Y14)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P14/1e3 >= CP.Points.Y14(i));
    CP.lat.Y14(i) = AT_data.elev_End_Clip_LAT.P2014(CP.AT_start(1)); 
    CP.lon.Y14(i) = AT_data.elev_End_Clip_LON.P2014(CP.AT_start(1));
    CP.Z.Y14(i) = AT_data.elev_End_Clip.P2014(CP.AT_start(1));
    CP.AT.Y14(i) = AT_data.elev_End_Clip_PASS.P2014(CP.AT_start(1))/1e3;
end

for i = 1:length(CP.Points.Y17)
    CP.AT_start = find(AT_data.Btrack_End_Clip.P17/1e3 >= CP.Points.Y17(i));
    CP.lat.Y17(i) = AT_data.elev_End_Clip_LAT.P2017(CP.AT_start(1)); 
    CP.lon.Y17(i) = AT_data.elev_End_Clip_LON.P2017(CP.AT_start(1));
    CP.Z.Y17(i) = AT_data.elev_End_Clip.P2017(CP.AT_start(1));
    CP.AT.Y17(i) = AT_data.elev_End_Clip_PASS.P2017(CP.AT_start(1))/1e3;
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
CP.export.P13 = cat(2, CP.lon.Y13, CP.lat.Y13, CP.Z.Y13, CP.Z_short.Y13,  CP.Z_long.Y13, CP.Z_mean.Y13, CP.AT.Y13, CP.width.Y13); 
CP.export.P14 = cat(2, CP.lon.Y14, CP.lat.Y14, CP.Z.Y14, CP.Z_short.Y14,  CP.Z_long.Y14, CP.Z_mean.Y14, CP.AT.Y14, CP.width.Y14);  
CP.export.P17 = cat(2, CP.lon.Y17, CP.lat.Y17, CP.Z.Y17, CP.Z_short.Y17,  CP.Z_long.Y17, CP.Z_mean.Y17, CP.AT.Y17, CP.width.Y17); 

% Change directory to save new transposed files to
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Crevasse_points_files';

%% NEW HEADER METHOD: Header Array of strings and vertically concatenate to data 
cheader = {'Lons', 'Lats', 'Depth','min_depth', 'Max_depth','Mean_depth', 'AT Position', 'width'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

%write header to file 2013
fid = fopen('P2_Crevasse_points_13.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P2_Crevasse_points_13.csv', CP.export.P13, '-append');

%write header to file 2014
fid = fopen('P2_Crevasse_points_14.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P2_Crevasse_points_14.csv', CP.export.P14, '-append');

%write header to file 2017
fid = fopen('P2_Crevasse_points_17.csv','w');
fprintf(fid,'%s\n',textHeader);
fclose(fid);
dlmwrite('P2_Crevasse_points_17.csv', CP.export.P17, '-append');

%% Sectoin 2 - Basal Crevasse height Change Plotting
    % Change working directory to the CSV_export_files path
    file_path = cd('C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Crevasse_points_files\');
    myfiles = pwd;
    % Define years 
    
    Radar_data = dir(fullfile(myfiles,'P2*.csv')); 
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
%%
%for k = 1:numel(Radar_data)
    for i = 1:length(Radar_data(i).Depth)
        DH_points(i).points = [Radar_data(1).Depth(i) Radar_data(2).Depth(i) Radar_data(3).Depth(i)];
        DH_pass(i).points = [Radar_data(1).Pass(i) Radar_data(2).Pass(i) Radar_data(3).Pass(i)];
        DH_Min_points(i).points = [Radar_data(1).elev_min(i) Radar_data(2).elev_min(i) Radar_data(3).elev_min(i)];
        DH_Max_points(i).points = [Radar_data(1).elev_max(i) Radar_data(2).elev_max(i) Radar_data(3).elev_max(i)];
        DH_Mean_points(i).points = [Radar_data(1).elev_mean(i) Radar_data(2).elev_mean(i) Radar_data(3).elev_mean(i)];
        DH_width_points(i).points = [Radar_data(1).width(i) Radar_data(2).width(i) Radar_data(3).width(i)];
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
       sizes = abs(((DH_width_points(i).points - total_mean)/total_std)*100);
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