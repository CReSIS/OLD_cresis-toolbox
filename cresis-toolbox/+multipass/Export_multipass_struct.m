% load pass struct data into new struct for manipulation
clear pos export_data i j

for i = 1:numel(pass)
    pos(i).x_pos = pass(i).proj_x;
    pos(i).y_pos = pass(i).proj_y;
    pos(i).lat = pass(i).lat;
    pos(i).lon = pass(i).lon;
    pos(i).surf = pass(i).layers(1).layer_elev;
    pos(i).bed = pass(i).layers(2).layer_elev;
    pos(i).thickness = pass(i).layers(1).layer_elev - pass(i).layers(2).layer_elev;
    pos(i).along_track = pass(i).along_track;
    pos(i).velocity_mag = pass(i).vel;
    pos(i).date = pass(i).param_pass.day_seg;
    pos(i).longest_array = length(pass(i).lat);
end

%% Interpolate fields with non-matching array sizes using largest array
% To change master array ID with master pass, change last input on the interp1 function. 
for i = 1:numel(pass)
    pos(i).x_pos = interp1(pos(i).x_pos, pos(i).x_pos, pos(3).x_pos, 'makima', 'extrap');
    pos(i).y_pos = interp1(pos(i).y_pos, pos(i).y_pos, pos(3).y_pos, 'makima', 'extrap');
    pos(i).lat = interp1(pos(i).lat, pos(i).lat, pos(3).lat, 'makima', 'extrap');
    pos(i).lon = interp1(pos(i).lon, pos(i).lon, pos(3).lon, 'makima', 'extrap');
    pos(i).along_track = interp1(pos(i).along_track, pos(i).along_track, pos(3).along_track, 'makima', 'extrap');
%     pos(i).surf = interp1(pos(i).surf, pos(i).surf, pos(2).surf, 'makima','extrap');
%     pos(i).bed = interp1(pos(i).bed, pos(i).bed, pos(2).bed, 'makima','extrap');
%     pos(i).thickness = interp1(pos(i).thickness, pos(2).thickness, pos(2).thickness, 'makima','extrap');
%     pos(i).velocity_mag = interp1(pos(i).velocity_mag, pos(i).velocity_mag, pos(2).velocity_mag, 'makima','extrap');
end

% Interpolate using linspace setting  
% for i = 1:numel(pass)
%     pos(i).x_pos = interp1(pos(i).x_pos, pos(i).x_pos, linspace(pos(i).x_pos(1), pos(i).x_pos(end),length(pos(1).x_pos)), 'makima', 'extrap');
%     pos(i).y_pos = interp1(pos(i).y_pos, pos(i).y_pos, linspace(pos(i).y_pos(1), pos(i).y_pos(end),length(pos(1).y_pos)), 'makima', 'extrap');
%     pos(i).lat = interp1(pos(i).lat, pos(i).lat, linspace(pos(i).lat(1), pos(i).lat(end),length(pos(1).lat)), 'makima', 'extrap');
%     pos(i).lon = interp1(pos(i).lon, pos(i).lon, linspace(pos(i).lon(1), pos(i).lon(end),length(pos(1).lon)), 'makima', 'extrap');
%     pos(i).along_track = interp1(pos(i).along_track, pos(i).along_track, linspace(pos(i).along_track(1), pos(i).along_track(end),length(pos(1).along_track)), 'makima', 'extrap');
%     pos(i).surf = interp1(pos(i).surf, pos(i).surf, linspace(pos(i).surf(1), pos(i).surf(end),length(pos(1).surf)), 'makima','extrap');
%     pos(i).bed = interp1(pos(i).bed, pos(i).bed, linspace(pos(i).bed(1), pos(i).bed(end),length(pos(1).bed)), 'makima','extrap');
%     pos(i).thickness = interp1(pos(i).thickness, pos(i).thickness, linspace(pos(i).thickness(1), pos(i).thickness(end),length(pos(1).thickness)), 'makima','extrap');
%     pos(i).velocity_mag = interp1(pos(i).velocity_mag, pos(i).velocity_mag, linspace(pos(i).velocity_mag(1), pos(i).velocity_mag(end),length(pos(1).velocity_mag)), 'makima','extrap');
% end

% Concatenate all fields together for each season and save into new export
% data structure with dataframes for each season. Fields surf, bed, and
% thickness are not tranposed since they are flipped compared to the other
% remaining fields
for j = 1:numel(pos)
    export_data(j).data_frame = cat(2, pos(j).x_pos.', pos(j).y_pos.',...
        pos(j).lon.', pos(j).lat.', pos(j).surf, pos(j).bed,...
        pos(j).thickness, pos(j).along_track.', pos(j).velocity_mag.');
end

%% Define Header Array of strings and vertically concatenate to data 
cheader = {'X', 'Y','Lons', 'Lats','Surface', 'Bed', 'Thickness', 'AT Distance', 'Velocity Mag'}; % header
commaHeader = [cheader;repmat({','},1,numel(cheader))];
commaHeader = commaHeader(:)';
textHeader = cell2mat(commaHeader);

% change folder
cd 'C:\Users\c262b531\Documents\scripts\cresis-toolbox\cresis-toolbox\+multipass\Pass_struct_CSV\'

% Change file name per file prior to exporting, comment out excess file
% names if there are fewer years than 5 being exported to CSV

for j = 1:numel(export_data)
    name(j).file = append('79N_EXTEND_pass_E_makima_', string(pos(j).date),'.csv');
    fid = fopen(name(j).file,'w');
    fprintf(fid,'%s\n',textHeader);    
    fclose(fid);
    dlmwrite(name(j).file, export_data(j).data_frame, '-append', 'precision','%.7g');
end

%% Old Hard Code
%write header to file 2013
% fid = fopen('P1_pass_11_E.csv','w');
% fprintf(fid,'%s\n',textHeader);
% fclose(fid);
% dlmwrite('P1_pass_11_E.csv', export_data(1).data_frame, '-append');
% 
% %write header to file 2014
% fid = fopen('P1_pass_14_E.csv','w');
% fprintf(fid,'%s\n',textHeader);
% fclose(fid);
% dlmwrite('P1_pass_14_E.csv', export_data(2).data_frame, '-append');
% 
% %write header to file 2017
% fid = fopen('P1_pass_17_E.csv','w');
% fprintf(fid,'%s\n',textHeader);
% fclose(fid);
% dlmwrite('P1_pass_17_E.csv', export_data(3).data_frame, '-append');
% 
% %write header to file 2017
% fid = fopen('P1_pass_17.csv','w');
% fprintf(fid,'%s\n',textHeader);
% fclose(fid);
% dlmwrite('P1_pass_17.csv', export_data(4).data_frame, '-append');
% 
% % %write header to file 2017
% fid = fopen('P1_pass_18.csv','w');
% fprintf(fid,'%s\n',textHeader);
% fclose(fid);
% dlmwrite('P1_pass_18.csv', export_data(5).data_frame, '-append');
