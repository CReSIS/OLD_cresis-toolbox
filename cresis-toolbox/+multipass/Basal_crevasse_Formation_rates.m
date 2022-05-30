%% Basal Crevasse Formation Rate Computer
% This script computes the rate of formation for basal crevasses pulled
% from CReSIS OIB profiles that have been aligned using the
% multipass.multipass framework. Identified crevasses are compiled per
% flight line grouping and pass to produce annual number of new crevasses. 
% We then use the first 2km of MEaSUREs velocity data to calculate the
% number of new crevasses formed between epochs. This script is designed to
% run in individual modules, each correponding to a specific flight line
% grouping. 

% 1) step 1 clip all profiles and corresponding velocities to their ground line locations
% 2) clip all profiles to the youngest grounding line location, remove all data down flow of the ground line location 
% 3) automatically compute, or semi-automatically compute number of new crevasse apexes
% 4) count the number of new apexes and save annual count into new field
% 5) divide number of new crevasse by the time difference between profiles

% code string corresponding to glacier of interest
input_string = 'Petermann_line1';

% Petermann 1 module
if strcmpi(input_string,'Petermann_line1')
    % load input data from script of paths (tidal + system delay time corrected)
    multipass.basal_crevasse_paths
    grounding_line_11 = XY_datas(i).AT_vel_corrected(XY_datas(i).AT_vel_corrected <= ...
        XY_datas(1).AT_vel_corrected(282));
    
    
end

if strcmpi(input_string,'Petermann_line2')
    
end

if strcmpi(input_string,'Petermann_line4')
    
end
%%

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

% Clip to Grounding line
%%
figure(999)
plot(XY_datas(1).AT_dist/1e3, XY_datas(1).Bed);
hold on
plot(XY_datas(2).AT_dist/1e3, XY_datas(2).Bed);
plot(XY_datas(3).AT_dist/1e3, XY_datas(3).Bed);
plot(XY_datas(4).AT_dist/1e3, XY_datas(4).Bed);
xlim([0 60]);
ylim([-600 100]);

figure(998)
plot((XY_datas(1).AT_dist + XY_datas(1).Vel_mag)/1e3, XY_datas(1).Bed);
hold on
plot((XY_datas(2).AT_dist + XY_datas(2).Vel_mag)/1e3, XY_datas(2).Bed);
plot((XY_datas(3).AT_dist + XY_datas(3).Vel_mag)/1e3, XY_datas(3).Bed);
plot((XY_datas(4).AT_dist + XY_datas(4).Vel_mag)/1e3, XY_datas(4).Bed);
xlim([0 60]);
ylim([-600 100]);
legend('2011', '2014', '2015', '2017', '2018', 'Location','southeast');