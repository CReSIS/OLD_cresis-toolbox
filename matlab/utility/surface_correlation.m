function surface_correlation(depthData,accumData)
%
%   surface_correlation(depthData,accumData)
%
% This script reads in individual layer files, both from the accumulation 
% and bed depth data, and correlates the relationship between the values
% that correspond to one another. This script performs the following steps:
%
%   (1) Extracts the bed values (meters) from the selected layer file(s).
%   (2) Extracts the accumulation values (meters) from the selected layer
%       files(s).
%   (3) Gets rid of all NaN values (as these will result in a null
%       coefficient value).
%   (4) Uses the interp1 function to interpolate the accum values on top of
%       the bed values, as the GPS times have different offsets.
%   (5) Plots each set of values
%   (6) Runs the corrcoef function to produce an r-value to quantify the
%       relationship between the two datasets. A value of +1/-1 would
%       indicate perfectly correlated data, and a value of 0 would indicate
%       a perfectly random data. 
%
%   Output
%       Two plots within MATLAB, and the r-value, which appears in the
%       command line.
%% Automated section

% Extract GPS data from accum layer file(s)
surface = {};
surf_gps = {};
for surf_file_idx = 1:length(accumData)
  tmp = load(accumData{surf_file_idx});
  surface{surf_file_idx} = tmp.layerData{1}.value{2}.data;
  surf_gps{surf_file_idx} = tmp.GPS_time;
end
surf_all = cell2mat(surface);
surf_gps_all = cell2mat(surf_gps);

% Extract GPS data from the mcords layer file(s)
depth = {};
dep_gps = {};
for dep_file_idx = 1:length(depthData)
  depth{dep_file_idx} = load(depthData{dep_file_idx});
  dep_gps{dep_file_idx} = depth{dep_file_idx}.GPS_time;
end
dep_gps_all = cell2mat(dep_gps)-1;


% Grab ice bed values from layer files 
bed = depth{1}.layerData{2}.value{2}.data;
indx1 = find(surf_gps_all>dep_gps_all(500),1,'first');
indx2 = find(surf_gps_all<dep_gps_all(1600),1,'last');
dep_gps_all(find(isnan(bed)==1))= [];
bed(find(isnan(bed)==1)) = [];
depth_i = interp1(dep_gps_all,bed,surf_gps_all(indx1:indx2));

subplot(1,2,1);plot(surf_gps_all(indx1:indx2)-surf_gps_all(indx1),surf_all(indx1:indx2)*1e6);axis ij;grid
xlabel('relative gps time(sec)')
ylabel('two-way propagation time(us)');
title('ice surface')
subplot(1,2,2);plot(surf_gps_all(indx1:indx2)-surf_gps_all(indx1),depth_i*1e6);axis ij;grid;
xlabel('relative gps time(sec)')
ylabel('two-way propagation time(us)');
title('ice bed')
R = corrcoef(depth_i,surf_all(indx1:indx2));
disp(R);


