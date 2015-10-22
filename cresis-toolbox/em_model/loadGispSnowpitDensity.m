function [depth,density] = loadGispSnowpitDensity(files)
% [depth,density] = loadGispSnowpitDensity(files)
%
% Loads all the gisp snow pit data and creates an average
% density profile from this.
%
% files is optional, and allows the user to specify only a subset of the snow pits
%
% depth is in m
% density is in kg/m^3

global gRadar;

if ~exist('files','var')
  files = {fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_13.txt'),
    fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_15.txt'),
    fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_31.txt'),
    fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_37.txt'),
    fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_44.txt'),
    fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_571.txt'),
    fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_57.txt'),
    fullfile(gRadar.data_support_path,'em_model_data','gisp_snowpit_73.txt')};
end

% Loads all the density data into a 2-D matrix
%  each column contains the density data from a single snowpit
for index=1:size(files)
  % fprintf('Loading density file: %s\n',files{index});
  d = load(files{index});
  tmp_depth(1:size(d),index) = d(:,1)/100;
  tmp_dens(1:size(d),index) = d(:,3)/1000;
end

% Create new density profile with 1 cm resolution
depth = [0:0.01:max(max(tmp_depth))];

for index=1:size(tmp_depth,2)
  [max_depth,end_ind] = max(tmp_depth(:,index));
  indexes = find(depth <= max_depth);
  tmp_dens2(indexes,index) = ...
    interp1(tmp_depth(1:end_ind,index),tmp_dens(1:end_ind,index),depth(indexes),'linear','extrap').';
end

for index=1:size(tmp_dens2,1)
  density(index) = mean(tmp_dens2(index,find(tmp_dens2(index,:) ~= 0)));
end

return;

% Example to run
[depth,density] = loadGispSnowpitDensity
plot(depth,density)
xlabel('Depth (m)');
ylabel('Density (kg/m^3)');

