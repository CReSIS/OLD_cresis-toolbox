% script: CSA_surface_extractor_prep
% By John Paden and Mingze Xu, July 2016

% fn_dir: Directory where files are at
fn_dir = './';
fn_dir = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_CSA_music/20140401_03/';

% img_01: right looking (positive theta)
% img_02: nadir looking (theta == 0)
% img_03: left looking (negative theta)
nadir_img = 2;

%% Automated loading section
% =========================================================================

if 0
  % Compile C++ functions
  % If you get a C++11 option error, you may be using pre-G++ 4.7. You can
  % check the g++ version with system('g++ --version');
  % To fix this, add -v option to mex function and look for a line like this:
  %   Options file: ~/.matlab/R2015b/mex_C++_glnxa64.xml
  % Replace -std=c++11 with -std=c++0x (should occur in two places)
  % Reference: http://stackoverflow.com/questions/14674597/cc1plus-error-unrecognized-command-line-option-std-c11-with-g
  mex -largeArrayDims fuse.cpp
  mex -largeArrayDims train_params.cpp
  mex -largeArrayDims detect.cpp
  mex -largeArrayDims extract.cpp
end

% Load Data
mdata = {};
for img=1:3
    fn = fullfile(fn_dir,sprintf('Data_img_%02.0f_20140401_03_044.mat',img));
    mdata{img} = load(fn);
end

% Check for twtt DEM data
if ~isfield(mdata{1},'twtt')
    Nx = length(mdata{1}.GPS_time);
    Nsv = length(mdata{1}.param_combine.array_param.theta);
    twtt = NaN*zeros(Nsv,Nx);
end

% Interpolate images onto a common propagation time axis
for img  = [1 3]
    mdata{img}.Data = interp1(mdata{img}.Time,mdata{img}.Data,mdata{nadir_img}.Time);
    mdata{img}.Topography.img = interp1(mdata{img}.Time,mdata{img}.Topography.img,mdata{nadir_img}.Time);
    mdata{img}.Time = mdata{nadir_img}.Time;
end

% Convert Surface and Bottom variables from propagation time into image pixels
Surface_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Surface);
Surface_Mult_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), 2 * mdata{1}.Surface);
Bottom_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.Bottom);
twtt_bin = interp1(mdata{1}.Time, 1:length(mdata{1}.Time), mdata{1}.twtt);
ice_mask = mdata{1}.ice_mask;
Bottom_bin(isnan(Bottom_bin)) = -1;
twtt_bin(isnan(twtt_bin)) = -1;

%% Automated labeling section
% =========================================================================
% Specify which range lines to browse
skip = 1;
rlines = 1:skip:size(mdata{1}.Topography.img,3);

% Template size
ms = 11;

% Training parameters and preparing fusion slices
mu = [];
sigma = [];
fusion_slices = [];
tic;
for rline = rlines
    fprintf('rline: %d (%f sec)\n', rline, toc);
    fusion = tomo.fuse(double(db(mdata{1}.Topography.img(:,:,rline))), ...
        double(db(mdata{2}.Topography.img(:,:,rline))), ...
        double(db(mdata{3}.Topography.img(:,:,rline))));
    fusion(fusion>27) = 27;
    fusion = reshape(fusion, size(db(mdata{1}.Topography.img(:,:,rline))));
    [m, s] = tomo.train_params(fusion, double(twtt_bin(:,rline)));
    fusion_slices(:,:,rline) = fusion;
    mu = [mu; m];
    sigma = [sigma; s];
end

mdata_combined = mdata{nadir_img};
mdata_combined.Topography.img = fusion_slices;
mdata_combined.Topography.mu = mu;
mdata_combined.Topography.sigma = sigma;
mdata_combined.ice_mask = ice_mask;
fn_combined = fullfile(fn_dir,sprintf('Data_20140401_03_044.mat'));
save(fn_combined,'-struct','mdata_combined');

%% Automated cycling section
% =========================================================================
% 3D surface
bottom_surface = [];

for rline = rlines
    fprintf('rline: %d (%f sec)\n', rline, toc);
    labels = tomo.detect(fusion_slices(:,:,rline), double(twtt_bin(:,rline)), ...
        double(Bottom_bin(rline)), [], double(ice_mask(:,rline)), double(mean(mu)), double(mean(sigma)));
    
    % If empty, assign random values
    if isempty(labels)
        bottom_surface = [bottom_surface; randi([1 700], 1, size(bottom_surface, 2))];
    else
        bottom_surface = [bottom_surface; labels];
    end
    
%     set(h_image(1),'CData',fusion_slices(:,:,rline));
%     set(h_title(1),'String',sprintf('%d',rline));
%     set(h_surf_plot(1),'YData',Surface_bin(rline))
%     set(h_bot_plot(1),'YData',Bottom_bin(rline))
%     set(h_surf_mult_plot(1),'YData',Surface_Mult_bin(rline))
%     set(h_dem_plot(1),'YData',twtt_bin(:,rline))
%     set(h_ice_bed_plot(1),'YData',labels)
%     ylim(Surface_bin(rline)+[-25 500]);
%     set(hplot,'XData',[rline rline]);
    
    % keyboard
end

save(fn_combined,'-append','bottom_surface');


