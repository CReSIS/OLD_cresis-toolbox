% script run_create_movie
%
% Run script for tomo.create_movie
%
% Authors: Victor Berger, John Paden
%
% See also: tomo.create_movie.m

%% User parameters
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');
% params.cmd.generic = 1;
% params.cmd.frms = [];

options.geotiff_fn = ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90.tif');
options.sys        = 'rds';
options.location   = 'arctic';

% Set to true if crossovers should be shown
options.showCrossovers = true;

% Set to true if map should be shown
options.showMAP = true;
options.pathMAP = ct_filename_gis(params(1),fullfile('canada','Landsat-7','Canada_90m.tif'));

% DEM_source standard value = 'DEM'
options.DEM_source = 'DEM';

% Set to true if a video should be created
options.createVideo = true;

% Set location to save video 
% Has no effect if video is not being created
options.videoPath = 'movies_TESTING';

% Set video format
% Recommended: 'MPEG-4'
options.videoFormat = 'MPEG-4';

% Set video size
% Recommended (720p): [1280 720]
options.videoSize = [1280 720];

% frameRate standard value = 30
% Changes video "speed"
% Has no effect if video is not being created
options.frameRate = 30;

% frameSkip standard value = 1
% Changes video "speed"
options.frameSkip = 3;

% zoom standard value = .25
% Changes video magnification
options.zoom = .25;

%% Automated section
tomo.create_movie(params, options);


