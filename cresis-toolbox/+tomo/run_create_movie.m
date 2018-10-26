% script run_create_movie
%
% Run script for tomo.create_movie
%
% Authors: Victor Berger, John Paden
%
% See also: tomo.create_movie.m

%% User parameters
params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'20140325_07','post');
params.cmd.generic = 1;
params.cmd.frms = 1 : 4;

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

% Set main frame oversampling rate
options.main_oversampling_rate = 3;

% Set crossover oversampling rate
options.CO_oversampling_rate = 2;

% Set to true if a video should be created
options.createVideo = true;

% Set location to save video 
% Has no effect if video is not being created
options.videoPath = 'movies_HQ';

% Set video format
% Recommended: 'MPEG-4'
% options.videoFormat = 'MPEG-4';
% options.videoFormat = 'Motion JPEG 2000';

% Set video size
% Recommended (720p): [1280 720]
% Recommended (1080p): [1920 1080] for high-res
options.videoSize = [1920 1080]; 

% frameRate standard value = 30
% Changes video "speed"
% Has no effect if video is not being created
options.frameRate = 30; 

% frameSkip standard value = 1
% Changes video "speed"
options.frameSkip = 10;

% zoom standard value = .25
% Changes video magnification
options.zoom = .25;

% Set to true if want to save a MAT file
%  corresponding to the video generated.
%  See memory dumping options
%  Default: true
options.saveMAT = true; % needs createVideo == 1

% Options for dumping memory 
% May run out of RAM if not done.
% Video MAT is split into separate parts
% Defaults: true, 1000
options.memory_dump = true;
options.memory_dump_size = 1000; % Arbitrary. 1000 works well

% Colormap limits
% Leave empty to use the entire dynamic range of each frame
options.cmaplim = [-250 2250]; 

%% Automated section
tomo.create_movie(params, options);


