% script run_load_data_by_gps_time
%
% This script prepares data for publication. It is an example script for
% setting up parameters and calling load_data_by_gps_time which does the
% real work.
%
% To use this function, go through and set the parameters highlighted by:
%   "<== CHANGE HERE"
%
% Author: John Paden

%% User Settings
% ====================================================================
clear param echo_param;

% data_load_method: string containing "file" or "arbitrary"
%   file: Loads a data frame and plots the whole data frame
%   arbitrary: Allows a specific GPS date range to be specified
data_load_method = 'file'; % <== CHANGE HERE

if strcmpi(data_load_method,'arbitrary')
  %% Load data for an arbitrary time period
  
  param.start.hour = 13; % <== CHANGE HERE
  param.start.minute = 17; % <== CHANGE HERE
  param.start.sec = 27; % <== CHANGE HERE
  
  param.stop.hour = 13; % <== CHANGE HERE
  param.stop.minute = 26; % <== CHANGE HERE
  param.stop.sec = 46; % <== CHANGE HERE
  
  param.start.year = 2009; % <== CHANGE HERE
  param.start.month = 4; % <== CHANGE HERE
  param.start.day = 9; % <== CHANGE HERE
  param.stop.year = param.start.year;
  param.stop.month = param.start.month;
  param.stop.day = param.start.day;
  
  param.start.gps_time = datenum_to_epoch(datenum(param.start.year,param.start.month,param.start.day,param.start.hour,param.start.minute,param.start.sec));
  param.stop.gps_time = datenum_to_epoch(datenum(param.stop.year,param.stop.month,param.stop.day,param.stop.hour,param.stop.minute,param.stop.sec));
  
elseif strcmpi(data_load_method,'file')
  %% Load data associated with a frame
  load('X:\ct_data\rds\2009_Greenland_TO\CSARP_post\CSARP_mvdr\20090409_01\Data_20090409_01_008.mat','GPS_time') % <== CHANGE HERE
  param.start.gps_time = GPS_time(1);
  param.stop.gps_time = GPS_time(end);
  
else
  error('Unsupported data_load_method.\n');
end

% param.radar_name: 'accum', 'kaband', 'kuband', 'rds', or 'snow'
param.radar_name = 'rds'; % <== CHANGE HERE

param.season_name = '2009_Greenland_TO'; % <== CHANGE HERE

% echo_param.elev_comp: Elevation compensation 0=none, 1=relative, 2=surface flattened, 3=WGS-84
echo_param.elev_comp = 3; % <== CHANGE HERE

% param.out: output data product to use. For example:
%   'qlook', 'standard', 'mvdr', 'CSARP_post/standard', 'CSARP_post/mvdr'
param.out = 'CSARP_post/mvdr'; % <== CHANGE HERE

% param.img_name: output data product image. For example:
%   '': combined product, 'img_01_', , 'img_02_'
param.img_name = '';
  
if echo_param.elev_comp == 3
  echo_param.depth = '[min(Surface_Elev)-1500 max(Surface_Elev)+20]'; % <== CHANGE HERE
  %echo_param.depth = '[1594.9 1635.1]';
else
  echo_param.depth = '[min(Surface_Depth)-20 max(Surface_Depth)+200]'; % <== CHANGE HERE
  %echo_param.depth = '[-5 120]';
end
echo_param.depth_offset = 0;
echo_param.time_offset = 0;
echo_param.caxis = [];
echo_param.er_ice = 3.15;

% echo_param.axis_type: Units representation 'bars' or 'standard'
echo_param.axis_type = 'standard';

% param.create_map_en: set to true to produce a map
param.create_map_en = false;

% param.post.map.type: 'combined', 'singles', 'contour', 'contour-singles', or 'ASCAT'
param.post.map.type = 'combined';

% param.post.map.location: map choice 'Greenland', 'Antarctica', 'Arctic', 'Norway', or 'Canada'
param.post.map.location = 'Greenland';

% Enable for detrending
if 0
  echo_param.detrend.mode = 'polynomial';
  echo_param.detrend.poly_order = 4;
  echo_param.detrend.depth = '[min(Surface_Depth)+1 max(Surface_Depth)+20.1]';
  %echo_param.filter = inline('wiener2(N,[3 31])');
end
% echo_param.filter = inline('fir_dec(N,ones([1 11])/11,1)');

% param.use_master_surf: set to true if you want to use the previous
% surface loaded by this function (this is useful when comparing to
% different systems and you want the second system loaded to use the
% surface from the first system loaded)
param.use_master_surf = 0;

%% Automated Section
% ====================================================================
% ====================================================================

param.save_files = false;
out_fn_dir = fullfile('~/',datestr(epoch_to_datenum(param.start.gps_time),'yyyymmdd'));

%% Determine the output figure handle
% ====================================================================
[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

if strcmp(output_dir,'accum')
  echo_param.fig_hand = 2;
elseif strcmp(param.radar_name,'kaband')
  echo_param.fig_hand = 2;
elseif strcmp(param.radar_name,'kuband')
  echo_param.fig_hand = 4;
elseif strcmp(param.radar_name,'rds')
  echo_param.fig_hand = 5;
elseif strcmp(param.radar_name,'snow')
  echo_param.fig_hand = 6;
else
  error('Unsupported radar type %s', param.radar_name);
end

%% Prepare variables and call load_data_by_gps_time.m
% ====================================================================

physical_constants;

ds = load_data_by_gps_time(param);

Data = ds.Data;
Latitude= ds.Latitude;
Longitude = ds.Longitude;
Elevation = ds.Elevation;
GPS_time = ds.GPS_time;
Surface = ds.Surface;
Time = ds.Time;

% Create output directory
if param.save_files && ~exist(out_fn_dir,'dir')
  fprintf('Creating output directory %s\n', out_fn_dir);
  mkdir(out_fn_dir);
end

% Create time-stamp strings used in filenames
start_time_stamp_str = datestr(epoch_to_datenum(param.start.gps_time),'yyyymmdd_HHMMSS');
stop_time_stamp_str = datestr(epoch_to_datenum(param.stop.gps_time),'yyyymmdd_HHMMSS');

%% Create map template using publish_map.m
% =========================================================================
if ~strcmpi(param.post.map.type,'none') && param.create_map_en
  fprintf('Creating map template (%s)\n', datestr(now));
  
  % --------------------------------------------------
  % Create base map to publish in with publish_map
  map_param = [];
  map_param.type = param.post.map.type;
  map_param.location = param.post.map.location;
  map_param.fig_hand = 1;
  map_param.map_title = '';
  map_param.decimate_seg = false;
  map_param.resample = false;
  % structure fields necessary for ASCAT maps
  if strcmpi(param.post.map.type, 'ASCAT')
    map_param.year = param.start.year;
    map_param.month = param.start.month;
    map_param.day = param.start.day;
  end
  map_info = publish_map('setup',map_param);
  
  % --------------------------------------------------
  % Create flight path variables for publish_map
  vectors.lat = Latitude;
  vectors.lon = Longitude;
  
  [x,y] = projfwd(map_info.proj,vectors.lat,vectors.lon);
  day_segs_idx(1) = 1;
  day_seg_x{1} = x/1e3;
  day_seg_y{1} = y/1e3;
  
  % --------------------------------------------------
  % Plot flightlines over geotiff using publish_map
  map_param.day_seg_x = day_seg_x{1}*1000;
  map_param.day_seg_y = day_seg_y{1}*1000;
  [frame_X,frame_Y] = projfwd(map_info.proj,Latitude,Longitude);
  map_param.frame_X = {frame_X};
  map_param.frame_Y = {frame_Y};
  map_info = publish_map('delete',map_param,map_info);
  map_info = publish_map('plot',map_param,map_info);
  
  % --------------------------------------------------
  % Save map file
  time_stamp_str = datestr(epoch_to_datenum(param.start.gps_time),'yyyymmdd_HHMMSS');
  out_fn = fullfile(out_fn_dir, sprintf('%s_%s_0map.jpg',start_time_stamp_str, ...
    stop_time_stamp_str));
  if param.save_files
    fprintf('Saving to %s\n', out_fn);
    saveas(map_param.fig_hand,out_fn);
  end
end

%% Uniform sampling in along-track of echogram
% ===================================================================
along_track = geodetic_to_along_track(Latitude,Longitude,0*Elevation);
if along_track(end)-along_track(1) > 100
  along_track = geodetic_to_along_track(Latitude,Longitude,0*Elevation,50);
end
Nx = numel(along_track);
dx = (along_track(end)-along_track(1)) / (Nx-1);
new_along_track = dx*(0:Nx-1);
Data = interp1(along_track,Data.',new_along_track,'linear','extrap').';
Latitude = interp1(along_track,Latitude,new_along_track,'linear','extrap').';
Longitude = interp1(along_track,Longitude,new_along_track,'linear','extrap').';
Elevation = interp1(along_track,Elevation,new_along_track,'linear','extrap').';
GPS_time = interp1(along_track,GPS_time,new_along_track,'linear','extrap').';
Surface = interp1(along_track,Surface,new_along_track,'linear','extrap').';

%% Create echogram plot
% ===================================================================

mdata.Data = Data;
mdata.Latitude= Latitude;
mdata.Longitude = Longitude;
mdata.Elevation = Elevation;
mdata.GPS_time = GPS_time;
mdata.Surface = Surface;
mdata.Time = Time;

echo_param.num_x_tics = 4;
echo_param.frm_id = '';

lay.GPS_time = GPS_time;
lay.Elevation = Elevation;
lay.layerData{1}.value{1}.data = NaN*zeros(size(Surface));
if param.use_master_surf
  % Use master surface (useful when comparing multiple radar outputs
  % which each have their own surface variable and a common surface should be used
  if exist('master_surf_filt','var')
    lay.layerData{1}.value{2}.data = interp1(master_surf_gps_time,master_surf_filt,lay.GPS_time,'linear','extrap');
    if master_surf_gps_time(1)-5 > lay.GPS_time(1) ...
        || master_surf_gps_time(end)+5 < lay.GPS_time(end)
      warning('Using master surface on nonmatching GPS times');
      keyboard
    end
  else
    master_surf_filt = Surface;
    master_surf_gps_time = lay.GPS_time;
    lay.layerData{1}.value{2}.data = Surface;
  end
else
  lay.layerData{1}.value{2}.data = Surface;
end
lay.layerData{2}.value{1}.data = NaN*zeros(size(Surface));
lay.layerData{2}.value{2}.data = NaN*zeros(size(Surface));
echo_info = publish_echogram(echo_param,mdata,lay);
season_name = param.season_name;
season_name(season_name=='_') = ' ';
title(sprintf('%s %s %s to %s', param.radar_name, season_name, ...
  datestr(epoch_to_datenum(param.start.gps_time)), datestr(epoch_to_datenum(param.stop.gps_time),'HH:MM:SS')));
% set(echo_info.echo_title,'Visible', 'off');
set(echo_info.h_surf,'Visible','off');
set(echo_info.h_bot,'Visible','off');

%% Save echogram file
if param.save_files
  out_fn = fullfile(out_fn_dir, sprintf('%s_%s_%s.jpg',start_time_stamp_str, ...
    stop_time_stamp_str, param.radar_name));
  fprintf('Saving to %s\n', out_fn);
  saveas(echo_param.fig_hand,out_fn);
end

return;
