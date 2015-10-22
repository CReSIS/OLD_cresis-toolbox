function make_layers_from_grid(param,proj,X,Y,layers)
% make_layers_from_grid(param,proj,X,Y,layers)
%
% Makes layer for OPS from grid map. This function calls
% ops_create_layer_points.m
% Input:
%   param: a structure the controls how to create layers 
%     .param_fn = file name of the parameter spreadsheets for a field season,
%                 set the generic to 1 for the selected segments and frames.
%     .post_dir = set to CSARP_post to operate in the post directory, otherwise 
%                 leave this empty
%     .combine_chan_wf_input = which radar echograms will be used to create the
%                              layers (usually csarp-combined, standard, qlook, or mvdr)
%     .plot = set to true for displaying echograms with layers, set to
%             false for creating layers faster (at the same time, set param.ops_en to true 
%     .ops_en = set to true to insert layer points to database; set to
%               false (at the same time to set param.plot to true) for
%               debugging or checking the layers created and the layer
%               points will not insert into database
%   proj: the projection structure of the grid map
%   X:  the vector of the x coordinates of the grid
%   Y:  the vector of the x coordinates of the grid
%   layers: a structure with layer information saved into database
%     .username = a cell array of string tells who creates or uses the layers
%     .name = a cell array of string gives the names for the layers to be created and inserted
%     .value = a three dimesntional array of layer elevation in meters; the
%      first dimention is the size of X, the second the size of Y and the
%      third the size of layers.name
%     .layers.surface_flag = an integer array tells a layer is the surface
%      ( = 1) or in ice ( = 0), the information needded to convert the
%      elevation to 2-way propagation time
%
% Output: layer points in OPS
%
% Examples:
%  Test code at bottom of this file
%
% Author: Jilu Li,John Paden
%
% See also: ops_create_layer_points

physical_constants;
if ~isfield(param,'post_dir')
  param.post_dir = '';
end
params = read_param_xls(param.param_fn);
param.radar_name = params(1).radar_name;
param.season_name = params(1).season_name;
[fn_pathstr, fn_name, fn_ext] = fileparts(param.param_fn);

% assign constant layer properties
if strfind(fn_name,'rds')
  sys = 'rds';
elseif strfind(fn_name,'accum')
  sys = 'accum';
elseif strfind(fn_name,'snow')
  sys = 'snow';
elseif strfind(fn_name,'kuband')
  sys = 'kuband';
else
  error('radar type not supported');
end
if strfind(param.season_name,'Greenland')
  param_lyr.properties.location = 'arctic';
elseif strfind(param.season_name,'Antarctica')
  param_lyr.properties.location = 'antarctic';
else
  error('location not supported');
end

param.day_seg = [];
for param_idx = 1:length(params)
  day_seg = params(param_idx).day_seg;
  if ~params(param_idx).cmd.generic
    continue;
  end
  
  param.day_seg = day_seg;
  input_dir = ct_filename_out(param,'', ...
    fullfile(param.post_dir,['CSARP_' param.combine_chan_wf_input]),'');
  
  if ~exist(input_dir,'dir')
    continue;
  end
  fns = get_filenames(input_dir,sprintf('Data_%s',day_seg(1:4)),'','.mat');

  img_01_only = false;
  if isempty(fns)
    fns = get_filenames(input_dir,sprintf('Data_img_01_%s',day_seg(1:4)),'','.mat');
    img_01_only = true;
  end
  
  frm_list = [];
  
  % Using the standard naming convention, extract the frame number
  % from each file
  for fn_idx = 1:length(fns)
    frm_list(fn_idx) = str2double(fns{fn_idx}(end-6:end-4));
  end
  [frm_list file_list] = unique(frm_list);
  
  % Create each layer for selected frames
  for frm_idx = 1:length(frm_list)
    frm = frm_list(frm_idx);
    if ~isempty(params(param_idx).cmd.frms) && isempty(find(frm==params(param_idx).cmd.frms))
      continue;
    end
    frm_id = [day_seg,sprintf('_%03d',frm)];
    
    % Load frame data file
    if param.plot
    load(fns{file_list(frm_idx)},'GPS_time','Latitude','Longitude','Elevation','Surface','Time','Data');
    else
    load(fns{file_list(frm_idx)},'GPS_time','Latitude','Longitude','Elevation','Surface');
    end
    % create layer from the grid by interpolation
    [x,y] = projfwd(proj,Latitude,Longitude);
    x_idxs = find(X>min(x)-proj.PixelScale(1) & X<max(x)+proj.PixelScale(1));
    y_idxs = find(Y>min(y)-proj.PixelScale(2) & Y<max(y)+proj.PixelScale(2));
    if isempty(x_idxs) | isempty(y_idxs)
      fprintf( 'Out of the grid, no layers for frame %s\n', frm_id);
      continue
    end
    X_loc = X(x_idxs);
    Y_loc = Y(y_idxs);
    [X_loc_mesh,Y_loc_mesh] = meshgrid(X_loc,Y_loc);
    sf_i = [];
    for lyr_idx = 1:length(layers.name)
      layer_loc = layers.value(x_idxs,y_idxs,lyr_idx);
      layer_i = interp2(X_loc_mesh,Y_loc_mesh,layer_loc',x,y);   
      if layers.surface_flag(lyr_idx)
        sf_elev = layer_i;
        layer_i =  2*(Elevation-layer_i)/c;
        sf_i = layer_i;
      else
        if ~isempty(sf_i)
          layer_i = sf_i + 2*(sf_elev - layer_i)/(c/sqrt(er_ice));
        else
          layer_i = Surface + 2*(Elevation - Surface*c/2 - layer_i)/(c/sqrt(er_ice));
        end
      end
      
      good_idxs = find(~isnan(layer_i) == 1);
      if isempty(good_idxs)
        fprintf('No layer %s for frame %s\n', layers.name{lyr_idx},frm_id);
        continue
      else
        GPS_time = GPS_time(good_idxs);
        Longitude = Longitude(good_idxs);
        Latitude = Latitude(good_idxs);
        Elevation = Elevation(good_idxs);
        layer_i = double(layer_i(good_idxs));
      end
      % ---------------------------------------------------------------------
      % Assign layer properties
      param_lyr.geometry.coordinates = [Longitude',Latitude',Elevation'];
      param_lyr.properties.username = layers.username{lyr_idx};
      param_lyr.properties.segment = day_seg;
      param_lyr.properties.gps_time = GPS_time; 
      param_lyr.properties.twtt = layer_i;
      param_lyr.properties.type = 2*ones(size(layer_i));          % auto
      param_lyr.properties.quality = 3*ones(size(layer_i));       % derived
      param_lyr.properties.lyr_name = layers.name{lyr_idx};
      % Create layer points
      fprintf('Creating layer ''%s'' of frame %s\n', layers.name{lyr_idx},frm_id);
      if param.plot
        figure(2);imagesc(GPS_time,Time,lp(Data));
        title([frm_id(1:8),'\',frm_id(9:11),'\',frm_id(12:15)]);
        figure(2);hold on;plot(GPS_time,layer_i,'o');hold off
        pause(0.5)
      end
      if param.ops_en
        [status,message] = ops_create_layer_points(sys,param_lyr);
      end
    end
  end
end
return;

%% Example 1: mass conservation grid
% Use parameters spreadsheet to select segment and frame list for creating layers
% Set the generic to 1 for the selected segments and frames
param.param_fn = '/users/jliwestc/scripts/matlab/personal/rds_param_2013_Greenland_P3.xls';

% Load the grid
grid_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/mass_conservation/MCdatset-2013-12-03.nc';
gridmap = ncinfo(grid_fn)
proj = geotiffinfo('/cresis/snfs1/dataproducts/GIS_data/greenland/Landsat-7/Greenland_natural.tif');
proj.PixelScale = [gridmap.Attributes(14).Value;gridmap.Attributes(14).Value;0];
X = ncread(grid_fn,'x');
Y = ncread(grid_fn,'y');
X = single(X);
Y = single(Y);

% layer information from the grid
layers.username = {'jliwestc','jliwestc'};    %  the name of the person who creates or uses the layers
layers.name = {'sf_mass_cv','bt_mass_cv'};
layers.value = single(zeros(length(X),length(Y),length(layers.name)));
layers.value(:,:,1) = single(ncread(grid_fn,'surface')); 
tmp = single(ncread(grid_fn,'bed'));
tmp(tmp == -9999) = NaN;
layers.value(:,:,2) = tmp;
layers.surface_flag = [1,0];

% param.post_dir: occasionally you may want to operate in the post directory
%   set this to CSARP_post, otherwise leave empty
param.post_dir = '';


% .combine_chan_wf_input = which radar echograms will be used to create the
% layers (usually csarp-combined, standard, qlook, or mvdr)
param.combine_chan_wf_input = 'standard';

% plot inserted layers for debug or monitoring
param.plot = true; 
param.ops_en = false;

make_layers_from_grid(param,proj,X,Y,layers);


%% Example 2: Joe Plummer's Jakobshavn grid
% Use parameters spreadsheet to select segment and frame list for creating layers
% Set the generic to 1 for the selected segments and frames
param.param_fn = '/users/jliwestc/scripts/matlab/personal/rds_param_2006_Greenland_TO.xls';

% Load the grid
grid_fn = '/cresis/snfs1/dataproducts/GIS_data/greenland/plummer_jakobshavn/jak_grid.tif';
proj = geotiffinfo(grid_fn);
[bt_elev,R] = geotiffread(grid_fn);
bt_elev(bt_elev == min(min(bt_elev))) = NaN;
X = R.XLimWorld(1) + [R.XLimIntrinsic(1):R.XLimIntrinsic(2)-1]'*R.DeltaX;
Y = R.YLimWorld(2) + [R.YLimIntrinsic(1):R.YLimIntrinsic(2)-1]'*R.DeltaY;

% layer information from the grid
layers.username = {'jliwestc'};     %the name of the person who creates the layers
layers.name = {'bt_plummer'};
layers.value = zeros(length(X),length(Y),length(layers.name));
layers.value(:,:,1) = bt_elev;
layers.surface_flag = 0;


% param.post_dir: occasionally you may want to operate in the post directory
%   set this to CSARP_post, otherwise leave empty
param.post_dir = '';


% .combine_chan_wf_input = which radar echograms will be used to create the
% layers (usually csarp-combined, standard, qlook, or mvdr)
param.combine_chan_wf_input = 'standard';

% plot inserted layers for debug or monitoring
param.plot = true; 
param.ops_en = false;

make_layers_from_grid(param,proj,X,Y,layers);



