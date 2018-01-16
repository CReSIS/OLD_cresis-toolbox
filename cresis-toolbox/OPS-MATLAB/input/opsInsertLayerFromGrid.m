function opsInsertLayerFromGrid(grid_param,layers)
% opsInsertLayerFromGrid(grid_param,layers)
%
% Makes layer for OPS from grid map using interp2. This function calls
%
% Input:
% grid_param: a structure the controls how to create layers
%  .param_fn = file name of the parameter spreadsheets for a field season,
%    set the generic column of the cmd worksheet to 1 for the
%    selected segments and frames.
%  .save_sources = string array containing layer source to update. Can be 'echogram',
%    'layerData', or 'ops'.
%  .echogram_source = string to source for echogram source option.
%    For example 'CSARP_post/qlook', 'standard', etc
%  .layerdata_source = string to source for layerData source option.
%    For example 'CSARP_post/layerData', 'layerData', etc
%  .plot_en = set to true for displaying echograms with layers, set to
%    false for creating layers faster
%   proj: the projection structure of the grid map
%   save_mode: string containing 'overwrite' or 'merge' (merge makes all
%     NaN in the new data source be unchanged)
% layers: structure array with layer information
%  .username = string containing username for this layer;
%  .name = string with layer name;
%  .group = string with layer group;
%  .ref_name = string with reference layer name (optional),
%    often set to 'surface'
%  .ref_eval = string with command to evaluate that combines new layer
%     with reference.  Variables available to the command are:
%       data: this is the data passed in after it is interpolated to flightline
%       ref_data: this is the reference layer data (twtt, sec)
%       elev: platform elevation (m)
%     Examples:
%       ref_eval = 'data = ref_data + data;'
%       ref_eval = 'data = data - (elev - ref_data*149896229)/84456957 + ref_data;'
%     not used.  If not empty, eval(layers(lyr_idx).ref_data) will be run
%     where ref_data contains ref layer data and data contains grid data
%     interpolated to the line.
%  .echogram_layer_name = string containing the variable name to store
%    into the echogram file (usually 'Surface' or 'Bottom')
%  .layerdata_layer_name = string that indicates which layer to save the
%    data to. For example, 'surface' (idx 1), 'bottom' (idx 2).
%  .x = x-vector for data matrix
%  .y = y-vector for data matrix
%  .data = twtt data matrix (or data matrix which will be used by ref_data)
%
% Output: layer points in OPS
%
% Examples:
%  Examples at bottom of this file
%
% Author: Jilu Li, John Paden
%
% See also: opsCreateLayerPoints.m, opsInsertLayerFromPointCloud.m

%% Load param worksheet
params = read_param_xls(grid_param.param_fn,[],'post');

param.day_seg = [];
for param_idx = 1:length(params)
  param = params(param_idx);
  
  %% Determine if this segment should be processed or not
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  % Load frames file
  load(ct_filename_support(param,'','frames'));
  
  %% Determine which frames need to be processed
  if isempty(param.cmd.frms)
    param.cmd.frms = 1:length(frames.frame_idxs);
  end
  % Remove frames that do not exist from param.cmd.frms list
  [valid_frms,keep_idxs] = intersect(param.cmd.frms, 1:length(frames.frame_idxs));
  if length(valid_frms) ~= length(param.cmd.frms)
    bad_mask = ones(size(param.cmd.frms));
    bad_mask(keep_idxs) = 0;
    warning('Nonexistent frames specified in param.cmd.frms (e.g. frame "%g" is invalid), removing these', ...
      param.cmd.frms(find(bad_mask,1)));
    param.cmd.frms = valid_frms;
  end
  
  if any(strcmpi(grid_param.save_sources,'ops'))
    %% Get all the frames for this segment
    sys = ct_output_dir(param.radar_name);
    ops_param = struct('properties',[]);
    ops_param.properties.season = param.season_name;
    ops_param.properties.segment = param.day_seg;
    [status,ops_seg_data] = opsGetSegmentInfo(sys,ops_param);
  end
  
  %% Update each of the frames
  for frm_idx = 1:length(param.cmd.frms)
    frm = param.cmd.frms(frm_idx);
    
    fprintf('  Updating %s frame %03d (%d of %d) (%s)\n', param.day_seg, ...
      frm, frm_idx, length(param.cmd.frms), datestr(now,'HH:MM:SS'));
    data_fn = fullfile(ct_filename_out(param,grid_param.echogram_source,''), ...
      sprintf('Data_%s_%03d.mat', param.day_seg, frm));
    if grid_param.plot_en
      plot_vars = load(data_fn,'Time','Data');
    end
    
    if any(strcmpi(grid_param.save_sources,'echogram'))
      %% 1. Open the specific echogram data file
      mdata = load(data_fn,'GPS_time','Latitude','Longitude','Elevation');
      
      %% 2. Interpolate the grid onto this echogram
      [x,y] = projfwd(grid_param.proj,mdata.Latitude,mdata.Longitude);
      
      %% 3. Save the layer back to the file under grid_param.echogram_layer_name
      for layer_idx = 1:length(layers)
        data = interp2(reshape(layers(layer_idx).x,[1 numel(layers(layer_idx).x)]), ...
          reshape(layers(layer_idx).y,[numel(layers(layer_idx).y) 1]),layers(layer_idx).data,x,y);
        if isfield(layers(layer_idx),'ref_name')
          tmp = load(data_fn,layers(layer_idx).ref_name);
          ref_data = tmp.(layers(layer_idx).ref_name);
          ref_data = interp_finite(ref_data);
          elev = mdata.Elevation;
          eval(layers(layer_idx).ref_eval);
        end
        
        if strcmpi(grid_param.save_mode,'merge')
          data(isnan(data)) = mdata.(layers(layer_idx).echogram_layer_name)(isnan(data));
        end
        mdata.(layers(layer_idx).echogram_layer_name) = data;
        
        if grid_param.plot_en
          figure(1); clf;
          imagesc([],plot_vars.Time,lp(plot_vars.Data));
          colormap(1-gray(256));
          hold on;
          plot(mdata.(layers(layer_idx).echogram_layer_name));
          hold off;
          keyboard
        end
        
        save(data_fn,'-append','-struct','mdata',layers(layer_idx).echogram_layer_name);
      end
    end
    
    if any(strcmpi(grid_param.save_sources,'layerdata'))
      %% 1. Open the specific layer data file
      layer_fn = fullfile(ct_filename_out(param,grid_param.layerdata_source,''), ...
        sprintf('Data_%s_%03d.mat', param.day_seg, frm));
      if ~exist(layer_fn,'file')
        error('Layer file %s does not exist', layer_fn);
      end
      % Load the layerData file
      lay = load(layer_fn);
      
      %% 2. Interpolate the grid onto this layerData
      [x,y] = projfwd(grid_param.proj,lay.Latitude,lay.Longitude);
      
      for layer_idx = 1:length(layers)
        data = interp2(reshape(layers(layer_idx).x,[1 numel(layers(layer_idx).x)]), ...
          reshape(layers(layer_idx).y,[numel(layers(layer_idx).y) 1]),layers(layer_idx).data,x,y);
        if isfield(layers(layer_idx),'ref_name')
          if strcmpi(layers(layer_idx).ref_name,'surface')
            ref_idx = 1;
          elseif strcmpi(layers(layer_idx).ref_name,'bottom')
            ref_idx = 2;
          else
            found = false;
            for ref_idx = 1:length(lay.layerData)
              if isfield(lay.layerData{ref_idx},'name') ...
                  && strcmpi(lay.layerData{ref_idx}.name,layers(layer_idx).ref_name)
                found = true;
                break;
              end
            end
            if ~found
              error('Reference layer %s not found\n', layers(layer_idx).ref_name);
            end
          end
          ref_data = lay.layerData{ref_idx}.value{2}.data;
          ref_data = interp_finite(ref_data);
          elev = lay.Elevation;
          eval(layers(layer_idx).ref_eval);
        end
        
        % Determine which index into lay.layerData needs to be updated
        if strcmpi(layers(layer_idx).name,'surface')
          lay_idx = 1;
          found = true;
        elseif strcmpi(layers(layer_idx).name,'bottom')
          lay_idx = 2;
          found = true;
        else
          found = false;
          for lay_idx = 1:length(lay.layerData)
            if isfield(lay.layerData{lay_idx},'name') ...
                && strcmpi(lay.layerData{lay_idx}.name,layers(layer_idx).name)
              found = true;
              break;
            end
          end
          if ~found
            % Add a new layer if layer does not exist
            lay_idx = length(lay.layerData) + 1;
            % Create manual points
            lay.layerData{lay_idx}.value{1}.data = NaN*zeros(size(data));
            % Create auto points
            lay.layerData{lay_idx}.value{2}.data = NaN*zeros(size(data));
            % Set quality to good
            lay.layerData{lay_idx}.quality.data = 1*ones(size(data));
          end
        end
        
        % Update automated points
        if strcmpi(grid_param.save_mode,'merge')
          data(isnan(data)) = lay.layerData{lay_idx}.value{2}.data(isnan(data));
        end
        lay.layerData{lay_idx}.value{2}.data = data;
        
        if grid_param.plot_en
          figure(1); clf;
          imagesc([],plot_vars.Time,lp(plot_vars.Data));
          colormap(1-gray(256));
          hold on;
          plot(lay.layerData{lay_idx}.value{2}.data);
          hold off;
          keyboard
        end
        
        %% 3. Save the layer back to the file under grid_param.layerdata_layer_name
        %    If layer names do not exist, then assume 'surface' is 1, 'bottom'
        %    is 2, and if a layer does not exist then a new layer is created.
        save(layer_fn,'-append','-struct','lay','layerData');
      end
    end
    
    if any(strcmpi(grid_param.save_sources,'ops'))
      for layer_idx = 1:length(layers)
        start_gps = ops_seg_data.properties.start_gps_time(frm);
        stop_gps = ops_seg_data.properties.stop_gps_time(frm);
        
        %% 1. Check to see if layer exists
        [status,data] = opsGetLayers(sys);
        if ~any(strcmpi(data.properties.lyr_name,layers(layer_idx).name))
          % Create the layer if it does not exist
          ops_param = [];
          ops_param.properties.lyr_name = layers(layer_idx).name;
          ops_param.properties.lyr_group_name = layers(layer_idx).group;
          ops_param.properties.lyr_description = layers(layer_idx).description;
          ops_param.properties.public = true;
          
          [status,ops_data] = opsCreateLayer(sys,ops_param);
        end
        
        %% 2. OPS query to get all the point path ID's
        ops_param = struct('properties',[]);
        ops_param.properties.location = param.post.ops.location;
        ops_param.properties.season = param.season_name;
        ops_param.properties.start_gps_time = start_gps;
        ops_param.properties.stop_gps_time = stop_gps;
        ops_param.properties.nativeGeom = true;
        [status,ops_data] = opsGetPath(sys,ops_param);
        
        %% 3. Interpolate the data onto the point path ID's
        [x,y] = projfwd(grid_param.proj,ops_data.properties.Y,ops_data.properties.X);
        data = interp2(reshape(layers(layer_idx).x,[1 numel(layers(layer_idx).x)]), ...
          reshape(layers(layer_idx).y,[numel(layers(layer_idx).y) 1]),layers(layer_idx).data,x,y);
        quality = round(interp2(reshape(layers(layer_idx).x,[1 numel(layers(layer_idx).x)]), ...
          reshape(layers(layer_idx).y,[numel(layers(layer_idx).y) 1]),layers(layer_idx).quality,x,y));
        
        %% 4. Load reference data
        if isfield(layers(layer_idx),'ref_name')
          ops_param = struct('properties',[]);
          ops_param.properties.location = param.post.ops.location;
          ops_param.properties.season = param.season_name;
          ops_param.properties.segment = param.day_seg;
          ops_param.properties.start_gps_time = start_gps;
          ops_param.properties.stop_gps_time = stop_gps;
          ops_param.properties.lyr_name = layers(layer_idx).ref_name;
          [status,ref_data] = opsGetLayerPoints(sys,ops_param);
          % Interpolate ref data onto point path
          [ref_data.properties.gps_time,sort_idxs] = sort(ref_data.properties.gps_time);
          ref_data.properties.twtt = ref_data.properties.twtt(sort_idxs);
          ref_data = interp1(ref_data.properties.gps_time,ref_data.properties.twtt, ...
            ops_data.properties.gps_time);
          ref_data = interp_finite(ref_data);
          elev = ops_data.properties.elev;
          eval(layers(layer_idx).ref_eval);
        end

        if strcmpi(grid_param.save_mode,'merge')
          % Remove NaN layer points
          ops_data.id = ops_data.id(~isnan(data));
          quality = quality(~isnan(data));
          data = data(~isnan(data));
        end
        
        ops_param = struct('properties',[]);
        ops_param.properties.point_path_id = ops_data.properties.id;
        ops_param.properties.twtt = double(data);
        ops_param.properties.type = 2*ones(size(ops_param.properties.twtt));
        ops_param.properties.quality = double(quality);
        ops_param.properties.lyr_name = layers(layer_idx).name;
        
        if grid_param.plot_en
          figure(1); clf;
          imagesc([],plot_vars.Time,lp(plot_vars.Data));
          colormap(1-gray(256));
          hold on;
          plot(ops_param.properties.twtt);
          hold off;
          figure(2); clf;
          imagesc(layers.x,layers.y,layers.data); set(gca,'ydir','normal');
          hold on;
          plot(x,y,'k')
          keyboard
        end
        
        opsCreateLayerPoints(sys,ops_param);
      end
    end
    
  end
  
end

return;

%% Example 1: mass conservation grid
% Use parameters spreadsheet to select segment and frame list for creating layers
% Set the generic to 1 for the selected segments and frames

param = [];
param.param_fn = ct_filename_param('rds_param_2009_Greenland_TO.xls');

% grid_fn = '/cresis/snfs1/scratch/paden/mass_conservation/HelheimStream-2014-11-18.nc';
grid_fn = '/cresis/snfs1/scratch/paden/mass_conservation/KangerdlugssuaqStream-2014-11-18.nc';

% Create geotiff projection structure
% 1. Grab from a geotiff file with the same projection
%param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));
% 2. Create an mstruct with the same projection
%[x,y,mstruct] = geodetic_to_stereographic(70,-45);
% 3. Hand construction of geotiff projection OR
param.proj = [];
param.proj.CornerCoords = [];
param.proj.Ellipsoid = [];
param.proj.PM = [];
param.proj.PMLongToGreenwich = [];
param.proj.Zone = [];
geoid = almanac('earth',ncreadatt(grid_fn,'polar_stereographic','ellipsoid'),'m');
param.proj.SemiMajor = geoid;
param.proj.SemiMinor = minaxis(geoid);
param.proj.ProjParm = zeros(7,1);
param.proj.ProjParm(1) = ncreadatt(grid_fn,'polar_stereographic','standard_parallel');
param.proj.ProjParm(2) = ncreadatt(grid_fn,'polar_stereographic','straight_vertical_longitude_from_pole');
param.proj.ProjParm(3) = 0;
param.proj.ProjParm(4) = 0;
param.proj.ProjParm(5) = 1;
param.proj.ProjParm(6) = ncreadatt(grid_fn,'polar_stereographic','false_easting');
param.proj.ProjParm(7) = ncreadatt(grid_fn,'polar_stereographic','false_northing');

param.proj.CTProjection = 'CT_PolarStereographic'; % Choose from list in toolbox/map/mapproj/private/projcode
param.proj.GeoTIFFCodes.CTProjection = 15; % CT_PolarStereographic from toolbox/map/mapproj/private/projcode

% Default GeoTIFFCode values
param.proj.GeoTIFFCodes.Model = int16(1);
param.proj.GeoTIFFCodes.PCS = int16(32767);
param.proj.GeoTIFFCodes.GCS = int16(32767);
param.proj.GeoTIFFCodes.UOMAngle = int16(32767);
param.proj.GeoTIFFCodes.Datum = int16(32767);
param.proj.GeoTIFFCodes.PM = int16(32767);
param.proj.GeoTIFFCodes.ProjCode = int16(32767);
param.proj.GeoTIFFCodes.Projection= int16(32767);
param.proj.GeoTIFFCodes.MapSys = int16(32767);

% Default units to meters for now
param.proj.GeoTIFFCodes.UOMLength = int16(9001);

layers = [];
layers(1).username = 'paden';
layers(1).name = 'mc_bottom';
layers(1).group = 'standard';
layers(1).description = 'Mass conservation ice bottom';
layers(1).ref_name = 'surface';
layers(1).ref_eval = 'data = ref_data + data;';
layers(1).echogram_layer_name = 'mc_bottom';
layers(1).layerdata_layer_name = 'mc_bottom';
layers(1).x = ncread(grid_fn,'x');
layers(1).y = ncread(grid_fn,'y');
physical_constants;
layers(1).data = ncread(grid_fn,'thickness').' / (c/2/sqrt(er_ice));

param.plot_en = false;
param.echogram_source = 'mvdr';
param.layerdata_source = 'layerData';
% save_sources: cell array of strings indicating which layer sources
%   should be updated (options are ops, layerData, and echogram)
param.save_sources = {'ops'};
param.save_mode = 'overwrite';
opsInsertLayerFromGrid(param,layers);

%% Example 2: Gogineni's Jakobshavn grid
% Use parameters spreadsheet to select segment and frame list for creating layers
% Set the generic to 1 for the selected segments and frames

param = [];
param.param_fn = ct_filename_param('rds_param_2009_Greenland_TO.xls');

grid_fn = '/cresis/snfs1/scratch/paden/mass_conservation/Jakobshavn_2006_2009_Composite/Jakobshavn_2006_2009_Composite/grids/jakobshavn_2006_2009_composite_thickness.txt';

% Load grid
fid = fopen(grid_fn,'r');
fseek(fid,0,-1);
grid = [];
for idx=1:6
  fields = textscan(fid,'%s%f',1);
  grid.(fields{1}{1}) = fields{2};
end
thickness = textscan(fid,'%f');
fclose(fid);
thickness = reshape(thickness{1},[grid.ncols grid.nrows]).';
thickness(thickness == grid.NODATA_value) = NaN;
imagesc(thickness);
x_axis = grid.xllcorner + grid.cellsize*(0:grid.ncols-1);
y_axis = grid.yllcorner + grid.cellsize*(0:grid.nrows-1);
y_axis = y_axis(end:-1:1);

% Create geotiff projection structure
% 1. Grab from a geotiff file with the same projection
param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));

layers = [];
layers(1).username = 'paden';
layers(1).name = 'gogineni2014';
layers(1).group = 'standard';
layers(1).description = 'Gogineni JofG 2014 grid';
layers(1).ref_name = 'surface';
layers(1).ref_eval = 'data = ref_data + data;';
layers(1).echogram_layer_name = 'gogineni2014';
layers(1).layerdata_layer_name = 'gogineni2014';
layers(1).x = x_axis;
layers(1).y = y_axis;
physical_constants;
layers(1).data = thickness / (c/2/sqrt(er_ice));

param.plot_en = false;
param.echogram_source = 'mvdr';
param.layerdata_source = 'layerData';
% save_sources: cell array of strings indicating which layer sources
%   should be updated (options are ops, layerData, and echogram)
param.save_sources = {'ops'};
param.save_mode = 'overwrite';
opsInsertLayerFromGrid(param,layers);

%% Example 3: Joe Plummer's Jakobshavn grid
% Use parameters spreadsheet to select segment and frame list for creating layers
% Set the generic to 1 for the selected segments and frames

param = [];
param.param_fn = ct_filename_param('rds_param_2009_Greenland_TO.xls');

grid_fn = ct_filename_gis([],'greenland/plummer_jakobshavn/jak_grid.tif');

param.proj = geotiffinfo(grid_fn);

% Load the grid
proj = geotiffinfo(grid_fn);
[bt_elev,R] = geotiffread(grid_fn);
bt_elev(bt_elev == min(min(bt_elev))) = NaN;
X = R.XLimWorld(1) + [R.XLimIntrinsic(1):R.XLimIntrinsic(2)-1]'*R.DeltaX;
Y = R.YLimWorld(2) + [R.YLimIntrinsic(1):R.YLimIntrinsic(2)-1]'*R.DeltaY;

layers = [];
layers(1).username = 'paden';
layers(1).name = 'plummer';
layers(1).group = 'standard';
layers(1).description = 'Plummer Jakobshavn Grid';
layers(1).ref_name = 'surface';
layers(1).ref_eval = sprintf('data = (elev - ref_data*%.0f)/%.0f - data + ref_data;',c/2,c/2/sqrt(er_ice));
layers(1).echogram_layer_name = 'plummer';
layers(1).layerdata_layer_name = 'plummer';
layers(1).x = X;
layers(1).y = Y;
physical_constants;
layers(1).data = bt_elev/(c/2/sqrt(er_ice));

param.plot_en = false;
param.echogram_source = 'mvdr';
param.layerdata_source = 'layerData';
% save_sources: cell array of strings indicating which layer sources
%   should be updated (options are ops, layerData, and echogram)
param.save_sources = {'ops'};
param.save_mode = 'overwrite';
opsInsertLayerFromGrid(param,layers);

