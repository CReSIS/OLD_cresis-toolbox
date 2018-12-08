% Class dem_class
%
% Class for loading and working with digital elevation models. Normally
% stored in the global variable gdem to avoid multiple loads of the GIS
% files.
%
% Examples:
% See run_dem_class.m
%
% Author: John Paden

classdef dem_class < handle
  properties
    % dem_info: Struct array of DEMs
    dem_info
    % res: Resolution (m)
    res
    
    % DEM properties
    dem_all
    x_all
    y_all
    
    % Region of interest fields
    % Truncated DEM structure: dem, x, y
    dem
    % ocean shape file structure: shp_all, bb_all, shp, bb
    ocean
    % mean sea level structure: lat, lon, elev
    msl
    
    % Geographic vector
    name % name of vector (usually day_seg)
    x % vector query point x-value, N elements
    y % vector query point y-value, N elements
    di % DEM index for each (x,y) vector query point, N elements
    lat % vector query point latitude-value, N elements
    lon % vector query point longitude-value, N elements
    min_x
    max_x
    min_y
    max_y
    min_lat
    max_lat
    mean_lon
    max_lon
    min_lon
  end
  
  methods
    %% constructor
    function obj = dem_class(res)
      
      % Input checks
      % ===================================================================
      if ~exist('res','var') || isempty(res)
        obj.res = 1000;
      end
      
      % Setup DEM List
      % ===================================================================
      obj.dem_info = [];
      
      % Arctic DEM
      % -------------------------------------------------------------------
      % http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/10m/11_41/11_41_10m_v3.0.tar.gz
      % http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/2m/48_63/48_63_1_1_2m_v3.0.tar.gz
      % http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/2m/48_63/48_63_2_2_2m_v3.0.tar.gz
      di = 1;
      obj.dem_info(di).url = 'http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v3.0/';
      obj.dem_info(di).res = [2 10 32 100 500 1000];
      obj.dem_info(di).res_str = {'2m' '10m' '32m' '100m' '500m' '1km'};
      obj.dem_info(di).tile_en = [1 1 1 0 0 0];
      obj.dem_info(di).subtile_en = [1 0 0 0 0 0];
      obj.dem_info(di).mosaic_fn_fh = @(res_str) sprintf('%s/arcticdem_mosaic_%s_v3.0.tif',res_str,res_str);
      obj.dem_info(di).tile_fn_fh = @(x,y,res_str) sprintf('%s/%d_%d/%d_%d_%s_v3.0.tar.gz',res_str,y,x,y,x,res_str);
      obj.dem_info(di).subtile_fn_fh = @(x,y,x_sub,y_sub,res_str) sprintf('%s/%d_%d/%d_%d_%d_%d_%s_v3.0.tar.gz',res_str,y,x,y,x,y_sub,x_sub,res_str);
      obj.dem_info(di).acknowledge = 'Geospatial support for this work provided by the Polar Geospatial Center under NSF-OPP awards 1043681 and 1559691. DEMs provided by the Polar Geospatial Center under NSF-OPP awards 1043681, 1559691, and 1542736.';
      obj.dem_info(di).citation = 'Porter, Claire; Morin, Paul; Howat, Ian; Noh, Myoung-Jon; Bates, Brian; Peterman, Kenneth; Keesey, Scott; Schlenk, Matthew; Gardiner, Judith; Tomko, Karen; Willis, Michael; Kelleher, Cole; Cloutier, Michael; Husby, Eric; Foga, Steven; Nakamura, Hitomi; Platson, Melisa; Wethington, Michael, Jr.; Williamson, Cathleen; Bauer, Gregory; Enos, Jeremy; Arnold, Galen; Kramer, William; Becker, Peter; Doshi, Abhijit; D?Souza, Cristelle; Cummens, Pat; Laurier, Fabien; Bojesen, Mikkel, 2018, "ArcticDEM", https://doi.org/10.7910/DVN/OHHUKH, Harvard Dataverse, V1, [Date Accessed]';
      obj.dem_info(di).y_tile_origin = 41;
      obj.dem_info(di).x_tile_origin = 41;
      obj.dem_info(di).y_tile_size = 100e3;
      obj.dem_info(di).x_tile_size = 100e3;
      obj.dem_info(di).y_subtile_size = 50e3;
      obj.dem_info(di).x_subtile_size = 50e3;
      obj.dem_info(di).no_data = -9999;
      obj.dem_info(di).proj = geotiffinfo(ct_filename_gis([],'greenland/DEM/GIMP/gimpdem_90m.tif'));
      obj.dem_info(di).out_path = '/cresis/snfs1/dataproducts/GIS_data/arctic/ArcticDEM/';
      
      % REMA
      % -------------------------------------------------------------------
      % http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/8m/09_38/09_38_8m.tar.gz
      di = 2;
      obj.dem_info(di).url = 'http://data.pgc.umn.edu/elev/dem/setsm/REMA/mosaic/v1.1/';
      obj.dem_info(di).res = [8 100 200 1000];
      obj.dem_info(di).res_str = {'8m' '100m' '200m' '1000m'};
      obj.dem_info(di).subtile_en = [0 0 0 0];
      obj.dem_info(di).tile_fn_fh = @(north,east,res_str) sprintf('%s/%d_%d/%d_%d_%s_v3.0.tar.gz',res_str,north,east,north,east,res_str);
      obj.dem_info(di).acknowledge = 'Geospatial support for this work provided by the Polar Geospatial Center under NSF-OPP awards 1043681 and 1559691. DEMs provided by the Byrd Polar and Climate Research Center and the Polar Geospatial Center under NSF-OPP awards 1543501, 1810976, 1542736, 1559691, 1043681, 1541332, 0753663, 1548562, 1238993 and NASA award NNX10AN61G. Computer time provided through a Blue Waters Innovation Initiative. DEMs produced using data from DigitalGlobe, Inc.';
      obj.dem_info(di).citation = 'Howat, Ian; Morin, Paul; Porter, Claire; Noh, Myong-Jong, 2018, "The Reference Elevation Model of Antarctica", https://doi.org/10.7910/DVN/SAIK8B, Harvard Dataverse, V1';
      obj.dem_info(di).y_tile_origin = 31;
      obj.dem_info(di).x_tile_origin = 31;
      obj.dem_info(di).y_tile_size = 100e3;
      obj.dem_info(di).x_tile_size = 100e3;
      obj.dem_info(di).no_data = -9999;
      obj.dem_info(di).proj = geotiffinfo(ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif'));
      
      % Prepare DEM fields
      % -------------------------------------------------------------------
      obj.clear();
      obj.reset();
    end
    
    %% destructor
    % =====================================================================
    function delete(obj)
    end
    
    %% set_vector: Set Geographic Vector
    % =====================================================================
    function obj = set_vector(obj,lat,lon,name)
      
      % Update class fields and determine which DEM each point belongs to
      % -------------------------------------------------------------------
      if exist('name','var')
        obj.name = name;
      end
      if isempty(lat)
        return;
      else
        obj.lat = lat;
        obj.lon = lon;
        obj.x = nan(size(lat));
        obj.y = nan(size(lat));
        obj.di = nan(size(lat));
        % Northern hemisphere
        obj.di(lat>0) = 1;
        % Southern hemisphere
        obj.di(lat<0) = 2;
        
        for di = unique(obj.di)
          mask = obj.di == di;
          [obj.x(mask),obj.y(mask)] = projfwd(obj.dem_info(di).proj,lat(mask),lon(mask));
        end
      end
      
      % Determine bounds of region of interest
      % -------------------------------------------------------------------
      obj.min_x = min(obj.x);
      obj.max_x = max(obj.x);
      obj.min_y = min(obj.y);
      obj.max_y = max(obj.y);
      
      obj.min_lat = min(obj.lat);
      obj.max_lat = max(obj.lat);
      % Handle longitude in a special way because it wraps around.
      obj.mean_lon = angle(mean(exp(1i*obj.lon/180*pi)))*180/pi;
      obj.max_lon = obj.mean_lon + max(angle(exp(1i*(obj.lon-obj.mean_lon)/180*pi)))*180/pi;
      obj.min_lon = obj.mean_lon + min(angle(exp(1i*(obj.lon-obj.mean_lon)/180*pi)))*180/pi;
      
      % Load GIS data as needed
      % -------------------------------------------------------------------
      
      % If msl is not big enough, it is reloaded to fit vector
      obj.load_msl();
      
      % If ocean is not big enough or has not been truncated yet, a new
      % truncated version is created
      obj.load_ocean();
      
      % If the whole region fits into a single tile, then this tile is
      % loaded and truncated to the region of interest
      obj.load_dem();
    end
    
    %% load_msl: Load MSL
    % =====================================================================
    function obj = load_msl(obj)
      if isempty(obj.x)
        return;
      end
      if ~isempty(obj.msl)
        if obj.msl.lon(1) <= obj.min_lon && obj.msl.lon(end) >= obj.max_lon ...
            && obj.msl.lat(1) <= obj.min_lat && obj.msl.lat(end) >= obj.max_lat
        end
      end
      
      % Mean Sea Level (msl)
      % -------------------------------------------------------------------
      if 0
        % EGM-96
        obj.msl.fn = ct_filename_gis([],'world\egm96_geoid\WW15MGH.DAC');
        points = [];
        [obj.msl.lat,obj.msl.lon,obj.msl.elev] = egm96_loader(obj.msl.fn);
        obj.msl.lon = [obj.msl.lon 360];
        obj.msl.elev = [obj.msl.elev obj.msl.elev(:,1)];
        [obj.msl.lon,obj.msl.lat] = meshgrid(obj.msl.lon,obj.msl.lat);
      else
        % Load DTU mean sea level
        obj.msl.fn = ct_filename_gis([],fullfile('world','dtu_meansealevel','DTU10MSS_1min.nc'));
        obj.msl.lat = ncread(obj.msl.fn,'lat');
        obj.msl.lon = ncread(obj.msl.fn,'lon');
        dlat = obj.msl.lat(2)-obj.msl.lat(1);
        lat_idxs = find(obj.msl.lat >= obj.min_lat-2*dlat & obj.msl.lat <= obj.max_lat+2*dlat);
        dlon = obj.msl.lon(2)-obj.msl.lon(1);
        rel_lon = obj.mean_lon + angle(exp(1i*(obj.msl.lon - obj.mean_lon)/180*pi))*180/pi;
        lon_idxs = find(rel_lon >= obj.min_lon-2*dlon & rel_lon <= obj.max_lon+2*dlon);
        break_idx = find(diff(lon_idxs)~=1);
        obj.msl.lat = obj.msl.lat(lat_idxs);
        obj.msl.lon = rel_lon(lon_idxs);
        % Transpose elev because "x" axis (which is longitude) must be on the
        % column dimension for interp2.
        % Convert to single because interp2 requires single or double type
        % and single is smaller yet has enough precision.
        if isempty(break_idx)
          obj.msl.elev = single(ncread(obj.msl.fn,'mss', ...
            [lon_idxs(1) lat_idxs(1)],[length(lon_idxs) length(lat_idxs)]).');
        else
          obj.msl.elev = single(ncread(obj.msl.fn,'mss', ...
            [lon_idxs(break_idx+1) lat_idxs(1)],[length(lon_idxs)-break_idx length(lat_idxs)]).');
          obj.msl.elev = [obj.msl.elev, single(ncread(obj.msl.fn,'mss', ...
            [1 lat_idxs(1)],[break_idx length(lat_idxs)]).')];
          obj.msl.lon = obj.msl.lon([break_idx+1:end,1:break_idx]);
        end
        [obj.msl.lon,unique_idxs] = unique(obj.msl.lon);
        obj.msl.elev = obj.msl.elev(:,unique_idxs);
      end
    end
    
    %% load_ocean: Load ocean mask
    % =====================================================================
    function obj = load_ocean(obj)
      if isempty(obj.x)
        return;
      end
      
      % Load ocean mask shape file (-180 to +180 lon)
      % -------------------------------------------------------------------
      if ~isfield(obj.ocean,'shp_all')
        ocean_mask_fn = ct_filename_gis([],fullfile('world','land_mask','Land_Mask_IDL_jharbeck','GSHHS_f_L1.shp'));
        warning off;
        obj.ocean.shp_all = shaperead(ocean_mask_fn);
        warning on;
        obj.ocean.bb_all = [obj.ocean.shp_all(:).BoundingBox];
      end
      
      % Get all ocean shapes within the data segment bounding box. Shape is
      % not included if any of the following holds:
      %  - Bottom of the shape is above the top of the segment, >max_lat
      %  - Top of the shape is below the bottom of the segment, <min_lat
      %  - Left side of the shape is to the right of the segment, >max_lon
      %  - Right side of the shape is to the left of the segment, <min_lon
      % Handle longitude in a special way because it wraps around.
      if isempty(obj.ocean.bb_all)
        obj.ocean.shp = [];
      else
        rel_min_lon = obj.mean_lon + angle(exp(1i*(obj.ocean.bb_all(1,1:2:end) - obj.mean_lon)/180*pi))*180/pi;
        rel_max_lon = obj.mean_lon + angle(exp(1i*(obj.ocean.bb_all(2,1:2:end) - obj.mean_lon)/180*pi))*180/pi;
        bb_good_mask = ~(obj.ocean.bb_all(1,2:2:end)>obj.max_lat | obj.ocean.bb_all(2,2:2:end)<obj.min_lat ...
          | rel_min_lon>obj.max_lon | rel_max_lon<obj.min_lon);
        obj.ocean.shp = obj.ocean.shp_all(bb_good_mask);
        % All bounding boxes of every shape
        obj.ocean.bb = [obj.ocean.shp(:).BoundingBox];
      end
      
      if 0
        % Debug code to check bounding box code
        figure(1); clf;
        for idx=1:length(obj.ocean.shp)
          if length(obj.ocean.shp(idx).X) > 2000
            plot(obj.ocean.shp(idx).X(1:5:end),obj.ocean.shp(idx).Y(1:5:end))
          else
            plot(obj.ocean.shp(idx).X,obj.ocean.shp(idx).Y)
          end
          hold on;
        end
        plot(obj.lon(1:100:end),obj.lat(1:100:end),'k.')
      end
      
    end
    
    %% load_dem: Load land DEM
    % =====================================================================
    function obj = load_dem(obj)

      for di = unique(obj.di)
        mask = obj.di == di;
        if all(~mask)
          continue
        end
        ri = obj.get_ri(di);
        
        if ~obj.dem_info(di).tile_en(ri)
          if ~isempty(obj.dem.x{di}) ...
              && obj.min_x >= min(obj.dem.x{di}) && obj.max_x <= max(obj.dem.x{di}) ...
              && obj.min_y >= min(obj.dem.y{di}) && obj.max_y <= max(obj.dem.y{di})
            continue;
          end
          
          % Single large DEM file
          % ---------------------------------------------------------------
          url = fullfile(obj.dem_info(di).url, obj.dem_info(di).mosaic_fn_fh(obj.dem_info(di).res_str{ri}));
          [~,url_name,url_ext] = fileparts(url);
          tif_fn = fullfile(obj.dem_info(di).out_path, [url_name,url_ext]);
          
          if ~exist(tif_fn,'file')
            cmd = sprintf('wget -P %s %s', obj.dem_info(di).out_path, url);
            fprintf('  %s\n', cmd);
            system(cmd);
            if ~exist(tif_fn,'file')
              error('Failed to download file.');
            end
          end
          
          if isempty(obj.x_all{di})
            % Load the DEM file
            % ---------------------------------------------------------------
            % Adapted from NSIDC website:
            % Retrieve the image info
            GeoKeys = imfinfo(tif_fn);
            N = GeoKeys.Width;
            M = GeoKeys.Height;
            dx = GeoKeys.ModelPixelScaleTag(1);
            dy = GeoKeys.ModelPixelScaleTag(2);
            minx = GeoKeys.ModelTiepointTag(4);
            maxy = GeoKeys.ModelTiepointTag(5);
            
            % Generate obj.x/obj.y pixel location vectors
            obj.x_all{di} = minx + dx/2 + ((0:N-1).*dx);
            obj.y_all{di} = maxy - dy/2 - ((M -1:-1:0).*dy);
            
            % Read the image data, flip upside down if necessary
            obj.dem_all{di} = imread(tif_fn);
            
            if GeoKeys.Orientation
              obj.y_all{di} = obj.y_all{di}(end:-1:1);
            end
          end
          
          % Truncate the DEM to region of interest
          % ---------------------------------------------------------------
          dx = abs(obj.x_all{di}(2) - obj.x_all{di}(1));
          dy = abs(obj.y_all{di}(2) - obj.y_all{di}(1));
          x_mask = obj.x_all{di} >= obj.min_x-2*dx & obj.x_all{di} <= obj.max_x+2*dx;
          y_mask = obj.y_all{di} >= obj.min_y-2*dy & obj.y_all{di} <= obj.max_y+2*dy;
          
          obj.dem.dem{di} = obj.dem_all{di}(y_mask,x_mask);
          obj.dem.dem{di}(obj.dem.dem{di}==obj.dem_info(di).no_data) = NaN;
          obj.dem.x{di} = obj.x_all{di}(x_mask);
          obj.dem.y{di} = obj.y_all{di}(y_mask);
          
          if 0
            figure(1); clf;
            imagesc(obj.dem.x{di}/1e3,obj.dem.y{di}/1e3,obj.dem.dem{di});
            set(gca,'YDir','normal')
            hold on;
            plot(obj.x/1e3,obj.y/1e3,'k','LineWidth',2);
            cc = caxis;
          end
        end
      end
    end
    
    %% get_vector_dem: Get DEM corresponding to geographic vector
    % =====================================================================
    function [land_dem,msl,ocean_mask] = get_vector_dem(obj)
      %  When a tile is accessed and it does not exist, then it is downloaded and
      %  loaded and stored
      
      % MSL
      % -------------------------------------------------------------------
      msl = interp2(obj.msl.lon,obj.msl.lat,obj.msl.elev,obj.lon,obj.lat);
      
      % DEM and ocean_mask
      % -------------------------------------------------------------------
      ocean_mask = true(size(obj.lat));
      for di = unique(obj.di)
        mask = obj.di == di;
        if all(~mask)
          continue
        end
        
        % DEM
        % -----------------------------------------------------------------
        ri = obj.get_ri(di);
        if ~obj.dem_info(di).tile_en(ri)
          % Single DEM file
          land_dem(mask) = interp2(obj.dem.x{di},obj.dem.y{di},obj.dem.dem{di},obj.x(mask),obj.y(mask));
          
        else
          % DEM is broken into multiple tile files
          land_dem = nan(size(obj.x));
          
          % Determine which tiles are needed
          tile_x = obj.dem_info(di).x_tile_origin + obj.x/obj.dem_info(di).x_tile_size;
          tile_y = obj.dem_info(di).y_tile_origin + obj.y/obj.dem_info(di).y_tile_size;
          tiles = unique(floor([tile_x.' tile_y.']),'rows');
          
          % Download and untar files
          % ===============================================================
          for tiles_idx = 1:size(tiles,1)
            url_list = {};
            xi_list = [];
            yi_list = [];
            if obj.dem_info(di).subtile_en(ri)
              for sub_x_idx = 1:2
                for sub_y_idx = 1:2
                  if any(tile_x >= tiles(tiles_idx,1) & tile_x <= tiles(tiles_idx,1)+1 ...
                      & tile_y >= tiles(tiles_idx,2) & tile_y <= tiles(tiles_idx,2)+1)
                    url_list{end+1} = fullfile(obj.dem_info(di).url,obj.dem_info(di).subtile_fn_fh(tiles(tiles_idx,1),tiles(tiles_idx,2),sub_x_idx,sub_y_idx,obj.dem_info(di).res_str{ri}));
                    xi_list(end+1) = 2*tiles(tiles_idx,1) + sub_x_idx - 2;
                    yi_list(end+1) = 2*tiles(tiles_idx,2) + sub_y_idx - 2;
                    x_tile_size = obj.dem_info(di).x_subtile_size;
                    y_tile_size = obj.dem_info(di).y_subtile_size;
                    x_tile_origin = obj.dem_info(di).x_tile_origin - 1;
                    y_tile_origin = obj.dem_info(di).y_tile_origin - 1;
                  end
                end
              end
            else
              url_list{end+1} = fullfile(obj.dem_info(di).url,obj.dem_info(di).tile_fn_fh(tiles(tiles_idx,1),tiles(tiles_idx,2),obj.dem_info(di).res_str{ri}));
              xi_list(end+1) = tiles(tiles_idx,1);
              yi_list(end+1) = tiles(tiles_idx,2);
              x_tile_size = obj.dem_info(di).x_tile_size;
              y_tile_size = obj.dem_info(di).y_tile_size;
              x_tile_origin = obj.dem_info(di).x_tile_origin;
              y_tile_origin = obj.dem_info(di).y_tile_origin;
            end
            
            for url_idx = 1:length(url_list)
              url = url_list{url_idx};
              xi = xi_list(url_idx);
              yi = yi_list(url_idx);
              fprintf('%s: %s\n', datestr(now), url);
              
              [~,url_name,url_ext1] = fileparts(url);
              [~,url_name,url_ext2] = fileparts(url_name);
              fn = fullfile(obj.dem_info(di).out_path,[url_name,url_ext2,url_ext1]);
              tif_fn = fullfile(obj.dem_info(di).out_path,[url_name '_reg_dem.tif']);
              
              if ~exist(fn,'file') && ~exist(tif_fn,'file')
                cmd = sprintf('wget -P %s %s', obj.dem_info(di).out_path, url);
                fprintf('  %s\n', cmd);
                system(cmd);
                if ~exist(fn,'file')
                  error('Failed to download file.');
                end
              end
              
              if ~exist(tif_fn,'file')
                fprintf('Untar %s\n  %s\n', fn, obj.dem_info(di).out_path);
                untar(fn, obj.dem_info(di).out_path);
                if ~exist(tif_fn,'file')
                  error('Failed to untar file or file not found in tar archive.');
                end
              end
              
              if exist(fn,'file')
                delete(fn);
              end
              
              if size(obj.dem_all{di},1) < xi || size(obj.dem_all{di},2) < yi ...
                  || isempty(obj.dem_all{di}{xi,yi})
                % Adapted from NSIDC website:
                % Retrieve the image info
                GeoKeys = imfinfo(tif_fn);
                N = GeoKeys.Width;
                M = GeoKeys.Height;
                dx = GeoKeys.ModelPixelScaleTag(1);
                dy = GeoKeys.ModelPixelScaleTag(2);
                minx = GeoKeys.ModelTiepointTag(4);
                maxy = GeoKeys.ModelTiepointTag(5);
                
                % Generate obj.x/obj.y pixel location vectors
                obj.x_all{di}{xi,yi} = minx + dx/2 + ((0:N-1).*dx);
                obj.y_all{di}{xi,yi} = maxy - dy/2 - ((M -1:-1:0).*dy);
                
                % Read the image data
                obj.dem_all{di}{xi,yi} = imread(tif_fn);
                
                % Set bad values to NaN
                obj.dem_all{di}{xi,yi}(obj.dem_all{di}{xi,yi}==obj.dem_info(di).no_data) = NaN;
                
                % Flip upside down if necessary
                if GeoKeys.Orientation
                  obj.y_all{di}{xi,yi} ...
                    = obj.y_all{di}{xi,yi}(end:-1:1);
                end
              end
              
              x_min_tile = (xi-x_tile_origin) * x_tile_size;
              x_max_tile = x_min_tile + x_tile_size;
              y_min_tile = (yi-y_tile_origin) * y_tile_size;
              y_max_tile = y_min_tile + y_tile_size;
              tile_mask = obj.x >= x_min_tile & obj.x < x_max_tile ...
                & obj.y >= y_min_tile & obj.y < y_max_tile;
              land_dem(tile_mask) = interp2(obj.x_all{di}{xi,yi}, ...
                obj.y_all{di}{xi,yi}, obj.dem_all{di}{xi,yi}, obj.x(tile_mask), obj.y(tile_mask));
              
            end
          end
        end
        
        % Ocean Mask
        % -----------------------------------------------------------------
        min_x = min(obj.x);
        max_x = max(obj.x);
        min_y = min(obj.y);
        max_y = max(obj.y);
        
        min_lat = min(obj.lat(mask));
        max_lat = max(obj.lat(mask));
        % Handle longitude in a special way because it wraps around.
        mean_lon = angle(mean(exp(1i*obj.lon(mask)/180*pi)))*180/pi;
        max_lon = mean_lon + max(angle(exp(1i*(obj.lon(mask)-mean_lon)/180*pi)))*180/pi;
        min_lon = mean_lon + min(angle(exp(1i*(obj.lon(mask)-mean_lon)/180*pi)))*180/pi;
        
        % Restrict ocean mask features to our dataset (i.e. mask all features
        % whose bounding boxes fall outside our limits.
        if isempty(obj.ocean.shp)
          ocean_shp = [];
        else
          rel_min_lon = mean_lon + angle(exp(1i*(obj.ocean.bb(1,1:2:end) - mean_lon)/180*pi))*180/pi;
          rel_max_lon = mean_lon + angle(exp(1i*(obj.ocean.bb(2,1:2:end) - mean_lon)/180*pi))*180/pi;
          bb_good_mask = ~(obj.ocean.bb(1,2:2:end)>max_lat | obj.ocean.bb(2,2:2:end)<min_lat ...
            | rel_min_lon>max_lon | rel_max_lon<min_lon);
          ocean_shp = obj.ocean.shp(bb_good_mask);
        end
        
        % Create polygons, poly_x/poly_y, with all ocean shape features which
        % lie in the bounding box.
        % Further restrict the polygons by checking for bounding box overlap in
        % projected coordinates
        poly_x = cell(0);
        poly_y = cell(0);
        for shp_idx = 1:length(ocean_shp)
          % convert polygon to projected coordinates
          [x,y] = projfwd(obj.dem_info(di).proj,ocean_shp(shp_idx).Y,ocean_shp(shp_idx).X);
          % if polygon is within projected bounding box
          if min(x) < max_x && max(x) > min_x ...
              && min(y) < max_y && max(y)>min_y
            % add polygon
            poly_x{end+1} = [x,nan];
            poly_y{end+1} = [y,nan];
          end
        end
        
        % Create ocean mask to determine which points lie in the ocean
        mask_idxs = find(mask);
        for poly_idx = 1:length(poly_x)
          
          % Mask showing which DEM points are in polygon (on land)
          land_mask_tmp = inpolygon(obj.x(mask),obj.y(mask),[poly_x{poly_idx}(1:100:end)],[poly_y{poly_idx}(1:100:end)]);
          ocean_mask(mask_idxs(land_mask_tmp)) = false;
        end
        
        % Debug
        % -----------------------------------------------------------------
        if 0
          x_out = obj.x(mask);
          y_out = obj.y(mask);
          dem_out = land_dem(mask);
          dec_idxs = round(linspace(1,length(obj.x(mask)),1000));
          mask_dec_idxs = logical(interp1(1:length(obj.x(mask)), single(ocean_mask(mask)), dec_idxs, 'nearest', 'extrap'));
          along_track = geodetic_to_along_track(obj.lat(mask),obj.lon(mask));
          
          figure(2); clf;
          scatter(along_track(dec_idxs),dem_out(dec_idxs),[],dem_out(dec_idxs));
          hold on;
          scatter(along_track(dec_idxs),msl(dec_idxs),[],msl(dec_idxs));
          plot(along_track(dec_idxs(mask_dec_idxs)),msl(dec_idxs(mask_dec_idxs)),'k.');
          caxis([min(obj.dem.dem{di}(:)) max(obj.dem.dem{di}(:))]);
          
          figure(3); clf;
          scatter(x_out(dec_idxs),y_out(dec_idxs),[],dem_out(dec_idxs));
          hold on;
          plot(x_out(dec_idxs(mask_dec_idxs)),y_out(dec_idxs(mask_dec_idxs)),'k.');
          caxis([min(obj.dem.dem{di}(:)) max(obj.dem.dem{di}(:))]);
        end
      end
      
    end
        
    %% clear: Clear all loaded DEMs from memory
    % =====================================================================
    function clear(obj)
      obj.dem_all = cell(size(obj.dem_info));
      obj.x_all = cell(size(obj.dem_info));
      obj.y_all = cell(size(obj.dem_info));
    end
    
    %% reset: Reset to load a new segment
    % =====================================================================
    function reset(obj)
      if ~isempty(obj.di)
        for di = unique(obj.di)
          mask = obj.di == di;
          if all(~mask)
            continue
          end
          ri = obj.get_ri(di);
          
          if obj.dem_info(di).tile_en(ri)
            % Tiles were used, so clear all loaded DEMs since it is likely that
            % they will not be reused
            obj.clear();
          end
        end
      end
      
      % Geographic vector
      obj.di = [];
      obj.x = [];
      obj.y = [];
      obj.lat = [];
      obj.lon = [];
      
      % Region of interest fields
      obj.name = [];
      obj.dem.dem = cell(size(obj.dem_info));
      obj.dem.x = cell(size(obj.dem_info));
      obj.dem.y = cell(size(obj.dem_info));
      obj.msl = [];
      obj.ocean.shp = [];
      obj.ocean.bb = [];
    end
    
    %% set_res: Set resolution
    % =====================================================================
    function set_res(obj,res)
      if res ~= obj.res
        obj.clear();
        obj.reset();
        obj.res = res;
      end
    end
    
    %% get_ri: Get resolution index
    % =====================================================================
    function ri = get_ri(obj,di)
      ri = find(obj.dem_info(di).res <= obj.res,1,'last');
      if isempty(ri)
        warning('All resolutions are coarser obj.dem_info(di).res = %g m than the desired obj.param_dem_class.res = %g m. Downloading the finest resolution available.', obj.dem_info(ri).res(1), obj.param_dem_class.res);
        ri = 1;
      end
    end
    
  end
  
  
end
