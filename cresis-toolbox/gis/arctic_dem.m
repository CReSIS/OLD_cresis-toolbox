% Files are stored with names like:
% http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/33_34/33_34_1_1_5m_v2.0.tar
%
% To determine the files that you want, load the natural earth geotiff
% Then load the the arctic dem: ArcticDEM_Tile_Index_2017sept06.shp
% This will show all the tiles and you can get the range of tiles
% from there. Then specify the range in lat_tiles and lon_tiles.
%
% File name format: lat_lon_lonsub_latsub:
%   Ignore sub indices (all subtiles will be downloaded for each tile)
%
% lon_tile 33,1: -800000 (left side)
% lat_tile 31,1: -1000000 (lower side)
% Subtile size: 50000
% Tile size: 100000
% lat origin (tile 0): -4100000
% lon origin (tile 0): -4100000
% resolution: 5

wget_enable = 1; % wget ArcticDEM tiles
untar_enable = 1; % untar the tiles
decimate_enable = 1; % decimate the tiles by M
M = 8;
mosaic_enable = 1; % mosaic the decimated tiles

if 0
  % 20140325_0[5-7]
  lat_tiles = [32:33];
  lon_tiles = [31:33];
  mosaic_fn = '2014_Greenland_P3_20140325.tif';
end

if 0
  % 20140401_03
  lat_tiles = [31:34];
  lon_tiles = [33:38];
  mosaic_fn = '2014_Greenland_P3_20140401_03.tif';
end

% 20140506_01
if 0
  lat_tiles = [28:31];
  lon_tiles = [31:35];
  mosaic_fn = '2014_Greenland_P3_20140506_01.tif';
end

% Iceland
if 1
  lat_tiles = [14:19];
  lon_tiles = [50:55];
  mosaic_fn = 'Iceland.tif';
end

%% Automated section
if mod(10000/M,1) ~= 0
  factor(10000)
  error('M must be a factor of 10000');
end
cd('/cresis/snfs1/dataproducts/GIS_data/arctic/ArcticDEM');

if mosaic_enable
  Nx = 10000/M * 2 * length(lon_tiles);
  Ny = 10000/M * 2 * length(lat_tiles);
  mosaic_RGB = NaN*zeros(Ny,Nx,'single');
  mosaic_R = zeros(3,2,'double');
  mosaic_R(3,1) = -4100000 + lon_tiles(1)*100000 - 5*M/2;
  mosaic_R(3,2) = -4100000 + lat_tiles(1)*100000 - 5*M/2;
  mosaic_R(1,2) = 5*M;
  mosaic_R(2,1) = 5*M;
  mosaic_proj = [];
end

for lat_tiles_idx = 1:length(lat_tiles)
  for lon_tiles_idx = 1:length(lon_tiles)
    for sub_lon_idx = 1:2
      for sub_lat_idx = 1:2
        tif_fn = sprintf('%02d_%02d_%01d_%01d_5m_v2.0.tar', ...
          lat_tiles(lat_tiles_idx), lon_tiles(lon_tiles_idx), sub_lon_idx, sub_lat_idx);
        geotiff_fn = sprintf('%02d_%02d_%01d_%01d_5m_v2.0_reg_dem.tif', ...
          lat_tiles(lat_tiles_idx), lon_tiles(lon_tiles_idx), sub_lon_idx, sub_lat_idx);
        geotiff_dec_fn = sprintf('%02d_%02d_%01d_%01d_5m_v2.0_reg_dem_dec%d.tif', ...
          lat_tiles(lat_tiles_idx), lon_tiles(lon_tiles_idx), sub_lon_idx, sub_lat_idx, M);
        
        fprintf('%s: %s\n', datestr(now), tif_fn);
        
        if wget_enable
          if ~exist(tif_fn)
            cmd = sprintf('wget http://data.pgc.umn.edu/elev/dem/setsm/ArcticDEM/mosaic/v2.0/%02d_%02d/%s', ...
              lat_tiles(lat_tiles_idx), lon_tiles(lon_tiles_idx), tif_fn);
            fprintf('  %s\n', cmd);
            system(cmd);
          end
        end
        
        if untar_enable
          geotiff_fn = sprintf('%02d_%02d_%01d_%01d_5m_v2.0_reg_dem.tif', ...
            lat_tiles(lat_tiles_idx), lon_tiles(lon_tiles_idx), sub_lon_idx, sub_lat_idx);
          if ~exist(geotiff_fn,'file')
            if exist(tif_fn,'file')
              fprintf('  Untar\n');
              untar(tif_fn);
            end
          end
        end
        
        if decimate_enable
          if exist(geotiff_fn,'file') %&& ~exist(geotiff_dec_fn,'file')
            fprintf('  Reading %s\n', geotiff_fn);
            proj = geotiffinfo(geotiff_fn);
            
            % Read the image
            [RGB, R, tmp] = geotiffread(geotiff_fn);
            % Change no data to NaN
            RGB(RGB==-9999) = NaN;
            
            if 0
              mapshow(double(RGB), R);
              imagesc(RGB); colormap(jet(256));
              imagesc( R(3,1) + R(2,1)*(1:size(RGB,2)), R(3,2) + R(1,2)*(1:size(RGB,1)), RGB);
              
              % Corner coordinates
              R(3,1) + R(2,1)
              R(3,2) + R(1,2)
              R(3,1) + R(2,1)*size(RGB,2)
              R(3,2) + R(1,2)*size(RGB,1)
              
              % Tile size
              X_tile = R(2,1)*size(RGB,2)
              Y_tile = R(1,2)*size(RGB,1)
            end
            
            % Decimate Tile
            fprintf('  Decimate\n');
            RGB = nan_fir_dec(nan_fir_dec(RGB,M).',M).';
            
            R(3,1) = R(3,1) + R(2,1)*(M+1)/2 - R(2,1)*M;
            R(3,2) = R(3,2) + R(1,2)*(M+1)/2 - R(1,2)*M;
            R(2,1) = R(2,1) * M;
            R(1,2) = R(1,2) * M;
            
            fprintf('  Writing %s\n', geotiff_dec_fn);
            geotiffwrite(geotiff_dec_fn, RGB, R,  ...
              'GeoKeyDirectoryTag', proj.GeoTIFFTags.GeoKeyDirectoryTag);
          end
        end
        
        if mosaic_enable
          if exist(geotiff_dec_fn,'file')
            fprintf('  Reading %s\n', geotiff_dec_fn);
            mosaic_proj = geotiffinfo(geotiff_dec_fn);
            lon_idxs = 2*10000/M*(lon_tiles_idx-1) + 10000/M*(sub_lon_idx-1) + (1:10000/M);
            lat_idxs = 2*10000/M*(lat_tiles_idx-1) + 10000/M*(sub_lat_idx-1) + (10000/M:-1:1);
            [mosaic_RGB(lat_idxs,lon_idxs), ~, ~] = geotiffread(geotiff_dec_fn);
            imagesc(mosaic_RGB)
            set(gca,'YDir','normal')
            drawnow;
          end
        end
        
      end
    end
  end
end

if mosaic_enable
  fprintf('  Writing %s\n', mosaic_fn);
  geotiffwrite(mosaic_fn, mosaic_RGB, mosaic_R,  ...
    'GeoKeyDirectoryTag', mosaic_proj.GeoTIFFTags.GeoKeyDirectoryTag);
  if 0
    plot_geotiff(mosaic_fn)
  end
end
