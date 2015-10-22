% Script to seek the appropriate ASCAT sea ice backscatter file from the
% IFREMER FTP site and convert to a GeoTIFF for the flight maps in the
% processing output.
% Another possible source is: ftp://ftp.scp.byu.edu/data/ascat/2010/sir/msfa/Arc/111/a/msfa-a-Arc10-111-112.sir.gz
%   http://search.scp.byu.edu/
%
% GeoTiffs produced by ascat_to_geotiff are used in publish_map.m when
% param.type = ASCAT (set from run_load_data_by_gps_time). 

% 06/29/12
% Author: Ben Panzer, Trey Stafford

% =======================================================================
%% User Settings
% =======================================================================

param_fn = ct_filename_param('snow_param_2012_Greenland_P3.xls');
% param_fn = ct_filename_param('snow_param_2010_Antarctica_DC8.xls');
overwrite = true;

% =======================================================================
%% Automated Section
% =======================================================================

% Save current directory
initial_directory = pwd;

% Read user-defined parameter spreadsheet 
params = read_param_xls(param_fn,[],'post');

% Pre-allocate some variables. 
old_year = 0;
old_month = 0;
old_day = 0;

% Get ASCAT data for each segment marked by the generic column of the
% params spreadsheet. 
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~param.cmd.generic
    continue;
  end
  
  fprintf('Getting ASCAT data for %s\n', param.day_seg);
  
  % Seperate day_seg to year, month, day
  year = str2num(param.day_seg(1:4));
  month = str2num(param.day_seg(5:6));
  day = str2num(param.day_seg(7:8));

  % Check to see if geotiff was previously generated.
  if year == old_year && month == old_month && day == old_day
    fprintf('  Just grabbed this day... so skipping\n');
    continue;
  end
  old_year = year;
  old_month = month;
  old_day = day;
    
  % Define filenames and output path.
  fn_ftp = sprintf('ASCAT_%04d%02d%02d.nc.bz2',year,month,day);
  fn_local = sprintf('ASCAT_%04d%02d%02d.nc',year,month,day);
  output_fn_name = sprintf('ASCAT_%04d%02d%02d.tif',year,month,day);
 
  % Set the correct output path and cd.
  location = param.post.map.location;
  if strcmpi(param.post.map.location, 'Arctic')
    output_dir = fullfile(ct_filename_gis(param,''),'arctic','ASCAT')
  elseif strcmpi(param.post.map.location, 'Antarctica')
    output_dir = fullfile(ct_filename_gis(param,''),'antarctica','ASCAT')
  end
  cd(output_dir);
  
  fn = fullfile(output_dir,fn_local);
  
  % Check if the geotiff has already been created
  if exist(fullfile(output_dir,output_fn_name)) && overwrite == true ...
      || ~exist(fullfile(output_dir,output_fn_name))
    % load correct data (arctic vs antarctic) based on location from FTP site.
    fprintf('\n Loading ASCAT Data from FTP site...\n\n');
    if strcmpi(param.post.map.location, 'Arctic')
      command = sprintf('wget ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/psi-backscatter/data/ascat/arctic/netcdf/%04d/ASCAT_%04d%02d%02d.nc.bz2',year,year,month,day);
      system(command);
      system(['bunzip2 ', fn_ftp]);
      % Arctic: EPSG:3413, polar stereographic, -45 deg lon, -70 lat
      geo_tiff_params_fn = fullfile(ct_filename_gis(param,''),'arctic','ASCAT','params','arctic_geo_tiff_params.mat');
      load(geo_tiff_params_fn);
      is_arctic = 1;
    elseif strcmpi(param.post.map.location, 'Antarctica')
      command = sprintf('wget ftp://ftp.ifremer.fr/ifremer/cersat/products/gridded/psi-backscatter/data/ascat/antarctic/netcdf/%04d/ASCAT_%04d%02d%02d.nc.bz2',year,year,month,day);
      system(command);
      system(['bunzip2 ', fn_ftp]);
      % Antarctica: EPSG:3031, polar stereographic, 0 deg lon, -71 lat
      geo_tiff_params_fn = fullfile(ct_filename_gis(param,''),'antarctica','ASCAT','params','antarctic_geo_tiff_params.mat');
      load(geo_tiff_params_fn);
      proj.GeoTIFFTags.GeoKeyDirectoryTag.ProjNatOriginLatGeoKey = -71; % Match NSIDC
      proj.ProjParm(1) = -71; % Match NSIDC
      is_arctic = 0;
    else
      error('Location %s not supported',param.post.map.location);
    end

    % Load the backscatter data with mask and transpose
    sigma_not = ncread(fn,'sigma_40_mask');
    lat = ncread(fn,'latitude');
    lon = ncread(fn,'longitude');
    sigma_not = sigma_not';
    spatial_res = ncreadatt(fn,'/','spatial_resolution'); % should return 12.5 km
    spatial_res = str2num(spatial_res(1:4))*10^3;   % meters

    MY_sea_ice = -0.0034;

    fprintf('Processing Data...\n');

    % Define number of colors
    ncolors = 256;

    % Define the X and Y lims (NOTE: there is an error in lat/lon... the
    % bottom right hand corner is bad, so that we have to be careful
    % which x/y we choose)
    [x,y] = projfwd(proj,lat,lon);
    x = x([1 end],1);
    y = y(1,[1 end]);

    %
    img = zeros(size(sigma_not));
    img(isnan(sigma_not)) = 2; % Ocean
    img(sigma_not < -4.7e-4 & sigma_not > -4.9e-4) = 1; % Land
    % Scale the rest
    ice_mask = img==0;
    sigma_not = sigma_not - MY_sea_ice;
    sigma_not(sigma_not < 0) = 0;
    sigma_not = round(sigma_not / max(sigma_not(ice_mask)) * 253 + 3);
    img(ice_mask) = sigma_not(ice_mask);
    cmap(1,:) = [0.3906 0.3906 0.3906]; % land = gray
    cmap(2,:) = [0 0 0.5430];           % ocean = blue
    cmap(3:256,:) = colormap(jet(254)); % sea ice
    img = img - 1;
    
    imagesc(img);
    colormap(cmap);

    fprintf('Creating geotiff...\n');
    % Set the size of the raster after processing
    % Create the MapRasterReference
    R = spatialref.MapRasterReference;
    R.XLimWorld = [x(1) x(2)];
    R.YLimWorld = [y(2) y(1)];
    R.RasterSize = size(img);
    R.ColumnsStartFrom = 'north';
    R.RowsStartFrom = 'west';
    R.RasterInterpretation = 'cells';

    % Write the geotiff
    if exist(fullfile(output_dir,output_fn_name)) == 2
      command = sprintf('rm %s',fullfile(output_dir,output_fn_name));
      system(command);
    end

    output_fn = fullfile(output_dir,output_fn_name);
    fprintf('Writing geotiff output %s\n', output_fn);
    if strcmpi(param.post.map.location, 'Arctic')
      warning('Known bug in Matlab for Arctic ASCAT geotiff, proj.GeoTIFFTags.GeoKeyDirectoryTag.ProjNatOriginLatGeoKey is not written/read properly. To use this geotiff you must load the projection object and then set ProjParm(1) = 70.');
    end
    geotiffwrite(output_fn,ind2rgb(uint8(img),cmap), R,...
      'GeoKeyDirectoryTag', proj.GeoTIFFTags.GeoKeyDirectoryTag);

    if strcmpi(param.post.map.location, 'Arctic')
      warning('Demonstration of Arctic ASCAT geotiff error\n');
      proj2 = geotiffinfo(output_fn);
      proj2.ProjParm(1)
      [x,y] = projfwd(proj,69,-135)
      fprintf('  It does not match now:\n');
      [x,y] = projfwd(proj2,69,-135)
      fprintf('  It does match now:\n');
      proj2.ProjParm(1) = 70;
      [x,y] = projfwd(proj,69,-135)
    end
    
    % Clean up temporary netcdf files.
    if exist(fullfile(output_dir,fn_local))
      command = sprintf('rm %s',fullfile(output_dir,fn_local));
      system(command);
    end
    if exist(fullfile(output_dir,fn_ftp))
      command = sprintf('rm %s',fullfile(output_dir,fn_ftp));
      system(command);
    end
  else
    fprintf('Geotiff for %s already exists...skipping.\n', param.day_seg);
  end
  
end

% Re-set directory
cd(initial_directory);

fprintf('Done!\n');
return
