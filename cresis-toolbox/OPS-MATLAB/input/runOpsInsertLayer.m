% script runOpsInsertLayer.m
%
% Example script for running opsInsertLayer.m. Demonstrates a few of the
% most common operations to be performed with opsInsertLayer.
%
% Authors: John Paden

%% Example Operations (just choose one)

if 0
  %% Example 1: Gogineni's Jakobshavn points
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('rds_param_2008_Greenland_TO.xls'));
  
  points_fn = '/cresis/snfs1/scratch/paden/mass_conservation/Jakobshavn_2006_2009_Composite/Jakobshavn_2006_2009_Composite/flightlines/Jakobshavn_2006_2009_Composite_Flightlines.txt';
  
  % Load CSV points file
  fid = fopen(points_fn,'r');
  first_line = fgets(fid);
  headers = textscan(first_line,'%s','Delimiter',',');
  data = textscan(fid,'%f%f%f%s%s%f%f%f%f%f%f%f%f%f%f','Delimiter',',');
  fclose(fid);
  points = [];
  for header_idx = 1:length(headers{1})
    points.(headers{1}{header_idx}) = data{header_idx};
  end
  
  % Create geotiff projection structure
  % 1. Grab from a geotiff file with the same projection
  insert_param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));
  [points.x,points.y] = projfwd(insert_param.proj,points.LAT,points.LON);
  
  insert_param.eval.ref_source.name = 'surface';
  insert_param.eval.ref_source.source = 'ops';
  insert_param.eval.ref_gaps_fill.method = 'interp_finite';
  insert_param.eval.cmd = 'source = ref.twtt + source;';
  insert_param.x = points.x;
  insert_param.y = points.y;
  insert_param.data = (points.A_SURF-points.A_BED) / (c/2/sqrt(er_ice));;
  
  insert_param.type = 'point'; % Point data
  insert_param.layer_dest.name = 'gogineni2014_pnt';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'Gogineni JofG 2014 grid'; % For OPS layer_dest source
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'preserve_gaps';
  insert_param.gaps_fill.method_args = [300 60];
  opsInsertLayer(params, insert_param);
  
elseif 0
  %% Example 2: ATM Ramp Pass Data
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('kuband_param_2014_Greenland_P3.xls'));
  
  points_fn = '/scratch/metadata/2014_Greenland_P3/130407_truk_l12_kangramp_anthtrem_nopark.txt';
  [points.lat,points.lon,points.elev] = read_ramp_pass(points_fn);
  
  % Create geotiff projection structure
  % 1. Grab from a geotiff file with the same projection
  insert_param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));
  
  [points.x,points.y] = projfwd(insert_param.proj,points.lat,points.lon);
  
  insert_param.eval.ref_source.name = 'surface';
  insert_param.eval.ref_source.source = 'ops';
  insert_param.eval.ref_gaps_fill.method = 'interp_finite';
  insert_param.eval.cmd = 'source = (elev - source)/(c/2);';
  insert_param.x = points.x;
  insert_param.y = points.y;
  insert_param.data = points.elev;
  
  insert_param.type = 'point'; % Point data
  insert_param.layer_dest.name = 'kanger_ramp';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'ATM ramp pass Kangerlussuaq 2014'; % For OPS layer_dest source
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'preserve_gaps';
  insert_param.gaps_fill.method_args = [300 60];
  opsInsertLayer(params, insert_param);
  
elseif 0
  %% Example 3: Insert 3-D data
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'));
  
  load('jakob_frm1.mat','points');
  frm2 = load('jakob_frm2.mat','points');
  points.lat = cat(2,points.lat,frm2.points.lat);
  points.lon = cat(2,points.lon,frm2.points.lon);
  points.elev = cat(2,points.elev,frm2.points.elev);
  frm3 = load('jakob_2009frm1.mat','points');
  points.lat = cat(2,points.lat,frm3.points.lat);
  points.lon = cat(2,points.lon,frm3.points.lon);
  points.elev = cat(2,points.elev,frm3.points.elev);
  
  % Create geotiff projection structure
  % 1. Grab from a geotiff file with the same projection
  insert_param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));
  
  [points.x,points.y] = projfwd(insert_param.proj,points.lat,points.lon);
  good_mask = ~isnan(points.elev);
  points.x = points.x(good_mask);
  points.y = points.y(good_mask);
  points.elev = points.elev(good_mask);
  
  insert_param.eval.ref_source.name = 'surface';
  insert_param.eval.ref_source.source = 'ops';
  insert_param.eval.ref_gaps_fill.method = 'interp_finite';
  insert_param.eval.cmd = 'source = (elev - ref.twtt*c/2)/(c/2/sqrt(er_ice)) - source + ref.twtt;';
  insert_param.x = points.x;
  insert_param.y = points.y;
  insert_param.data = points.elev/(c/2/sqrt(er_ice));
  
  insert_param.type = 'point'; % Point data
  insert_param.layer_dest.name = 'jakobshavn3D';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'Basic 3D results for 20060530_08'; % For OPS layer_dest source
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'preserve_gaps';
  insert_param.gaps_fill.method_args = [300 60];
  opsInsertLayer(params, insert_param);
  
elseif 0
  %% Example 4: Merge Geoid data onto GIMP layer
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_Polar6.xls'));
  
  % Load in Geoid
  geoid_fn = ct_filename_gis([],'world\egm96_geoid\WW15MGH.DAC');
  points = [];
  [points.lat,points.lon,points.elev] = egm96_loader(geoid_fn);
  points.lon = [points.lon 360];
  points.elev = [points.elev points.elev(:,1)];
  [points.lon,points.lat] = meshgrid(points.lon,points.lat);
  
  % Create geotiff projection structure
  % 1. Grab from a geotiff file with the same projection
  insert_param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));
  
  [points.x,points.y] = projfwd(insert_param.proj,points.lat,points.lon);
  good_mask = ~isnan(points.elev);
  points.x = points.x(good_mask);
  points.y = points.y(good_mask);
  points.elev = points.elev(good_mask);
  
  insert_param.eval.ref_source.name = 'surface';
  insert_param.eval.ref_source.source = 'layerData';
  insert_param.eval.ref_source.layerdata_source = 'layerData';
  insert_param.eval.ref_gaps_fill.method = 'interp_finite';
  insert_param.eval.cmd = 'source = (elev - ref.twtt*c/2)/(c/2/sqrt(er_ice)) - source + ref.twtt;';
  insert_param.x = points.x;
  insert_param.y = points.y;
  insert_param.data = points.elev/(c/2/sqrt(er_ice));
  
  insert_param.type = 'point'; % Point data
  insert_param.layer_dest.name = 'GIMP';
  insert_param.layer_dest.source = 'layerData';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'GIMP Grid'; % For OPS layer_dest source
  insert_param.layer_dest.layerdata_source = 'layerData'; % For layerData layer_dest source
  insert_param.copy_method = 'fillgaps';
  insert_param.gaps_fill.method = 'interp_finite';
  opsInsertLayer(params, insert_param);
  
elseif 0
  %% Example 5: mass conservation grid
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
    
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_Polar6.xls'));
  
  % grid_fn = '/cresis/snfs1/scratch/paden/mass_conservation/HelheimStream-2014-11-18.nc';
  grid_fn = '/cresis/snfs1/scratch/paden/mass_conservation/KangerdlugssuaqStream-2014-11-18.nc';
  
  % Create geotiff projection structure
  % 1. Grab from a geotiff file with the same projection
  % insert_param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));
  % 2. Create an mstruct with the same projection
  %[x,y,mstruct] = geodetic_to_stereographic(70,-45);
  % 3. Hand construction of geotiff projection OR
  insert_param.proj = [];
  insert_param.proj.CornerCoords = [];
  insert_param.proj.Ellipsoid = [];
  insert_param.proj.PM = [];
  insert_param.proj.PMLongToGreenwich = [];
  insert_param.proj.Zone = [];
  geoid = almanac('earth',ncreadatt(grid_fn,'polar_stereographic','ellipsoid'),'m');
  insert_param.proj.SemiMajor = geoid;
  insert_param.proj.SemiMinor = minaxis(geoid);
  insert_param.proj.ProjParm = zeros(7,1);
  insert_param.proj.ProjParm(1) = ncreadatt(grid_fn,'polar_stereographic','standard_parallel');
  insert_param.proj.ProjParm(2) = ncreadatt(grid_fn,'polar_stereographic','straight_vertical_longitude_from_pole');
  insert_param.proj.ProjParm(3) = 0;
  insert_param.proj.ProjParm(4) = 0;
  insert_param.proj.ProjParm(5) = 1;
  insert_param.proj.ProjParm(6) = ncreadatt(grid_fn,'polar_stereographic','false_easting');
  insert_param.proj.ProjParm(7) = ncreadatt(grid_fn,'polar_stereographic','false_northing');
  
  insert_param.proj.CTProjection = 'CT_PolarStereographic'; % Choose from list in toolbox/map/mapproj/private/projcode
  insert_param.proj.GeoTIFFCodes.CTProjection = 15; % CT_PolarStereographic from toolbox/map/mapproj/private/projcode
  
  % Default GeoTIFFCode values
  insert_param.proj.GeoTIFFCodes.Model = int16(1);
  insert_param.proj.GeoTIFFCodes.PCS = int16(32767);
  insert_param.proj.GeoTIFFCodes.GCS = int16(32767);
  insert_param.proj.GeoTIFFCodes.UOMAngle = int16(32767);
  insert_param.proj.GeoTIFFCodes.Datum = int16(32767);
  insert_param.proj.GeoTIFFCodes.PM = int16(32767);
  insert_param.proj.GeoTIFFCodes.ProjCode = int16(32767);
  insert_param.proj.GeoTIFFCodes.Projection= int16(32767);
  insert_param.proj.GeoTIFFCodes.MapSys = int16(32767);
  
  % Default units to meters for now
  insert_param.proj.GeoTIFFCodes.UOMLength = int16(9001);
  
  insert_param.eval.ref_source.name = 'surface';
  insert_param.eval.ref_source.source = 'ops';
  insert_param.eval.ref_gaps_fill.method = 'interp_finite';
  insert_param.eval.cmd = 'source = ref.twtt + source;';
  insert_param.x = ncread(grid_fn,'x');
  insert_param.y = ncread(grid_fn,'y');
  insert_param.data = ncread(grid_fn,'thickness').' / (c/2/sqrt(er_ice));
  
  insert_param.type = 'raster'; % Raster data
  insert_param.layer_dest.name = 'mc_bottom';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'Mass conservation ice bottom'; % For OPS layer_dest source
  insert_param.layer_dest.existence_check = false; % Create layer if it does not exist
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'interp_finite';
  opsInsertLayer(params, insert_param);  
  
elseif 0
  %% Example 6: Gogineni's Jakobshavn grid
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'));
  
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
  insert_param.proj = geotiffinfo(ct_filename_gis([],'greenland\Landsat-7\Greenland_natural_90m.tif'));
    
  insert_param.eval.ref_source.name = 'surface';
  insert_param.eval.ref_source.source = 'ops';
  insert_param.eval.ref_gaps_fill.method = 'interp_finite';
  insert_param.eval.cmd = 'source = ref.twtt + source;';
  insert_param.x = x_axis;
  insert_param.y = y_axis;
  insert_param.data = thickness / (c/2/sqrt(er_ice));
  
  insert_param.type = 'raster'; % Raster data
  insert_param.layer_dest.name = 'gogineni2014';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'Gogineni JofG 2014 grid'; % For OPS layer_dest source
  insert_param.layer_dest.existence_check = false; % Create layer if it does not exist
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'interp_finite';
  opsInsertLayer(params, insert_param); 
  
elseif 0
  %% Example 7: Joe Plummer's Jakobshavn grid
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('rds_param_2009_Greenland_TO.xls'));
  
  grid_fn = ct_filename_gis([],'greenland/plummer_jakobshavn/jak_grid.tif');
  
  insert_param.proj = geotiffinfo(grid_fn);
  
  % Load the grid
  [bt_elev,R] = geotiffread(grid_fn);
  bt_elev(bt_elev == min(min(bt_elev))) = NaN;
  x_axis = R.XLimWorld(1) + [R.XLimIntrinsic(1):R.XLimIntrinsic(2)-1]'*R.DeltaX;
  y_axis = R.YLimWorld(2) + [R.YLimIntrinsic(1):R.YLimIntrinsic(2)-1]'*R.DeltaY;
    
  insert_param.eval.ref_source.name = 'surface';
  insert_param.eval.ref_source.source = 'ops';
  insert_param.eval.ref_gaps_fill.method = 'interp_finite';
  insert_param.eval.cmd = 'source = (elev - ref.twtt*c/2)/(c/2/sqrt(er_ice)) - source + ref.twtt;';
  insert_param.x = x_axis;
  insert_param.y = y_axis;
  insert_param.data = bt_elev/(c/2/sqrt(er_ice));
  
  insert_param.type = 'raster'; % Raster data
  insert_param.layer_dest.name = 'plummer';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'Plummer Jakobshavn Grid'; % For OPS layer_dest source
  insert_param.layer_dest.existence_check = false; % Create layer if it does not exist
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'interp_finite';
  opsInsertLayer(params, insert_param);
  
elseif 0
  %% Example 8: GIMP grid
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  %params = read_param_xls(ct_filename_param('snow_param_2015_Greenland_Polar6.xls'));
  params = read_param_xls(ct_filename_param('rds_param_2016_Greenland_G1XB.xls'),'','post');
  
  grid_fn = ct_filename_gis([],'greenland/DEM/GIMP/gimpdem_90m.tif');
  
  % Load the grid
  points = [];
  [points.elev,R] = geotiffread(grid_fn);
  points.elev = double(points.elev);
  x_axis = R.XLimWorld(1) + [R.XLimIntrinsic(1):R.XLimIntrinsic(2)-1]'*R.DeltaX;
  y_axis = R.YLimWorld(2) + [R.YLimIntrinsic(1):R.YLimIntrinsic(2)-1]'*R.DeltaY;
  
  insert_param.proj = geotiffinfo(grid_fn);
  
  insert_param.eval.cmd = 'source = (elev - source)/(c/2);';
  insert_param.x = x_axis;
  insert_param.y = y_axis;
  insert_param.data = points.elev;
  
  insert_param.type = 'raster'; % Raster data
  insert_param.layer_dest.name = 'GIMP';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'GIMP Grid'; % For OPS layer_dest source
  insert_param.layer_dest.layerdata_source = 'layerData'; % For layerData layer_dest source
  insert_param.layer_dest.existence_check = false; % Create layer if it does not exist
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'interp_finite';
  opsInsertLayer(params, insert_param);
  
elseif 1
  %% Example 9: BEDMAP2 grid
  % =====================================================================
  % =====================================================================
  % Use parameters spreadsheet to select segment and frame list for creating layers
  % Set the generic to 1 for the selected segments and frames
  
  physical_constants;
  insert_param = [];
  
  params = read_param_xls(ct_filename_param('rds_param_2009_Antarctica_TO_ndh_targetframes.xls'),'','post');
  
  grid_fn = ct_filename_gis([],'antarctica/DEM/BEDMAP2/original_data/bedmap2_tiff/bedmap2_surface.tif');
  
  % Load the grid
  points = [];
  [points.elev,R] = geotiffread(grid_fn);
  points.elev = double(points.elev);
  x_axis = R.XLimWorld(1) + [R.XLimIntrinsic(1):R.XLimIntrinsic(2)-1]'*R.DeltaX;
  y_axis = R.YLimWorld(2) + [R.YLimIntrinsic(1):R.YLimIntrinsic(2)-1]'*R.DeltaY;
  
  insert_param.proj = geotiffinfo(grid_fn);
  
  insert_param.eval.cmd = 'source = (elev - source)/(c/2);';
  insert_param.x = x_axis;
  insert_param.y = y_axis;
  insert_param.data = points.elev;
  
  insert_param.type = 'raster'; % Raster data
  insert_param.layer_dest.name = 'BEDMAP_surface';
  insert_param.layer_dest.source = 'ops';
  insert_param.layer_dest.username = 'paden'; % For OPS layer_dest source
  insert_param.layer_dest.group = 'standard'; % For OPS layer_dest source
  insert_param.layer_dest.description = 'BEDMAP 2 Surface'; % For OPS layer_dest source
  insert_param.layer_dest.layerdata_source = 'layerData'; % For layerData layer_dest source
  insert_param.layer_dest.existence_check = false; % Create layer if it does not exist
  insert_param.copy_method = 'overwrite';
  insert_param.gaps_fill.method = 'interp_finite';
  opsInsertLayer(params, insert_param);
   
end
