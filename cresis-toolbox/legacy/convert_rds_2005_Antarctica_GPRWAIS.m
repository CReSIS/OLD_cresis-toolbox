% Script convert_2005_Antarctica_GPRWAIS.m
%
% Author: John Paden

%% ======================================================================
% User Settings
% =======================================================================
clear;

base_dir = '/cresis/scratch2/mdce/sar/2005_Antarctica_GPRWAIS/MAT';

base_out_dir = '/cresis/scratch2/mdce/sar/2005_Antarctica_GPRWAIS/CSARP_standard/';

% INCOMPLETE: gps_dir = '/cresis/data1/SAR/2006_Antarctica/insar_gps/nmeafiles/';
% 20060103, 20060103b have no overlap
% 20060104.mat, 20060105.mat, 20060106.mat, 20060107.mat some overlap

gps_dir = '/cresis/data1/SAR/2006_Antarctica/insar2006/';

% William Blake confirmed that 3.15 and 3e8 were used.
%   169 m/us was confirmed in Laird publication (which is 3.1511...)
er_ice = 3.15;
c = 3e8;

% Found by comparing min/max ice thicknesses reported in publication with
% those in the files
system_delay_correction_in_meters = 34;

create_files = false;
post_processing = true;
plot_flag = false;

wais.lat = -79.467;
wais.lon = -112.085;
% Bottom of borehole location from Jeff:
wais.lat = -79.4689;
wais.lon = -112.0833;

%% ======================================================================
% Automated Section
% =======================================================================

[depth,er] = waisPerm(150e6,10001,0);
[TWtime,gain] = genPropProfileFromPerm(depth,er,150e6);

if create_files && exist(base_out_dir,'dir')
  rmdir(base_out_dir,'s');
end

fns = get_filenames(base_dir,'','','.mat');

figure(1); clf; hold on;

colors = {'k','r','g','c','b','m'};

minval = [];
maxval = [];
lat_all = [];
lon_all = [];
elev_all = [];
bottom_all = [];
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  [fn_path fn_name] = fileparts(fn);
  fprintf('Converting %s\n', fn_name);
  year = str2double(fn_name(1:4));
  month = str2double(fn_name(5:6));
  day = str2double(fn_name(7:8));
  fprintf('  %.0f %.0f %.0f\n', year, month, day);
  day_seg = [fn_name(1:8) '_01'];
  
  data = load(fn);
  data.Depths = data.Depths - system_delay_correction_in_meters;
  data.Bedrock_Estimate = data.Bedrock_Estimate - system_delay_correction_in_meters;
  
  Latitude = data.Latitude;
  Longitude = data.Longitude;
  Elevation = data.Elevation;

  gps_dir_day= fullfile(gps_dir, fn_name(3:end), 'fk_data');
  gps_fns = get_filenames(gps_dir_day,'','','fk_pos.mat');
  UTC_time = [];
  lat = [];
  lon = [];
  for gps_fns_idx = 1:length(gps_fns)
    gps_fn_name = gps_fns{gps_fns_idx};
    fprintf('    %s\n', gps_fn_name);
    gps = load(gps_fn_name);
    UTC_time = [UTC_time gps.UTC_time];
    lat = [lat gps.lat];
    lon = [lon gps.lon];
  end

  time_idx = zeros(size(Latitude));
  for idx = 1:length(Latitude)
    [val time_idx(idx)] = min(abs(Latitude(idx)-lat) + abs(Longitude(idx)-lon));
  end
  if plot_flag
    figure(2); clf;
    plot(lat)
    hold on;
    plot(Latitude,'rx')
    figure(3); clf;
    plot(lon)
    hold on;
    plot(Longitude,'rx')
    figure(4); clf;
    plot(diff(time_idx));
  end
  
  GPS_time = UTC_time(time_idx) + utc_leap_seconds(UTC_time(1));
  for idx = 2:length(GPS_time)
    if GPS_time(idx) <= GPS_time(idx-1)
      GPS_time(idx) = GPS_time(idx-1) + 1e-6;
    end
  end
  Surface = zeros(size(Latitude));
  Bottom = data.Bedrock_Estimate;
  Data = data.Data;
  
  Depth = data.Depths;
  Time = Depth / (c/2/sqrt(er_ice));
  fprintf('  dt = %f ns, BW = %f Hz\n', mean(diff(Time))*1e9, 1/mean(diff(Time)));
  Bottom = interp1(Depth,Time,Bottom);
  if plot_flag
    figure(1);
    col_idx = mod(fn_idx-1,length(colors)) + 1;
    plot(lon,lat,[colors{col_idx} '.']);
    plot(Longitude,Latitude,colors{col_idx});
    plot(Longitude(1),Latitude(1),[colors{col_idx}, 'o']);
    plot(Longitude(end),Latitude(end),[colors{col_idx}, 'x']);
  end
  
  Latitude = fir_dec(Latitude,10);
  Longitude = fir_dec(Longitude,10);
  Elevation = fir_dec(Elevation,10);
  GPS_time = fir_dec(GPS_time,10);
  Surface = fir_dec(Surface,10);
  Bottom = fir_dec(Bottom,10);
  Data = sgolayfilt(Data,3,25,[],2);
  Data = Data(:, 5 + 10*(0:length(Latitude)-1) );
  
  minval = [minval min(data.Bedrock_Estimate)];
  maxval = [maxval max(data.Bedrock_Estimate)];

  if create_files
    out_fn = fullfile(base_out_dir,day_seg,sprintf('Data_%s_001.mat',day_seg));
    out_dir = fileparts(out_fn);
    if ~exist(out_dir,'dir')
      mkdir(out_dir);
    end
    frm = 1;
    while exist(out_fn,'file')
      frm = frm + 1;
      out_fn = fullfile(base_out_dir,day_seg,sprintf('Data_%s_%03d.mat',day_seg,frm));
    end
    fprintf('  Output file %s\n', out_fn);
    save(out_fn,'-v6','Latitude','Longitude','Elevation','GPS_time','Surface','Bottom','Data','Time','Depth');
  end

  along_track = geodetic_to_along_track(Latitude,Longitude,Elevation);
  lat_all = [lat_all Latitude];
  lon_all = [lon_all Longitude];
  elev_all = [elev_all Elevation];
  bottom_all = [bottom_all interp1(TWtime,depth(2:end),Bottom)];
%   bottom_all = [bottom_all Bottom*c/2/sqrt(er_ice)];
end

if post_processing
  round([min(minval) max(maxval)])
  
  % 3460 m, 3467 m with 3.15
  
  global gRadar;
  proj = geotiffinfo(fullfile(ct_filename_gis(gRadar,'antarctica'),'Landsat-7','Antarctica_LIMA_480m.tif'));
  [x,y] = projfwd(proj,lat_all,lon_all);
  [wais.x,wais.y] = projfwd(proj,wais.lat,wais.lon);
  xi = linspace(min(x),max(x),1000).';
  yi = linspace(min(y),max(y),200);
  if 0
    bottom = griddata(x,y,bottom_all,xi,yi);
    wais.bottom = griddata(x,y,bottom_all,wais.x,wais.y)
  else
    dt = DelaunayTri(x.',y.');
    F = TriScatteredInterp(dt,bottom_all.','natural');
    [xi_mesh,yi_mesh] = meshgrid(xi,yi);
    bottom = F(xi_mesh,yi_mesh);
    wais.bottom = F(wais.x,wais.y)
  end
  figure(1); clf;
  imagesc(xi/1e3,yi/1e3,bottom);
  set(gca,'YDir','normal');
  contour(xi/1e3,yi/1e3,bottom,[3400:5:3500]);
  % colormap([1 1 1; flipud(jet(255))]);
  orig_caxis = caxis;
  caxis([orig_caxis(1)-4 orig_caxis(2)]);
  caxis([3400 3500]);
  colorbar;
  hold on;
  plot(wais.x/1e3,wais.y/1e3,'kx')
  plot(wais.x/1e3,wais.y/1e3,'ko')
  plot(x/1e3,y/1e3,'k.');
  hold off;
  axis equal;
  xlabel('X (km)');
  ylabel('Y (km)');
  
  DEM = bottom;
  eastAxis = xi;
  northAxis = yi;
  
  R = spatialref.MapRasterReference;
  R.XLimWorld = [eastAxis(1) eastAxis(end)];
  R.YLimWorld = [northAxis(1) northAxis(end)];
  R.RasterSize = size(DEM);
  R.RasterInterpretation = 'cells';
  R.ColumnsStartFrom = 'south';
  R.RowsStartFrom = 'west';
  R.RasterInterpretation = 'postings';
  
  key.GTModelTypeGeoKey  = 1;  % Projected Coordinate System (PCS)
  key.GTRasterTypeGeoKey = 2;  % PixelIsPoint
  %   key.ProjectedCSTypeGeoKey = 32622; % Russell
  key.ProjectedCSTypeGeoKey = 3031; % Antarctica LIMA
  filename = '/users/paden/WAIS_geotiff.tif';
  DEM(isnan(DEM)) = -9999;
  geotiffwrite(filename, int16(DEM), R, 'GeoKeyDirectoryTag', key);
end

return;

%                  Latitude: [1x3637 double]
%                 Longitude: [1x3637 double]
%                 Elevation: [1x3637 double]
%                  GPS_time: [1x3637 double]
%                   Surface: [1x3637 double]
%                    Bottom: [1x3637 double]
%                      Data: [1839x3637 double]
%                      Time: [1x1839 double]
%                     Depth: [1x1839 double]

