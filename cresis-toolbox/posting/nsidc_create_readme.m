% script nsidc_create_readme
%
% Creates outputs that are needed in the readme file, including:
% 1. Data coverage maps
% 2. Equations for beamwidth and loop sensitivity
%
% Author: John Paden
% 
% See also: type "nsidc_help.m"

proj = plot_geotiff(fullfile(ct_filename_gis([],'greenland'),'Landsat-7','Greenland_natural.tif'));
fns = get_filenames('Y:/mdce/csarp_support/gps/2011_Greenland_P3','gps','','.mat');
fns2 = get_filenames('Y:/mdce/csarp_support/gps/2011_Greenland_TO','gps','','.mat');
fns = cat(1,fns,fns2);
for fn_idx = 1:length(fns)
    fn = fns{fn_idx};
    gps = load(fn);
    [x,y] = projfwd(proj,gps.lat,gps.lon);
    figure(1);
    hold on;
    plot(x/1e3,y/1e3);
    hold off;
end

% ==================================================================
% User Settings
% ==================================================================

param.radar_name = 'mcords2';
param.season_name = '2011_Greenland_P3';
param.location = 'greenland'; % greenland or antarctica or canada
out_dir = '~/tmp/';
segment_list_mcords2_2011_greenland_P3;

post_dir = ct_filename_out(param,'','CSARP_post',1);
csv_dir = fullfile(post_dir,'csv');

coverage_image_fn = fullfile(out_dir,sprintf('coverage_%s_%s.jpg',param.radar_name,param.season_name));

for day_seg_idx = 1:length(param.segs)
  days{day_seg_idx} = param.segs{day_seg_idx}(1:8);
end
days = unique(days);

% ==================================================================
% Automated Section
% ==================================================================

if strcmpi(param.location,'antarctica')
  geotiff_fn = fullfile(ct_filename_gis(param,'antarctica'),'Landsat-7','Antarctica_LIMA.tif');
elseif strcmpi(param.location,'greenland')
  geotiff_fn = fullfile(ct_filename_gis(param,'greenland'),'Landsat-7','Greenland_natural.tif');
elseif strcmpi(param.location,'canada')
  geotiff_fn = fullfile(ct_filename_gis(param,'canada'),'Landsat-7','Canada.tif');
else
  error('Invalid param.location');
end

proj = plot_geotiff(geotiff_fn, [], [], 1);
xlabel('X (km)');
ylabel('Y (km)');

good_pnts = 0;
all_pnts = 0;
all_pnts_sec = 0;
all_gps_sec = 0;

for day_idx = 1:length(days)
  day = days{day_idx};
  csv_fns = get_filenames(csv_dir, sprintf('Data_%s',day),'','.csv');
  
  param.day_seg = [day '_01'];
  gps = load(ct_filename_support(param,'','gps',true));
  [gps.gps_time good_idxs] = unique(gps.gps_time);
  gps.lat = gps.lat(good_idxs);
  gps.lon = gps.lon(good_idxs);
  gps.elev = gps.elev(good_idxs);
  gps.mask = zeros(size(gps.lat));
  [gps.x,gps.y] = projfwd(proj,gps.lat,gps.lon);
  gps.x = gps.x / 1000;
  gps.y = gps.y / 1000;
  
  hold on;
  plot(gps.x,gps.y,'b');
  hold off;
  
  for day_seg_idx = 1:length(csv_fns)
    csv_fn = csv_fns{day_seg_idx}
    fid = fopen(csv_fn);
    C = textscan(fid,'%f%f%f%f%f%s%f%f%f','Delimiter',',','Headerlines',1);
    fclose(fid);
    lat = C{1};
    lon = C{2};
    bottom = C{8};
    utc_time_sod = C{3};
    frame_id = C{6};
    year = str2double(frame_id{1}(1:4));
    month = str2double(frame_id{1}(5:6));
    day = str2double(frame_id{1}(7:8));
    utc_time = datenum_to_epoch(datenum(year,month,day,0,0,utc_time_sod));
    gps_time = utc_time + utc_leap_seconds(utc_time(1));
    [x,y] = projfwd(proj,lat,lon);
    x = x / 1e3;
    y = y / 1e3;
    good_idxs = find(bottom ~= -9999);
    gps.mask(gps_time(1) <= gps.gps_time & gps.gps_time <= gps_time(end)) = 1;
    good_pnts = good_pnts + length(good_idxs);
    all_pnts = all_pnts + length(bottom);
    hold on;
    plot(x,y,'r');
    plot(x(good_idxs),y(good_idxs),'g.');
    hold off;
  end
  
  first_idx = find(gps.mask == 1, 1);
  last_idx = find(gps.mask == 1, 1, 'last');
  
  % Interpolate to 1 second intervals from first measurement of the
  % day to the last measurement of the day
  gps.mask_interp = interp1(gps.gps_time, gps.mask, gps.gps_time(first_idx):gps.gps_time(last_idx));
  
  all_pnts_sec = all_pnts_sec + sum(gps.mask_interp);
  all_gps_sec = all_gps_sec + length(gps.mask_interp);
  
end

coverage = all_pnts_sec / all_gps_sec * 100

bottom_good = good_pnts / all_pnts * 100

saveas(1,coverage_image_fn);

kt = 1.6;
physical_constants;
BW = 30e6;
kt*c/(2*BW*sqrt(er_ice))

H = 500;
T = 2000;
2*sqrt( (H + T/sqrt(er_ice)) * c * kt / BW )

% Number of tx elements
Ny = 7;
% Element spacing
dy = 0.5;
% Beamwidth
beta_y = 1/((Ny+1)*dy) * 180/pi

% Cross-track antenna array window
ky = 1.6;
2*(H + T/sqrt(er_ice)) * tan(beta_y/180*pi * ky / 2)

% Loop sensitivity
Pt = 75*Ny;
Nc = 7;
G = 4;
coh_ave = 32*100; % including SAR aperture
fc = 195e6;
lambda_fc = c/fc;
Tpd = 10e-6;
Fnoise = 2;

10*log10(Pt*coh_ave*BW*Tpd*(Ny*G*lambda_fc)^2/(4*pi) / (BoltzmannConst*290*BW*Fnoise) )







return;
