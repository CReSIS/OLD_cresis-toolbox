% script coverage_maps_with_param_file
%
% Creates data coverage maps (called by coverage_maps.m, which includes
%     user-defined settings.)
% 1. Plots blue line for all GPS trajectory
% 2. Computes coverage rate which is the percentage of time that the radar
%    is on starting from the time it is turned on to the time it is turned
%    off.
% 3. Plots red line on top of this blue line of when radar was on.
% 4. Plots green dots on red line where ice bottom was detected.
%
% Authors: John Paden
% Contributions by Steven Foga

% ==================================================================
% Automated Section
% ==================================================================
for param_idx = 1:length(param)
  clear gps_legend all_legend good_legend mod_legend det_legend;
  fig(param_idx) = figure(101); clf;
  set(fig(param_idx),'Position',fig_size);
    for geotiff_fns_idx = 1:length(geotiff_fns)
    proj = geotiffinfo(geotiff_fns{geotiff_fns_idx});
    [RGB{geotiff_fns_idx}, R{geotiff_fns_idx}, tmp] = geotiffread(geotiff_fns{geotiff_fns_idx});
    R{geotiff_fns_idx} = R{geotiff_fns_idx}/1e3;
    mapshow(RGB{geotiff_fns_idx}, R{geotiff_fns_idx});
    hold on;
  end
  hold on;
  axis([R{1}(3,1)+[0 R{1}(2,1)*(size(RGB{1},2)-1)] sort(R{1}(3,2)+[0 R{1}(1,2)*(size(RGB{1},1)-1)])]);
  
  post_dir = ct_filename_out(param(param_idx),'','CSARP_post',1);
  csv_dir = fullfile(post_dir,'csv');
  
  coverage_image_fn = fullfile(out_dir,sprintf('coverage_%s_%s.jpg',param(param_idx).radar_name,param.season_name));
  
  for day_seg_idx = 1:length(param.segs)
    days{day_seg_idx} = param(param_idx).segs{day_seg_idx}(1:8);
  end
  days = unique(days);
  
  clear good_pnts; clear all_pnts; clear all_pnts_sec; clear all_gps_sec;
  good_pnts = 0;
  all_pnts = 0;
  all_pnts_sec = 0;
  all_gps_sec = 0;
  
  for day_idx = 1:length(days)
    day = days{day_idx};
    csv_fns = get_filenames(csv_dir, sprintf('Data_%s',day),'','.csv');
    
    param(param_idx).day_seg = [day '_01'];
    gps_path = (ct_filename_support(param(param_idx),'','gps',true));
    if isempty(dir(gps_path))
      continue
    else
      gps = load(gps_path);
    end
    [gps.gps_time good_idxs] = unique(gps.gps_time);
    gps.lat = gps.lat(good_idxs);
    gps.lon = gps.lon(good_idxs);
    gps.elev = gps.elev(good_idxs);
    % Decimation
    along_track = geodetic_to_along_track(gps.lat,gps.lon,gps.elev);
    decim_idxs = get_equal_alongtrack_spacing_idxs(along_track,50);
    gps.lat = gps.lat(decim_idxs);
    gps.lon = gps.lon(decim_idxs);
    gps.elev = gps.elev(decim_idxs);
    gps.gps_time = gps.gps_time(decim_idxs);
    gps.mask = zeros(size(gps.lat));
    [gps.x,gps.y] = projfwd(proj,gps.lat,gps.lon);
    gps.x = gps.x / 1000;
    gps.y = gps.y / 1000;
    figure(fig(param_idx));
    h_gps = plot(gps.x,gps.y,'-b');
    
    for day_seg_idx = 1:length(csv_fns)
      csv_fn = csv_fns{day_seg_idx};
      fid = fopen(csv_fn);
      C = textscan(fid,'%f%f%f%f%f%s%f%f%f','Delimiter',',','Headerlines',1);
      fclose(fid);
      clear lat; clear lon; clear bottom elev; clear utc_time_sod; clear frame_id;
      clear year; clear month; clear day; clear utc_time; clear gps_time;
      clear x_good; clear y_good; clear qual moderate_qual_idxs detrend_qual_idxs;
      clear along_track_csv decim_csv_idxs;
      lat = C{1};
      lon = C{2};
      % Decimation
      along_track_csv = geodetic_to_along_track(lat,lon,(zeros(size(lat))));
      decim_csv_idxs = get_equal_alongtrack_spacing_idxs(along_track_csv,50);
      lat = lat(decim_csv_idxs);
      lon = lon(decim_csv_idxs);
      bottom = C{8};
      bottom = bottom(decim_csv_idxs);
      qual = C{9};
      qual = qual(decim_csv_idxs);
      utc_time_sod = C{3};
      utc_time_sod = utc_time_sod(decim_csv_idxs);
      frame_id = C{6};
      frame_id = frame_id(decim_csv_idxs);
      year = str2double(frame_id{1}(1:4));
      month = str2double(frame_id{1}(5:6));
      day = str2double(frame_id{1}(7:8));
      utc_time = datenum_to_epoch(datenum(year,month,day,0,0,utc_time_sod));
      gps_time = utc_time + utc_leap_seconds(utc_time(1));
      [x,y] = projfwd(proj,lat,lon);
      x = x / 1e3;
      y = y / 1e3;
      good_idxs = find(bottom ~= -9999);
      % For each segment, set the mask to 1 during this segment to indicate
      % that the radar was turned on (i.e. "coverage" variable below)
      gps.mask(gps_time(1) <= gps.gps_time & gps.gps_time <= gps_time(end)) = 1;
      good_pnts = good_pnts + length(good_idxs);
      all_pnts = all_pnts + length(bottom);
      moderate_qual_idxs = find(qual == 2);
      detrend_qual_idxs = find(qual == 3);
      [x_good,y_good] = projfwd(proj,lat(good_idxs),lon(good_idxs));
      x_good = x_good/1e3;
      y_good = y_good/1e3;

      [X_mod,Y_mod] = projfwd(proj,lat(moderate_qual_idxs),lon(moderate_qual_idxs));
      X_mod = X_mod/1e3;
      Y_mod = Y_mod/1e3;

      [X_det,Y_det] = projfwd(proj,lat(detrend_qual_idxs),lon(detrend_qual_idxs));
      X_det = X_det/1e3;
      Y_det = Y_det/1e3;
      
      first_idx = find(gps.mask == 1, 1);
      last_idx = find(gps.mask == 1, 1, 'last');
      
      % Interpolate to 1 second intervals from first measurement of the
      % day to the last measurement of the day
      gps.mask_interp = interp1(gps.gps_time, gps.mask, gps.gps_time(first_idx):gps.gps_time(last_idx));
      all_pnts_sec = all_pnts_sec + sum(gps.mask_interp);
      all_gps_sec = all_gps_sec + length(gps.mask_interp);

      figure(fig(param_idx));
      h_all = plot(x,y,'-r');

      figure(fig(param_idx));
      h_good = plot(x_good,y_good,'.g');

      figure(fig(param_idx));
      h_mod = plot(X_mod,Y_mod,'.y');

      figure(fig(param_idx));
      h_det = plot(X_det,Y_det,'.r');
    end
  end
  
  clear coverage; clear bottom_good;
  coverage = all_pnts_sec / all_gps_sec * 100;
  
  bottom_good = good_pnts / all_pnts * 100;
  % =========================================================================
  % =========================================================================
  if ~exist(out_dir,'dir')
    mkdir(out_dir);
  end
  
  stat_fn_name = sprintf('%s_stats.txt', param.season_name);
  stat_fn_name(stat_fn_name==' ') = '_';
  stat_fn = fullfile(out_dir,stat_fn_name);
  fid = fopen(stat_fn,'w');
  fprintf(fid,'%s\n',param.season_name);
  stat_txt = sprintf('  %.2f%% of data has good ice bottom.\n',bottom_good);
  fprintf(fid,'%s',stat_txt);
  fprintf('%s',stat_txt);
  stat_txt = sprintf('  %.2f%% of the time the radar was turned on.\n',coverage);
  fprintf(fid,'%s',stat_txt);
  fprintf('%s',stat_txt);
  fclose(fid);
  
% Add legend and  plot labels
  figure(fig);
  h_gps = plot(1,NaN,'-b');
  h_all = plot(1,NaN,'-r');
  h_good = plot(1,NaN,'g.');
  h_mod = plot(1,NaN,'y.');
  h_det = plot(1,NaN,'r.');
  h_legend = legend([h_gps h_all h_good h_mod h_det],'Entire Flight','All Radar Data', ...
        'Good Radar Data','Moderate Radar Data','Detrended Radar Data','Location','Best');
  if ~ispc % Linux fonts need to be enlarged for the legend
    set(h_legend,'FontSize',18);
  else
    set(h_legend,'FontSize',6);
  end
    
  season_name = param.season_name;
  season_name = strrep(season_name, '_',' ');
  title(season_name);
  xlabel('X (km)');
  ylabel('Y (km)');
  hold off;
  
  fprintf('%s', datestr(now,'HH:MM:SS\n'));
  
  if ~exist(out_dir,'dir')
    mkdir(out_dir);
  end
  
  if strcmpi(location,'Canada')
    image_fn_name = sprintf('%s_%s_coverage_map.jpg','Canada',param.season_name);
  else % Else Antarctica or Greenland
    image_fn_name = sprintf('%s_coverage_map.jpg', param.season_name);
  end
  
  image_fn = fullfile(out_dir,image_fn_name);
  print(fig(param_idx),'-djpeg','-r300',image_fn);
end
return;
