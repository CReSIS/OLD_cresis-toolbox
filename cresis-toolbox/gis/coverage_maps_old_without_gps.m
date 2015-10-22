% script coverage_maps_old_without_gps.m
%
% Creates data coverage maps for ever year before 2001 (since these years
% do not have dedicated GPS files.) They are: 1993, 1995, 1996, 1997, 1998,
% 1999, 2001 and 2002_Greenland_P3.
% 1. Plots blue line for all GPS trajectory
% 2. Computes coverage rate which is the percentage of time that the radar
%    is on starting from the time it is turned on to the time it is turned
%    off.
% 3. Plots red line on top of this blue line of when radar was on.
% 4. Plots green dots on red line where ice bottom was detected.

% For ever year before 1993 (since these years do not have dedicated GPS
% files.) They are: 1993, 1993, 1993, 1993, 1993, 1993, 1993,
% and 1993_Greenland_P3.

% Author: Steven Foga

%----------------------------------------
% Script
%----------------------------------------

for param_idx = 1:length(param)
  
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
  coverage_image_fn = fullfile(out_dir,sprintf('coverage_%s_%s.jpg',param.season_name));
  file_end = '.csv';
  param(1,length(param)).recursive = 1;
  flights = param.all_season_type;
  
  for file_idx = 1:length(flights)
    data_all = fopen(flights,'r');
    C = textscan(data_all, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
    clear lat; clear lon; clear bottom; clear utc_time_sod; clear frame_id;
    clear year; clear month; clear day; clear utc_time; clear gps_time;
    clear x_good; clear y_good;
    lat = C{1};
    lon = C{2};
    utc_time_sod = C{3};
    fclose(data_all);
    if utc_time_sod > 0
      gpsDiff = diff(utc_time_sod);
      gpsOff = gpsDiff >= 30; % If the gap is >= 30 seconds
      gpsOff2 = find(gpsOff==1);
      totalOff(file_idx) = (sum(gpsDiff(gpsOff2)));
      totalOn(file_idx) = max(utc_time_sod) - min(utc_time_sod);
      
    else
      totalOff(file_idx) = 0;
      totalOn(file_idx) = 0;
    end

    gpsOffTotal = sum(totalOff);
    coverage_stat = (sum(totalOn) / (gpsOffTotal + sum(totalOn)))*100;
    [x,y] = projfwd(proj,lat,lon);
    x = x / 1e3;
    y = y / 1e3;
    figure(fig(param_idx));
    h_all = plot(x,y,'-r');
  end
  
  % Good data pct
  clear all_data_fn data_all_csv
  all_data_fn = param(param_idx).all_season_type;
  data_all_csv = fopen(all_data_fn,'r');
  F = textscan(data_all_csv, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
  clear lat_all lon_all bottom_all good_idxs totalGood totalAll bottom_good_stat x_good y_good
  lat_all = F{1};
  lon_all = F{2};
  bottom_all = F{8};
  qual_all = F{9};
  fclose(data_all_csv);
  good_idxs = find(bottom_all ~= -9999);
  totalGood = length(good_idxs);
  totalAll = length(bottom_all);
  bottom_good_stat = ((totalGood/totalAll)*100);
  moderate_qual_idxs = find(qual_all == 2);
  detrend_qual_idxs = find(qual_all == 3);
  [x_good,y_good] = projfwd(proj,lat_all(good_idxs),lon_all(good_idxs));
  x_good = x_good/1e3;
  y_good = y_good/1e3;
  
  [X_mod,Y_mod] = projfwd(proj,lat_all(moderate_qual_idxs),lon_all(moderate_qual_idxs));
  X_mod = X_mod/1e3;
  Y_mod = Y_mod/1e3;
  
  [X_det,Y_det] = projfwd(proj,lat_all(detrend_qual_idxs),lon_all(detrend_qual_idxs));
  X_det = X_det/1e3;
  Y_det = Y_det/1e3;
  
  disp(bottom_good_stat); disp(coverage_stat);
  clear h_good
  figure(fig(param_idx));
  h_good = plot(x_good,y_good,'.g');

  figure(fig(param_idx));
  h_mod = plot(X_mod,Y_mod,'.y');

  figure(fig(param_idx));
  h_det = plot(X_det,Y_det,'.r');
  
  %--------------
  % Plot
  %--------------
  if ~exist(out_dir,'dir')
    mkdir(out_dir);
  end
  
  stat_fn_name = sprintf('%s_stats.txt', param.season_name);
  stat_fn_name(stat_fn_name==' ') = '_';
  stat_fn = fullfile(out_dir,stat_fn_name);
  fid = fopen(stat_fn,'w');
  fprintf(fid,'%s\n',param.season_name);
  stat_txt = sprintf('  %.2f%% of data has good ice bottom.\n',bottom_good_stat);
  fprintf(fid,'%s',stat_txt);
  fprintf('%s',stat_txt);
  stat_txt = sprintf('  %.2f%% of the time the radar was turned on.\n',coverage_stat);
  fprintf(fid,'%s',stat_txt);
  fprintf('%s',stat_txt);
  fclose(fid);

% Add single-season plot labels
figure(fig(param_idx));
h_all = plot(1,NaN,'-r');
h_good = plot(1,NaN,'g.');
h_mod = plot(1,NaN,'y.');
h_det = plot(1,NaN,'r.');
h_legend = legend([h_all h_good h_mod h_det],'All Radar Data','Good Radar Data',...
        'Moderate Radar Data','Detrended Radar Data','Location','Best');
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



