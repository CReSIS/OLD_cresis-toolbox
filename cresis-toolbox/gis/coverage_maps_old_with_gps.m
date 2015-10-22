% script coverage_maps_old_with_gps
%
% This script is called from coverage_maps
%
% It plots out the coverage maps and creates statistics files for
% all the old seasons that satisfy these confitions:
%   1. do not have parameter file
%   2. do have GPS files
%
% See also: coverage_maps, coverage_maps_old_without_gps,
%   coverage_maps_with_param_file
%
% Author: Steven Foga, John Paden

fprintf('===========================================================\n');
fprintf('coverage_maps_old_with_gps\n\n');

% ------------------------------------------------------------------------
% Plot coverage maps and create statistics files
% ------------------------------------------------------------------------

% Calculate statistics - Greenland
totalGood = {}; totalAll = {}; pctGood = {};
clear gps_lat gps_lon;

for season_idx = 1:length(param.stat_good)
  segs = make_segment_list(param.rdr_season_type);
  for seg_idx = 1:length(segs)
    segs{seg_idx} = segs{seg_idx}(1:8);
  end
  radar_days = unique(segs);
  flights = param.all_season_type;
  clear totalOff totalOn;
  for file_idx = 1
    segment_csv_fn = flights;
    season_segments = fopen(segment_csv_fn,'r');
    % Read in CSV file (using field 3, GPSTIMESOD)
    gps_time_sod = textscan(season_segments, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
    gps_time_sod = gps_time_sod{3};
    if gps_time_sod > 0
      gpsDiff = diff(gps_time_sod);
      gpsOff = gpsDiff >= 30; % If the gap is >= 30 seconds
      % Sum up all the time gaps greater than 30 seconds
      totalOff(file_idx) = (sum(gpsDiff(gpsOff)));
      % Get the total time on
      totalOn(file_idx) = max(gps_time_sod) - min(gps_time_sod);
    else
      totalOff(file_idx) = 0;
      totalOn(file_idx) = 0;
      fclose(season_segments);
    end
  end
  gpsOffTotal = sum(totalOff);
  pctOn{season_idx} = (sum(totalOn) / (gpsOffTotal + sum(totalOn)))*100;
  
  gre_good_fn = param.stat_good;
  gre_all_fn = param.stat_all;
  
  fid = fopen(gre_good_fn,'r');
  % Scan csv_good file (only good bottom data)
  tmp = textscan(fid, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
  fclose(fid);
  LAT_good{season_idx} = tmp{1};
  LON_good{season_idx} = tmp{2};
  BOTTOM_good{season_idx} = tmp{8};
  QUAL_good{season_idx} = tmp{9};
  
  fid = fopen(gre_all_fn,'r');
  % Scan csv file (all data)
  tmp = textscan(fid, '%f%f%f%f%f%s%f%f%f','headerlines',1,'delimiter',',');
  LAT{season_idx} = tmp{1};
  LON{season_idx} = tmp{2};
  TIME{season_idx} = tmp{3};
  FRAME{season_idx} = tmp{6};
  BOTTOM{season_idx} = tmp{8};
  fclose(fid);
  
  clear tmp;
  clear moderate_qual_idxs detrend_qual_idxs; 
  moderate_qual_idxs{season_idx} = find(QUAL_good{season_idx} == 2);
  detrend_qual_idxs{season_idx} = find(QUAL_good{season_idx} == 3);
  
  totalGood{season_idx} = length(BOTTOM_good{season_idx});
  totalAll{season_idx} = length(BOTTOM{season_idx});
  pctGood{season_idx} = ((totalGood{season_idx}/totalAll{season_idx})*100);
  
  % --------------------------------------------------------------------
  % Get the lat/lon for the entire trajectory for every day the radar
  % was turned on
  % --------------------------------------------------------------------
  gps_path = param.gps_all;
  gps_end_path = '.mat';
  % Load all the GPS files in a loop
  gps_set{season_idx} = get_filenames(gps_path,'','',gps_end_path,'recursive');
  for gps_idx = 1:length(gps_set{season_idx})
    gps_all{season_idx}{gps_idx} = load(gps_set{season_idx}{gps_idx});
    if isempty(strmatch(datestr(epoch_to_datenum(gps_all{season_idx}{gps_idx}.UTC_time(1)),'yyyymmdd'), radar_days))
      % If this is not a day that the radar flew, then skip it
      continue
    else
      gps_lat{season_idx}{gps_idx} = gps_all{season_idx}{gps_idx}.lat;
      gps_lon{season_idx}{gps_idx} = gps_all{season_idx}{gps_idx}.lon;
    end
  end
end

%% Plot Coverage Maps
fprintf('Statistics calculated (Greenland).\n');

% Download and load the geotiff based on the location selection.
fprintf('-----------------------------------------------------------------\n');
fprintf('Loading and Plotting GeoTIFF. May take a few minutes ... ');
fprintf('%s', datestr(now,'HH:MM:SS\n'));
fprintf('Loading GeoTIFF ... ');

% if run_antarctica
%   clear fig fig_paper_position
% else
% end
for plot_idx = 1
% for plot_idx = 1:length(param.stat_good);
  fig(plot_idx) = figure(101); clf;
  set(fig(plot_idx),'Position',fig_size);
  for geotiff_fns_idx = 1
%   for geotiff_fns_idx = 1:length(geotiff_fns)
    proj = geotiffinfo(geotiff_fns{geotiff_fns_idx});
    [RGB{geotiff_fns_idx}, R{geotiff_fns_idx}, tmp] = geotiffread(geotiff_fns{geotiff_fns_idx});
    R{geotiff_fns_idx} = R{geotiff_fns_idx}/1e3;
    mapshow(RGB{geotiff_fns_idx}, R{geotiff_fns_idx});
    hold on;
  end
  hold on;
  for geotiff_fns_idx = 1
%   for geotiff_fns_idx = 1:length(geotiff_fns)
    mapshow(RGB{geotiff_fns_idx}, R{geotiff_fns_idx});
  end
  axis([R{1}(3,1)+[0 R{1}(2,1)*(size(RGB{1},2)-1)] sort(R{1}(3,2)+[0 R{1}(1,2)*(size(RGB{1},1)-1)])]);
  
  fprintf('%s', datestr(now,'HH:MM:SS\n'));
  fprintf('Loading Flightlines ...\n');
  
  % Plot flightlines
  % Display csv data (all data)
  [X,Y] = projfwd(proj,LAT{plot_idx},LON{plot_idx});
  X = X/1e3;
  Y = Y/1e3;
  % Display csv_good data (bottom only data)
  [X_good,Y_good] = projfwd(proj,LAT_good{plot_idx},LON_good{plot_idx});
  X_good = X_good/1e3;
  Y_good = Y_good/1e3;
  
  clear mod_qual det_qual;
  mod_qual = moderate_qual_idxs{plot_idx};
  det_qual = detrend_qual_idxs{plot_idx};
  
  [X_mod,Y_mod] = projfwd(proj,LAT_good{plot_idx}(mod_qual),LON_good{plot_idx}(mod_qual));
  X_mod = X_mod/1e3;
  Y_mod = Y_mod/1e3;
  
  [X_det,Y_det] = projfwd(proj,LAT_good{plot_idx}(det_qual),LON_good{plot_idx}(det_qual));
  X_det = X_det/1e3;
  Y_det = Y_det/1e3;
  
  if ~exist(out_dir,'dir')
    mkdir(out_dir);
  end
  stat_fn_name = sprintf('%s_stats.txt', param.season_name);
  stat_fn_name(stat_fn_name==' ') = '_';
  stat_fn = fullfile(out_dir,stat_fn_name);
  fid = fopen(stat_fn,'w');
  fprintf(fid,'%s\n',param.season_name);
  stat_txt = sprintf('  %.2f%% of data has good ice bottom.\n',pctGood{plot_idx});
  fprintf(fid,'%s',stat_txt);
  fprintf('%s',stat_txt);
  stat_txt = sprintf('  %.2f%% of the time the radar was turned on.\n',pctOn{plot_idx});
  fprintf(fid,'%s',stat_txt);
  fprintf('%s',stat_txt);
  fclose(fid);
  
  % Display entire flight data
  clear X_gps; clear Y_gps;
  for gps_idx = 1:length(gps_lat{plot_idx})
    [X_gps{gps_idx},Y_gps{gps_idx}] = projfwd(proj,gps_lat{plot_idx}{gps_idx},gps_lon{plot_idx}{gps_idx});
    [len wid] = size(X_gps{gps_idx});
    if len > 1
      X_gps{gps_idx} = (X_gps{gps_idx}/1e3)';
      Y_gps{gps_idx} = (Y_gps{gps_idx}/1e3)';
    else
      X_gps{gps_idx} = X_gps{gps_idx}/1e3;
      Y_gps{gps_idx} = Y_gps{gps_idx}/1e3;
    end

    figure(fig(plot_idx));
    h_gps = plot(X_gps{gps_idx},Y_gps{gps_idx},'-b');
  end

  figure(fig(plot_idx));
  h_all = plot(X,Y,'-r');

  figure(fig(plot_idx));
  h_good = plot(X_good,Y_good,'.g');

  figure(fig(plot_idx));
  h_mod = plot(X_mod,Y_mod,'.y');

  figure(fig(plot_idx));
  h_det = plot(X_det,Y_det,'.r');
  
% Add single season plot labels
  figure(fig(plot_idx));
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

  print(fig(plot_idx),'-djpeg','-r300',image_fn);
end

return;
