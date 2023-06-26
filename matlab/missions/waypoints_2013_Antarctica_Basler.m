if 0
  
% Plot multiple flight lines
% Listbox of each flight line
% Listbox of active waypoints for selected flight line
% Ctrl-click to select active flight line
% Click to select active waypoint
% Right click to "delete","insert before", "insert after" current
% selection and can have multiple points selected. Insert brings
% up listbox to choose active flight line.
%
% Opens window with specified geotiff
% Right click in flight lines box to "Load", "New"
% Tool listbox to "Select","Insert"
% Scroll button zooms in/out on map
% Left click applies tool
%
% Read .csv file in lat,lon format
% Read .csv file in lat,lon,name format
%   Load in test flight 1
%   
% Read .shp file in lat,lon format
%   Load in all Whillans shapes

% lat,lon,name format

fns = {};
fns{end+1} = 'C:\Users\dangermo\Documents\Travel\Antarctica_2013\cresis_flight_plans\Test_Flight_1_2013.csv';
fns{end+1} = 'C:\Users\dangermo\Documents\Travel\Antarctica_2013\cresis_flight_plans\Test_Flight_2_2013.csv';

%
  fid = fopen(fn,'r');
  C = textscan(fid,'%f %f %s','Delimiter',',','HeaderLines',1);
  fclose(fid);
  [lat,lon,name] = C{:};

  save(way_point_fn,'lat','lon','name');

%   fig_h = 1;
%   figure(fig_h); clf;
%   proj = plot_geotiff('C:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif',lat,lon,fig_h);
  
  modify_waypoints('C:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif',lat,lon);
  return;
end

fp_dir = 'C:\Users\dangermo\Documents\Travel\Antarctica_2013\cresis_flight_plans';

fns = get_filenames(fp_dir,'','','.shp');

fig_h = 1;
colors = {'r','g','b'};
figure(fig_h); clf;
proj = plot_geotiff('C:\GIS_data\antarctica\Landsat-7\Antarctica_LIMA_480m.tif',[],[],fig_h);

for fn_idx = 3%1:length(fns)
  fn = fns{fn_idx};
  fprintf('%s\n', fn);
  
  S = shaperead(fn);
  S_info = shapeinfo(fn);
  S_geo = shaperead(fn,'UseGeoCoords',true);
  
  if isfield(S,'Flt_ln_ID')
    for idx = 1:length(S)
      S(idx).flt_ln_ID = S(idx).Flt_ln_ID;
    end
  end
  
  % Find unique days
  days = [];
  day_mapping = zeros(size(S));
  for idx = 1:length(S)
    if ~isnumeric(S(idx).Day)
      keyboard;
    end
    if ~any(S(idx).Day == days)
      days = [days S(idx).Day];
    end
  end
  day_mapping = cell2mat({S.Day});
  [days_sorted day_idxs] = sort(days);
  
  % Sort each day
  for day_idx = 1:8%length(days_sorted)
    day = days_sorted(day_idx);
    day_idxs = find(day_mapping == day);
    turn_idxs = [];
    turn_ids = [];
    straight_idxs = [];
    straight_ids = [];
    debug_plot = false;
    if debug_plot
      figure(1001); clf;
    end
    for pnt = 1:length(day_idxs)
      pnt_idx = day_idxs(pnt);
      % Look for decimal
      flt_ln_ID = S(pnt_idx).flt_ln_ID;
      flt_ln_dec_idx = find(flt_ln_ID == '.');
      flt_ln_num = str2double(flt_ln_ID);
      if ~isempty(flt_ln_dec_idx)
        % This is a turn
        flt_ln_num = str2double(flt_ln_ID(flt_ln_dec_idx+1:end));
        if debug_plot
          fprintf('(%s): This is a turn (%d)\n', flt_ln_ID, flt_ln_num);
        end
        turn_idxs = [turn_idxs,pnt_idx];
        turn_ids = [turn_ids flt_ln_num];
      elseif flt_ln_num >= 1
        if debug_plot
          fprintf('(%s): This is a straight (%d)\n', flt_ln_ID,flt_ln_num);
        end
        straight_idxs = [straight_idxs,pnt_idx];
        straight_ids = [straight_ids flt_ln_num];
      else
        if strcmpi(flt_ln_ID,'Start')
          start_idx = pnt_idx;
        else
          stop_idx = pnt_idx;
        end
        if debug_plot
          fprintf('(%s): This ingress/egress\n', flt_ln_ID);
        end
      end
      
      if debug_plot
        plot(S(pnt_idx).X,S(pnt_idx).Y,'r');
        hold on;
        pause;
        plot(S(pnt_idx).X,S(pnt_idx).Y);
      end
    end
    
    % Now plot in order
    [turn_ids sort_idxs] = sort(turn_ids);
    turn_idxs = turn_idxs(sort_idxs);
    [straight_ids sort_idxs] = sort(straight_ids);
    straight_idxs = straight_idxs(sort_idxs);
    if length(straight_idxs) > length(turn_idxs)
      combined_idxs = reshape([straight_idxs; [turn_idxs 0]],[1 2*length(straight_idxs)]);
      ordered_idxs = [start_idx combined_idxs(1:end-1) stop_idx];
    elseif length(straight_idxs) < length(turn_idxs)
      combined_idxs = reshape([turn_idxs; [straight_idxs 0]],[1 2*length(turn_idxs)]);
      ordered_idxs = [start_idx combined_idxs(1:end-1) stop_idx];
    else
      start_x = S(start_idx).X(~isnan(S(start_idx).X));
      start_y = S(start_idx).Y(~isnan(S(start_idx).Y));
      turn_dist = sqrt( (start_x(end) - S(turn_idxs(1)).X(1)).^2 ...
        + (start_y(end) - S(turn_idxs(1)).Y(1)).^2 );
      straight_dist = sqrt( (start_x(end) - S(straight_idxs(1)).X(1)).^2 ...
        + (start_y(end) - S(straight_idxs(1)).Y(1)).^2 );
      if straight_dist < turn_dist
        combined_idxs = reshape([straight_idxs; turn_idxs],[1 2*length(turn_idxs)]);
      else
        combined_idxs = reshape([turn_idxs; straight_idxs],[1 2*length(turn_idxs)]);
      end
      ordered_idxs = [start_idx combined_idxs stop_idx];
    end
    %type_mask = reshape([ones(size(straight_idxs)); [turn_idxs 0]],[1 2*length(straight_idxs)]);
    
    flight_length = 0;
    num_pnts = 0;
    X = [];
    Y = [];
    lat = [];
    lon = [];
    for pnt = 1:length(ordered_idxs)
      flight_length = flight_length + S(pnt_idx).Length;
      pnt_idx = ordered_idxs(pnt);
      title(sprintf('%s', S(pnt_idx).flt_ln_ID));
      X_data = S(pnt_idx).X(~isnan(S(pnt_idx).X));
      Y_data = S(pnt_idx).Y(~isnan(S(pnt_idx).X));
      Lat_data = S_geo(pnt_idx).Lat(~isnan(S(pnt_idx).X));
      Lon_data = S_geo(pnt_idx).Lon(~isnan(S(pnt_idx).X));
      if any(isnan(X_data))
        warning('NaN found');
      end
      if any(isnan(Y_data))
        warning('NaN found');
      end
      if any(isnan(Lat_data))
        warning('NaN found');
      end
      if any(isnan(Lon_data))
        warning('NaN found');
      end
      figure(fig_h);
      hold on;
      plot(X_data/1e3,Y_data/1e3,'r');
      %pause;
      plot(X_data/1e3,Y_data/1e3);
      X = [X,X_data];
      Y = [Y,Y_data];
      lat = [lat,Lat_data];
      lon = [lon,Lon_data];
      
    end
    along_track = geodetic_to_along_track(lat,lon,zeros(size(lat)));
    fprintf(' Day %d: Total length %f km (%f km)\n', day, flight_length, along_track(end)/1e3);
    
  end
end

return
