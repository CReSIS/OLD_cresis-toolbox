% script flight_tracker
%
% Flight tracker which reads in a NMEA stream (GPGGA) from a serial
% device and plots the result on a geotiff.
%
% Optionally plots a KML flight line (this part is customized
% to NASA ATM's flightplan). Leave kml_fn empty to turn off
% this feature.
%
% Close plotting window or press ctrl-C (multiple times rapidly!)
% to quit.  Serial device is stored in global variable,
% flight_tracker_serial_dev, so that you can close it if the program
% fails to do so.
%
% Requires Mapping Toolbox, read_xml.m, plot_geotiff.m
%
% Examples:
%   flight_tracker
%
% Author: John Paden

% =================================================================
% User Settings
% =================================================================
% geotiff_fn = 'C:\GIS_data\greenland\Landsat-7\mzl7geo_90m_lzw.tif';
% geotiff_fn = 'C:\GIS_data\arctic\NaturalEarth_Data\Arctic_NaturalEarth.tif';
% geotiff_fn = '/scratch/GIS_data/greenland/Landsat-7/mzl7geo_90m_lzw.tif'; % For Land Ice
% geotiff_fn = '/scratch/GIS_data/greenland/Landsat-7/Greenland_natural_90m.tif'; % For Land Ice
geotiff_fn = '/scratch/GIS_data/arctic/Landsat-7/arctic_natural_90m.tif'; % For Land Ice
% geotiff_fn = '/scratch/GIS_data/arctic/ArcticDEM/arcticdem_mosaic_500m_v3.0.tif'; % For Land Ice
% geotiff_fn = '/scratch/GIS_data/arctic/NaturalEarth_Data/Arctic_NaturalEarth.tif'; % For Sea Ice

dem_fn = '/scratch/GIS_data/arctic/ArcticDEM/arcticdem_mosaic_500m_v3.0.tif';
[dem_RGB, dem_R, ~] = geotiffread(dem_fn);
dem_proj = geotiffinfo(dem_fn);
dem_x_axis = dem_R(3,1) + dem_R(2,1)*(1:size(dem_RGB,2));
dem_y_axis = dem_R(3,2) + dem_R(1,2)*(1:size(dem_RGB,1));

gps_input_type = 'serial'; % file_mcords, file_accum, or serial

% You may have to run from bash shell: "sudo chmod a+rwx /dev/ttyS0" as root
serial_dev = '/dev/ttyUSB0';
% serial_dev = '/dev/ttyUSB1';
% serial_dev = '/dev/ttyS0';
BAUD_RATE = 115200;

gps_input_fn_dir = '/tmp/';
gps_input_fn_start = 'GPS';
% gps_input_fn_dir = 'E:\';
% gps_input_fn_start = 'GPS';
% gps_input_fn_dir = '\\172.18.1.33\accum\';
% gps_input_fn_dir = '/net/field1/landing/mcords/mcords5/';
gps_input_fn_skip = false; % Enables skipping reading old data, sometimes
                           % need to do this if files contain errors that
                           % cause program to crash

% For OIB, get kmz file from John Sonntag, then unzip the kmz file and use the "doc.kml" that is inside
% kml_fn = 'C:\Users\administrator\Desktop\doc.kml';
kml_fn = '/scratch/metadata/2022_Greenland_P3/flight_plans/Eureka_flight_coords_04_22_2022_GT2L.kml'; % Set this to empty if the file is not available
kml_mission_name = '';

enable_gps_record = false;
gps_fn_dir = '/scratch/metadata/2022_Greenland_P3/';


[year,month,day] = datevec(now);

% =================================================================
% Automated Section
% =================================================================

if isempty(kml_fn)
  kml_lon = [];
  kml_lat = [];
else
  if 1
    % Simple KML file
    xDoc = xmlread(kml_fn);
    document = read_xml(xDoc);
    pos = textscan(document.kml{1}.Document{1}.Placemark{1}.LineString{1}.coordinates{1}.text{1}.node_val, ...
      '%f%f%f','Delimiter',',');
    kml_lon = pos{1};
    kml_lat = pos{2};
  else
    % KML file with many named lines in it (allows for selecting one of the lines)
    kml_pos = kml_read_shapefile(kml_fn);
    if isempty(kml_mission_name)
      for idx=1:length(kml_pos)
        fprintf('%d: %s\n', idx, kml_pos(idx).name);
      end
      kml_mission_idx = input('Enter mission number: ');
    else
      kml_mission_idx = find(strcmpi(kml_mission_name,{kml_pos.name}),1);
      if isempty(kml_mission_idx)
        warning('Could not find mission "%s" in kml_pos.name .',kml_mission_name);
      end
    end
    if isempty(kml_mission_idx)
      kml_mission_name = '';
      kml_fn = '';
      kml_lon = [];
      kml_lat = [];
    else
      kml_lon = kml_pos(kml_mission_idx).X;
      kml_lat = kml_pos(kml_mission_idx).Y;
    end
  end
end

if strcmpi(gps_input_type,'serial')
  global flight_tracker_serial_dev;
  if isempty(flight_tracker_serial_dev)
    fprintf('Getting serial device %s handle\n', serial_dev);
    flight_tracker_serial_dev = serial(serial_dev,'BaudRate',BAUD_RATE);
  end
  
  if ~strcmpi(get(flight_tracker_serial_dev,'Status'),'open')
    fprintf('Opening serial device %s\n', serial_dev);
    fopen(flight_tracker_serial_dev);
  else
    fprintf('Serial device %s already open, just going to start reading\n', ...
      serial_dev);
  end
  
  pos_buf = NaN*zeros(1e5,2);
  time_buf = NaN*zeros(10,1);
  lat_buf = NaN*zeros(length(time_buf),1);
  lon_buf = NaN*zeros(length(time_buf),1);
  elev_buf = NaN*zeros(length(time_buf),1);
end

fprintf('Plotting geotiff\n');
[proj,fig_h] = geotiff_plot(geotiff_fn, kml_lat, kml_lon, 1,'r');
haxes = get(fig_h,'Children');

if any(strcmpi(gps_input_type,{'file_accum','file_mcords'}))
  % Look for, load, and plot all GPS files
  if strcmpi(gps_input_type,'file_accum')
    gps_input_fn_ext = '.gps';
  else
    gps_input_fn_ext = '.txt';
  end
  
  % Monitoring:
  % Look for the latest GPS file
  % If there is no file, skip steps
  % If newest file matches the current file, then keep current position
  % Load latest file and search for '$'
  if ~gps_input_fn_skip
    gps = struct('lat',[],'lon',[]);
    gps_in_fns = get_filenames(gps_input_fn_dir,'','',gps_input_fn_ext);
    if isempty(gps_in_fns)
      warning('No GPS files in %s\n', gps_input_fn_dir);
    else
      gps_param = struct('year',year,'month',month,'day',day,'time_reference','utc','format',3);
      for gps_in_fn_idx = 1:length(gps_in_fns)
        gps_in_fn = gps_in_fns{gps_in_fn_idx};
        fid = fopen(gps_in_fn,'r');
        while ~feof(fid)
          line_input = fgets(fid);
          if feof(fid)
            break;
          end
          A = textscan(line_input,'%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s%f%f%f%f','delimiter',', ','emptyvalue',NaN);
          try
            if all(~cellfun(@isempty,A([1 4 5 6]))) && strcmp(A{1},'$GPGGA') && ~isnan(A{3}) && any(strcmpi(A{4},{'N','S'})) && ~isnan(A{5}) && any(strcmpi(A{6},{'W','E'}))
              gps.lat(1,end+1) = ((A{4}=='N')*2-1) .* A{3};
              gps.lon(1,end+1) = ((A{6}=='E')*2-1) .* A{5};
            end
          end
        end
        fclose(fid);
      end
      gps.lat = fix(gps.lat/100) + (gps.lat/100 - fix(gps.lat/100))./60*100;
      gps.lon = fix(gps.lon/100) + (gps.lon/100 - fix(gps.lon/100))./60*100;
    end
  end
  
  gps_in_fn = '';
  gps_in_fn_pos = -inf;
  
  pos_buf = NaN*zeros(1e5,2);
  time_buf = NaN*zeros(10,1);
  lat_buf = NaN*zeros(length(time_buf),1);
  lon_buf = NaN*zeros(length(time_buf),1);
  elev_buf = NaN*zeros(length(time_buf),1);
  
  if ~gps_input_fn_skip
    [x,y] = projfwd(proj,gps.lat,gps.lon);
    x = x/1e3;
    y = y/1e3;
    idx_to_use = max(1,length(x)-size(pos_buf,1)+1) : length(x);
    pos_buf(1:length(idx_to_use),1) = fliplr(x(idx_to_use));
    pos_buf(1:length(idx_to_use),2) = fliplr(y(idx_to_use));
  end
end

hold on;
hline = plot(pos_buf(:,1),pos_buf(:,2),'b-','Parent',haxes);
hpos = plot(pos_buf(1,1),pos_buf(1,2),'rx','MarkerSize',10,'LineWidth',3,'Parent',haxes);
hold off;

if enable_gps_record
  gps_fn = fullfile(gps_fn_dir,sprintf('gps_%04d%02d%02d.csv',year,month,day));
  if exist(gps_fn,'file')
    try
      gps = read_gps_csv(gps_fn, struct('time_reference','utc'));
      [x,y] = projfwd(proj,gps.lat,gps.lon);
      x = x/1e3;
      y = y/1e3;
      idx_to_use = max(1,length(x)-size(pos_buf,1)+1) : length(x);
      pos_buf(1:length(idx_to_use),1) = fliplr(x(idx_to_use));
      pos_buf(1:length(idx_to_use),2) = fliplr(y(idx_to_use));
      set(hline,'XData',pos_buf(:,1));
      set(hline,'YData',pos_buf(:,2));
      set(hpos,'XData',pos_buf(1,1));
      set(hpos,'YData',pos_buf(1,2));
      drawnow;
      pos_buf(2:end,:) = pos_buf(1:end-1,:);
    catch
    end
  end
  
  % ===========================================================
  %% Opening GPS log file
  % ===========================================================
  fprintf('Opening GPS log file %s\n', gps_fn);
  [fid_out,msg] = fopen(gps_fn,'a');
  if fid_out <= 0
    error(msg);
  end
  if ftell(fid_out) == 0
    % Write header line if file is empty
    fprintf('  New file, writing CSV header line\n');
    fprintf(fid_out,'year,month,day,UTC_sod,latNdeg,lonEdeg,elevm\n');
  end
end

update_geotif_tstart = uint64(0);

try
  done = false;
  while ~done
    %% Update Map Axes
    xlim_orig = xlim(haxes);
    ylim_orig = ylim(haxes);
    axis normal;
    % Force axis to be equal
    xlim_new = xlim();
    ylim_new = ylim();
    axes_pos = get(haxes,'position'); aspect_ratio = axes_pos(3)/axes_pos(4);
    xlim_bigger = diff(xlim_new)/diff(ylim_new)/aspect_ratio > 1;
    if xlim_bigger
      ylim_new(1) = mean(ylim_new) - diff(xlim_new)/2*aspect_ratio;
      ylim_new(2) = mean(ylim_new) + diff(xlim_new)/2*aspect_ratio;
    else
      xlim_new(1) = mean(xlim_new) - diff(ylim_new)/2/aspect_ratio;
      xlim_new(2) = mean(xlim_new) + diff(ylim_new)/2/aspect_ratio;
    end
    if toc(update_geotif_tstart) > 30
      update_geotif_tstart = tic;
      if isfinite(pos_buf(1,1)) && isfinite(pos_buf(1,2))
        % Update to current platform position
        xlim_new = xlim_new + pos_buf(1,1) - mean(xlim_new);
        ylim_new = ylim_new + pos_buf(1,2) - mean(ylim_new);
      end
    end
    xlim(haxes,xlim_new);
    ylim(haxes,ylim_new);
    
    %% Load GPS data from serial
    if strcmpi(gps_input_type,'serial')
      try
        nmea_str = fscanf(flight_tracker_serial_dev);
      catch ME
        pause(0.5);
        continue;
      end
      if isempty(nmea_str)
        fprintf('Empty string\n');
        pause(0.5);
        continue;
      end
      A = textscan(nmea_str,'%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s','delimiter',',','emptyvalue',NaN);
    end
    
    %% Load GPS data from file
    if any(strcmpi(gps_input_type,{'file_accum','file_mcords'}))
      % Look for, load, and plot all GPS files
      if strcmpi(gps_input_type,'file_accum')
        gps_input_fn_ext = '.gps';
      else
        gps_input_fn_ext = '.txt';
      end
      % Check to see if we are looking at the most recent GPS file
      gps_in_fns = get_filenames(gps_input_fn_dir,gps_input_fn_start,'',gps_input_fn_ext);
      if ~isempty(gps_in_fns)
        if ~strcmpi(gps_in_fn,gps_in_fns{end})
          % New GPS file
          gps_in_fn = gps_in_fns{end};
          
          % Find the last dollar sign '$' in the file
          fid = fopen(gps_in_fn,'r');
          fseek(fid,0,1);
          dollar_found = false;
          while ~dollar_found
            if ftell(fid) == 0
              % Last GPS file does not contain '$'
              gps_in_fn_pos = -inf;
              break;
            end
            fseek(fid,-1,0);
            A = fread(fid,1,'char');
            if A == '$'
              gps_in_fn_pos = ftell(fid)-1;
              dollar_found = true;
            else
              fseek(fid,-1,0);
            end
          end
          fclose(fid);
        end
        
        if isfinite(gps_in_fn_pos)
          % We have a valid file position handle
          pause(0.5);
          fid = fopen(gps_in_fn,'r');
          fseek(fid,gps_in_fn_pos,-1);
          nmea_str = fgets(fid);
          A = fread(fid,1,'char');
          if isempty(A) || A(1) ~= '$'
            gps_in_fn_pos = ftell(fid)-1;
            A = {''};
          else
            gps_in_fn_pos = ftell(fid)-1;
            A = textscan(nmea_str,'%s%f%f%c%f%c%u%u%f%f%c%f%c%s%s%f%f%f%f','delimiter',',','emptyvalue',NaN);
          end
          fclose(fid);
        else
          % We do not have valid data
          A = {''};
        end
      end
    end
      
    if strcmpi(A{1},'$GPGGA')
      lat_deg = floor(A{3}/100);
      lat_min = A{3}-lat_deg*100;
      lat = lat_deg + lat_min/60;
      if A{4} ~= 'N'
        lat = -lat;
      end
      lon_deg = floor(A{5}/100);
      lon_min = A{5}-lon_deg*100;
      lon = lon_deg + lon_min/60;
      if A{6} ~= 'E'
        lon = -lon;
      end
      hour = floor(A{2}/1e4);
      minute = floor((A{2}-hour*1e4)/1e2);
      sec= A{2}-hour*1e4-minute*1e2;
      
      elev = A{10};
      
      utc_time = datenum(year,month,day,hour,minute,sec);
      day_start = datenum(year,month,day,0,0,0);
      utc_sod = (utc_time - day_start)*86400;
      
      [dem_x,dem_y] = projfwd(dem_proj,lat,lon);
      dem_elev = interp2(dem_x_axis,dem_y_axis,dem_RGB,dem_x,dem_y);

      title(sprintf('UTC %02d:%02d:%05.2f, %.2f SOD, %.0f m/%.0f ft, %.0f m/%.0f ft',hour,minute,sec, ...
        utc_sod, elev, elev*100/2.54/12, elev-dem_elev, (elev-dem_elev)*100/2.54/12));
      
      if isempty(lat)
        lat = NaN;
      end
      if isempty(lon)
        lon = NaN;
      end
      if isempty(elev)
        elev = NaN;
      end
      time_buf(2:end) = time_buf(1:end-1);
      time_buf(1) = utc_sod;
      lat_buf(2:end) = lat_buf(1:end-1);
      lat_buf(1) = lat;
      lon_buf(2:end) = lon_buf(1:end-1);
      lon_buf(1) = lon;
      elev_buf(2:end) = elev_buf(1:end-1);
      elev_buf(1) = elev;
      along_track = geodetic_to_along_track(lat_buf([1 end]),lon_buf([1 end]),elev_buf([1 end]));
      speed = along_track(2) / abs(diff(time_buf([1 end])));
      
      fprintf('%7.1f | %9.6f N %11.6f E | %6.1f m | %3.0f m/s %3.0f kn\n', utc_sod, lat, lon, elev, speed, speed/0.5144444);
      
      if enable_gps_record
        fprintf(fid_out,'%04d,%02d,%02d,%f,%f,%f,%f\n', year, month, day, utc_sod, lat, lon, elev);
      end
      
      [x,y] = projfwd(proj,lat,lon);
      x = x/1e3;
      y = y/1e3;
      pos_buf(1,1) = x;
      pos_buf(1,2) = y;
      set(hline,'XData',pos_buf(:,1));
      set(hline,'YData',pos_buf(:,2));
      set(hpos,'XData',pos_buf(1,1));
      set(hpos,'YData',pos_buf(1,2));
      drawnow;
      pos_buf(2:end,:) = pos_buf(1:end-1,:);
    end
  end
catch ME
  ME
  ME.stack(1)
end

if strcmpi(gps_input_type,'serial')
  fprintf('Closing serial device %s\n', serial_dev);
  fclose(flight_tracker_serial_dev);
end
if enable_gps_record
  fprintf('Closing gps log file %s\n', gps_fn);
  fclose(fid_out);
end

return;
