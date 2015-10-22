function plot_vectors(filename,bbox,geoTiff)
% plot_vectors(filename,bbox,geoTiff)
%
% filename of file created by create_seasonVector
%    (can also be a cell array of filenames to load)
% bbox = bounding boxes (4 by N). These are two (lat/lon) pairs in deg.
%   where N is the number of bounding boxes you want displayed
%   as long as opposite corners are specified it will work:
%   [corner1_lat; corner1_lon; corner2_lat; corner2_lon]
% lat,lon pair provides the start point
%
% Takes .mat file saved by Anthony's create_seasonVectorX.m functions
% and plots and makes a convenient graphical user interface
%
% Commands (only work when not in zoom/pan/special modes):
%   q = quit
%   left mouse button = get file numbers at this location
%   z = zoom mode
%   k = keyboard/debug mode
%
% Examples:
%
% See run_plot_vectors.m
%
% Author: John Paden
%
% See also: create_season_vectors_snow.m, create_season_vectors_mcords.m,
%   run_plot_vectors.m

global h1 h1_axis;
global h1_highlight;
global geobase;
geobase = [];

% Check arguments
if exist('bbox','var') && ~isempty(bbox)
  if bbox(1) > bbox(3)
    tmp = bbox(3);
    bbox(3) = bbox(1);
    bbox(1) = tmp;
  end
  if bbox(2) > bbox(4)
    tmp = bbox(4);
    bbox(4) = bbox(2);
    bbox(2) = tmp;
  end
else
  bbox = [];
end

if iscell(filename)
  load(filename{1});
  [fn_dir fn_name] = fileparts(filename{1});
  for lat_idx=1:length(vectors.lat)
    vectors.day_seg{lat_idx} = fn_name(end-10:end);
  end
  
  for idx = 2:length(filename)
    A = load(filename{idx});
    [fn_dir fn_name] = fileparts(filename{idx});
    for lat_idx=1:length(A.vectors.lat)
      vectors.day_seg{end+1} = fn_name(end-10:end);
    end
    vectors.lon = [vectors.lon A.vectors.lon];
    vectors.lat = [vectors.lat A.vectors.lat];
    vectors.elev = [vectors.elev A.vectors.elev];
    if isfield(vectors,'file')
      vectors.file = [vectors.file A.vectors.file];
    end
    vectors.fileNumber = [vectors.fileNumber A.vectors.fileNumber];
    vectors.gps_time = [vectors.gps_time A.vectors.gps_time];
  end
else
  load(filename);
  [fn_dir fn_name] = fileparts(filename);
  for lat_idx=1:length(vectors.lat)
    vectors.day_seg{lat_idx} = fn_name(end-10:end);
  end
end

% Rename vectors to geobase
geobase = vectors;
clear vectors;
geobase.Longitude = geobase.lon;
geobase = rmfield(geobase,'lon');
geobase.Latitude = geobase.lat;
geobase = rmfield(geobase,'lat');
geobase.line = 1;
geobase.plot_color = 'b-x';

if exist('geoTiff','var') && ~isempty(geoTiff)
  % Obtain the projection structure.
  geobase.proj = geotiffinfo(geoTiff);
  % Read the image
  [RGB, R, tmp] = geotiffread(geoTiff);
  if size(RGB,3) == 3 && strcmp(class(RGB),'uint16') && max(RGB(:)) <= 255
    RGB = uint8(RGB);
  end
  R = R/1e3;
  % Select all data within this many kilometers
  geobase.dist = 2;
else
  % Select all data within this many degrees
  geobase.dist = 0.1;
end

% Add file numbers only when they are separated by this many file numbers
geobase.dist_fn = 50;

% ===================================================================
% Plot geo data
h1 = figure(1000); clf;
zoom on;
zoom off;
if isfield(geobase,'proj')
  mapshow(RGB, R);
end

if isfield(geobase,'proj')
  [geobase.X,geobase.Y] = projfwd(geobase.proj,geobase.Latitude,geobase.Longitude);
  geobase.X = geobase.X/1e3;
  geobase.Y = geobase.Y/1e3;
  if ~isempty(bbox)
    [geobase.bbox.X,geobase.bbox.Y] = projfwd(geobase.proj,bbox([1 3]),bbox([2 4]));
    geobase.bbox.X = geobase.bbox.X/1e3;
    geobase.bbox.Y = geobase.bbox.Y/1e3;
  else
    geobase.bbox = [];
  end
else
  geobase.X = geobase.Longitude;
  geobase.Y = geobase.Latitude;
  if ~isempty(bbox)
    geobase.bbox.X = bbox([2 4]);
    geobase.bbox.Y = bbox([1 3]);
  else
    geobase.bbox = [];
  end
end
% Plot geocoordinates
figure(h1);
hold on;
geobase.h1 = plot(geobase.X,geobase.Y,geobase.plot_color);
hold off;
grid on;
axis tight;
h1_axis = gca;
if ~isempty(geobase.bbox)
  % Plot bbox
  figure(h1);
  hold on;
  plot(geobase.bbox.X([1 2 2 1 1]),geobase.bbox.Y([1 1 2 2 1]),'y');
  hold off;
end

% ===================================================================
% Set highlights
figure(h1);
zoom reset
geobase.axis = axis;
if isfield(geobase,'proj')
  xlabel('Easting (km)');
  ylabel('Northing (km)');
else
  xlabel('Longitude (deg)');
  ylabel('Latitude (deg)');
end
hold on;
h1_highlight = plot( ...
  geobase.X(geobase.line), ...
  geobase.Y(geobase.line),'ro');
set(h1_highlight, 'MarkerSize', 9);
set(h1_highlight, 'LineWidth', 1.5);
hold off;

% Set mouse pointer
set(h1,'Pointer','crosshair');

% Make sure special tools are off before assigning callbacks
zoom on;
zoom off;

% Set callback properties
set(h1,'Interruptible','off');
set(h1,'WindowButtonUpFcn',@gui_windowbuttonupfcn)
set(h1,'WindowButtonMotionFcn',[])
set(h1,'WindowKeyPressFcn',@gui_windowkeypressfcn);

set(h1,'UserData',0);
while 1
  waitfor(h1,'UserData');
  if ~ishandle(h1)
    % User closed the figure
    break;
  end
  if get(h1,'UserData') == 1
    axis tight;
    set(gca,'XLim',geobase.axis(1:2));
    set(gca,'YLim',geobase.axis(3:4));

    % Unset callback properties
    set(h1,'Interruptible','on');
    set(h1,'WindowButtonUpFcn',[]);
    set(h1,'WindowKeyPressFcn',[]);

    fprintf('Press any key to exit zoom mode...\n');
    zoom on;
    pause;
    fprintf('\n  done\n');
    geobase.axis = axis;
    zoom on;
    zoom off;

    set(h1,'UserData',0);
    
    % Set callback properties
    set(h1,'Interruptible','off');
    set(h1,'WindowButtonUpFcn',@gui_windowbuttonupfcn)
    set(h1,'WindowKeyPressFcn',@gui_windowkeypressfcn);
  elseif get(h1,'UserData') == 2
    % Quitting
    break;
  end
end
  
return;

% =====================================================================
% WindowButtonUpFcn call back function
% =====================================================================
function gui_windowbuttonupfcn(src,event)

global h1 h1_axis;
global h1_highlight;
global geobase;

[x,y,but] = get_mouse_info(h1,h1_axis);
%fprintf('1: x = %.2f, y = %.2f, but = %d\n', x, y, but);

% Find the closest flight trajectory
dist = sqrt((geobase.X-x).^2 + (geobase.Y-y).^2);
[tmp inds] = sort(dist);
geobase.line = inds(1);

val_file_numbers = geobase.fileNumber(inds(1));
clear val_file_info;
[year month day] = datevec(epoch_to_datenum(geobase.gps_time(inds(1))));
val_file_info{1} = sprintf('%s:%04d%02d%02d:%04d', geobase.day_seg{inds(1)}, ...
  year, month, day, geobase.fileNumber(inds(1)));
ind = 2;
while ind < length(inds) && dist(inds(ind)) < geobase.dist
  if ~too_close_check(val_file_numbers,geobase.fileNumber(inds(ind)),geobase.dist_fn)
    val_file_numbers(end+1) = geobase.fileNumber(inds(ind));
    [year month day] = datevec(epoch_to_datenum(geobase.gps_time(inds(ind))));
    val_file_info{end+1} = sprintf('%s:%04d%02d%02d:%04d',  geobase.day_seg{inds(ind)}, ...
      year, month, day, geobase.fileNumber(inds(ind)));
  end
  ind = ind + 1;
end
fileNumber_str = sprintf('%s; ', val_file_info{1:length(val_file_numbers)});
figure(h1);
if isfield(geobase,'file')
  title(sprintf('%s\n%.3f N %.3f E, GPS %s FILE %s', fileNumber_str, ...
    geobase.Latitude(inds(1)), geobase.Longitude(inds(1)), ...
    datestr(epoch_to_datenum(geobase.gps_time(inds(1))),'HH:MM:SS'), ...
    datestr(geobase.file(inds(1)).datenum,'dd HH:MM:SS') ), 'Interpreter','none');
  geobase.file(inds(1))
  fprintf('File: %s\n', datestr(geobase.file(inds(1)).datenum,'yyyymmdd HH:MM:SS'))
else
  title(sprintf('%s\n%.3f N %.3f E, %s GPS', fileNumber_str, ...
    geobase.Latitude(inds(1)), geobase.Longitude(inds(1)), ...
    datestr(epoch_to_datenum(geobase.gps_time(inds(1))),'HH:MM:SS')));
end

% ===================================================================
% Set highlights
figure(h1);
delete(h1_highlight);
hold on;
h1_highlight = plot( ...
  geobase.X(geobase.line), ...
  geobase.Y(geobase.line),'ro');
set(h1_highlight, 'MarkerSize', 7);
set(h1_highlight, 'LineWidth', 1.5);
hold off;

return;

% =====================================================================
% WindowKeyPressFcn call back function
% =====================================================================
function gui_windowkeypressfcn(src,event)

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if ~isempty(event.Character)

  global h1 h1_axis;
  global h1_highlight;
  global geobase;

  [x,y,but] = get_mouse_info(h1,h1_axis);

  if 0 % This code for debugging
    if ischar(event.Key)
      fprintf('x = %.2f, y = %.2f, key = %s\n', x, y, event.Key);
    else
      fprintf('x = %.2f, y = %.2f, key = %d\n', x, y, event.Key);
    end
    if ~isempty(event.Modifier)
      fprintf('  Modifiers ');
      for ind = 1:length(event.Modifier)
        fprintf('%s ', event.Modifier{ind});
      end
      fprintf('\n');
    end
  end

  % see event.Modifier for modifiers
  switch event.Character
    case 'z'
      set(h1,'UserData',1);
    case 'q'
      fprintf('Quitting\n');
      figure(h1);
      title('');
      set(h1,'UserData',2);
      % Set mouse pointer
      set(h1,'Pointer','arrow');

      % Set callback properties
      set(h1,'Interruptible','on');
      set(h1,'WindowButtonUpFcn',[]);
      set(h1,'WindowKeyPressFcn',[]);
      
    case 28 % Left-arrow
      dist = 0;
      if geobase.line > 1
        geobase.line = geobase.line - 1;
      end
      inds = geobase.line;
      
      val_file_numbers = geobase.fileNumber(inds(1));
      clear val_file_info;
      [year month day] = datevec(epoch_to_datenum(geobase.gps_time(inds(1))));
      val_file_info{1} = sprintf('%s:%04d%02d%02d:%04d', geobase.day_seg{inds(1)}, ...
        year, month, day, geobase.fileNumber(inds(1)));
      ind = 2;
      while ind < length(inds) && dist(inds(ind)) < geobase.dist
        if ~too_close_check(val_file_numbers,geobase.fileNumber(inds(ind)),geobase.dist_fn)
          val_file_numbers(end+1) = geobase.fileNumber(inds(ind));
          [year month day] = datevec(epoch_to_datenum(geobase.gps_time(inds(ind))));
          val_file_info{end+1} = sprintf('%s:%04d%02d%02d:%04d',  geobase.day_seg{inds(ind)}, ...
            year, month, day, geobase.fileNumber(inds(ind)));
        end
        ind = ind + 1;
      end
      fileNumber_str = sprintf('%s; ', val_file_info{1:length(val_file_numbers)});
      figure(h1);
      if isfield(geobase,'file')
        title(sprintf('%s\n%.3f N %.3f E, GPS %s FILE %s', fileNumber_str, ...
          geobase.Latitude(inds(1)), geobase.Longitude(inds(1)), ...
          datestr(epoch_to_datenum(geobase.gps_time(inds(1))),'HH:MM:SS'), ...
          datestr(geobase.file(inds(1)).datenum,'dd HH:MM:SS') ), 'Interpreter','none');
        geobase.file(inds(1))
        fprintf('File: %s\n', datestr(geobase.file(inds(1)).datenum,'yyyymmdd HH:MM:SS'))
      else
        title(sprintf('%s\n%.3f N %.3f E, %s GPS', fileNumber_str, ...
          geobase.Latitude(inds(1)), geobase.Longitude(inds(1)), ...
          datestr(epoch_to_datenum(geobase.gps_time(inds(1))),'HH:MM:SS')));
      end
      
      % ===================================================================
      % Set highlights
      figure(h1);
      delete(h1_highlight);
      hold on;
      h1_highlight = plot( ...
        geobase.X(geobase.line), ...
        geobase.Y(geobase.line),'ro');
      set(h1_highlight, 'MarkerSize', 7);
      set(h1_highlight, 'LineWidth', 1.5);
      hold off;
  
    case 29 % Right-arrow
      dist = 0;
      if geobase.line < length(geobase.Latitude)
        geobase.line = geobase.line + 1;
      end
      inds = geobase.line;
      
      val_file_numbers = geobase.fileNumber(inds(1));
      clear val_file_info;
      [year month day] = datevec(epoch_to_datenum(geobase.gps_time(inds(1))));
      val_file_info{1} = sprintf('%s:%04d%02d%02d:%04d', geobase.day_seg{inds(1)}, ...
        year, month, day, geobase.fileNumber(inds(1)));
      ind = 2;
      while ind < length(inds) && dist(inds(ind)) < geobase.dist
        if ~too_close_check(val_file_numbers,geobase.fileNumber(inds(ind)),geobase.dist_fn)
          val_file_numbers(end+1) = geobase.fileNumber(inds(ind));
          [year month day] = datevec(epoch_to_datenum(geobase.gps_time(inds(ind))));
          val_file_info{end+1} = sprintf('%s:%04d%02d%02d:%04d',  geobase.day_seg{inds(ind)}, ...
            year, month, day, geobase.fileNumber(inds(ind)));
        end
        ind = ind + 1;
      end
      fileNumber_str = sprintf('%s; ', val_file_info{1:length(val_file_numbers)});
      figure(h1);
      if isfield(geobase,'file')
        title(sprintf('%s\n%.3f N %.3f E, GPS %s FILE %s', fileNumber_str, ...
          geobase.Latitude(inds(1)), geobase.Longitude(inds(1)), ...
          datestr(epoch_to_datenum(geobase.gps_time(inds(1))),'HH:MM:SS'), ...
          datestr(geobase.file(inds(1)).datenum,'dd HH:MM:SS') ), 'Interpreter','none');
        geobase.file(inds(1))
        fprintf('File: %s\n', datestr(geobase.file(inds(1)).datenum,'yyyymmdd HH:MM:SS'))
      else
        title(sprintf('%s\n%.3f N %.3f E, %s GPS', fileNumber_str, ...
          geobase.Latitude(inds(1)), geobase.Longitude(inds(1)), ...
          datestr(epoch_to_datenum(geobase.gps_time(inds(1))),'HH:MM:SS')));
      end
      
      % ===================================================================
      % Set highlights
      figure(h1);
      delete(h1_highlight);
      hold on;
      h1_highlight = plot( ...
        geobase.X(geobase.line), ...
        geobase.Y(geobase.line),'ro');
      set(h1_highlight, 'MarkerSize', 7);
      set(h1_highlight, 'LineWidth', 1.5);
      hold off;
    
  end
end

return;

% The idea is to plot just one file number from each
% flight line that is nearby a clicked point. Since
% many file numbers may be near by, this function
% makes sure that only one of the file numbers will
% be included.
function too_close = too_close_check(fileNumbers,fileNumber,dist)

too_close = false;
for ind = 1:length(fileNumbers)
  if abs(fileNumbers(ind) - fileNumber) < dist
    too_close = true;
    return;
  end
end

return;
