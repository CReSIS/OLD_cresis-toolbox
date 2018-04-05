function create_frames(param, param_override)
% create_frames(param, param_override)
%
% Takes .mat file saved by create_records_mcords.m function
% and allows the user to create frames which will be used by the
% CSARP processor.
%
% When the program loads, the entire dataset is treated as one
% frame.
%
% To create a new frame, select a position on the flight line
% and then press the create button or the 'c' key.
%
% To remove a break, delete the frame after the break by selecting
% it and pressing the delete button or the "backspace" key.
%
% To delete the last created frame, press the "backspace" key.
%
% To auto-generate frames, fill out the auto-generate fields and
% click on the "auto-generate" button.
%
% To quit, close the window
%
% records_fn: filename of file created by create_records_mcords
%    (can also be a cell array of filenames to load)
% frames_fn: output filename where framing will be stored, if this
%   file already exists, it will be loaded in
% geotiff_fn: filename of geotiff file, optional, can leave empty 
%
% Inputs:
% param = struct with processing parameters
%         -- OR --
%         function handle to script with processing parameters
% param_override = parameters in this struct will override parameters
%         in param.  This struct must also contain the gRadar fields.
%         Typically global gRadar; param_override = gRadar;
%
% Examples:
%
% create_frames(param);

%% General Setup
% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=====================================================================\n');

%% Setup creation of frames
% =====================================================================

[output_dir,radar_type,radar_name] = ct_output_dir(param.radar_name);

records_fn = ct_filename_support(param,'','records');
records_ver = load(records_fn,'ver');
if isfield(records_ver,'ver')
  records_file.records = load(records_fn);
else
  records_file = load(records_fn, 'records');
end
if ~isfield(records_file.records,'lat') || isempty(records_file.records.lat)
  fprintf('No geographic data present\n');
  frames_fn = ct_filename_support(param,'','frames');
  if exist(frames_fn,'file')
    fprintf('  Frames file already exists, just exiting\n');
  else
    fprintf('  No frames file exists, creating single frame and exiting\n');
    frames.frame_idxs = 1;
    frames.nyquist_zone = NaN;
    frames.proc_mode = 0;
    fprintf('  Saving %s\n', frames_fn);
    save(frames_fn,'frames');
  end
  return;
end

global hui; % hui: user interface handles
hui = [];
global GB; % Geobase: contains all the geographic info
GB = [];

hui.fig.handle = figure;
set(hui.fig.handle,'Name',param.day_seg);

%% Set up GUI
% =====================================================================
% =====================================================================
figure(hui.fig.handle); clf;

% Create widgets of main table
hui.fig.ctrl_panel.handle = uipanel('Parent',hui.fig.handle);
set(hui.fig.ctrl_panel.handle,'Title','Frame Controls');
set(hui.fig.ctrl_panel.handle,'TitlePosition','CenterTop');
set(hui.fig.ctrl_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(hui.fig.ctrl_panel.handle,'ShadowColor',[0.6 0.6 0.6]);

hui.fig.map_axes.handle = axes;

% Setup main table
hui.fig.table.ui = hui.fig.handle;
row = 1; col = 1;
hui.fig.table.handles{row,col}   = hui.fig.ctrl_panel.handle;
hui.fig.table.width(row,col)     = 220;
hui.fig.table.height(row,col)    = inf;
hui.fig.table.false_height(row,col) = 0;
col = col + 1;
hui.fig.table.handles{row,col}   = hui.fig.map_axes.handle;
hui.fig.table.width(row,col)     = inf;
hui.fig.table.height(row,col)    = inf;
hui.fig.table.false_height(row,col) = 0;
hui.fig.table.width_margin(row,col) = 60;
hui.fig.table.height_margin(row,col) = 40;
clear row col
table_draw(hui.fig.table);

hui.fig.ctrl_panel.framesLB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.framesLB,'Style','listbox');
set(hui.fig.ctrl_panel.framesLB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.framesLB,'FontName','fixed');
set(hui.fig.ctrl_panel.framesLB,'Callback',@framesLB_callback);

hui.fig.ctrl_panel.startLabel = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.startLabel,'Style','Text');
set(hui.fig.ctrl_panel.startLabel,'String','Start Record');
set(hui.fig.ctrl_panel.startLabel,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.startTB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.startTB,'Style','Edit');
set(hui.fig.ctrl_panel.startTB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.startTB,'Enable','off');

hui.fig.ctrl_panel.stopLabel = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.stopLabel,'Style','Text');
set(hui.fig.ctrl_panel.stopLabel,'String','Stop Record');
set(hui.fig.ctrl_panel.stopLabel,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.stopTB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.stopTB,'Style','Edit');
set(hui.fig.ctrl_panel.stopTB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.stopTB,'Enable','off');

hui.fig.ctrl_panel.procLabel = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.procLabel,'Style','Text');
set(hui.fig.ctrl_panel.procLabel,'String','Process Mode');
set(hui.fig.ctrl_panel.procLabel,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.procTB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.procTB,'Style','Edit');
set(hui.fig.ctrl_panel.procTB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.procTB,'Callback',@procTB_callback);

hui.fig.ctrl_panel.createPB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.createPB,'Style','PushButton');
set(hui.fig.ctrl_panel.createPB,'String','Create');
set(hui.fig.ctrl_panel.createPB,'Callback',@createPB_callback);
set(hui.fig.ctrl_panel.createPB,'Enable','off');

hui.fig.ctrl_panel.deletePB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.deletePB,'Style','PushButton');
set(hui.fig.ctrl_panel.deletePB,'String','Delete');
set(hui.fig.ctrl_panel.deletePB,'Callback',@deletePB_callback);

hui.fig.ctrl_panel.frmTable.ui = []; % Parent is a table container
hui.fig.ctrl_panel.frmTable.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.fig.ctrl_panel.frmTable.height_margin = NaN*zeros(30,30);
hui.fig.ctrl_panel.frmTable.false_width = NaN*zeros(30,30);
hui.fig.ctrl_panel.frmTable.false_height = NaN*zeros(30,30);
row = 1; col = 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.startLabel;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
col = col + 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.startTB;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
col = 1; row = row + 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.stopLabel;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
col = col + 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.stopTB;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
col = 1; row = row + 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.procLabel;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
col = col + 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.procTB;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
col = 1; row = row + 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.createPB;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
col = col + 1;
hui.fig.ctrl_panel.frmTable.handles{row,col}   = hui.fig.ctrl_panel.deletePB;
hui.fig.ctrl_panel.frmTable.width(row,col)     = inf;
hui.fig.ctrl_panel.frmTable.height(row,col)    = 25;
hui.fig.ctrl_panel.frmTable.false_width(row,col) = 5;
clear row col

hui.fig.ctrl_panel.autoGenSpacer = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.autoGenSpacer,'Style','Text');
set(hui.fig.ctrl_panel.autoGenSpacer,'String','__________________');
set(hui.fig.ctrl_panel.autoGenSpacer,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.headingLabel = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.headingLabel,'Style','Text');
set(hui.fig.ctrl_panel.headingLabel,'String','Max Heading (deg)');
set(hui.fig.ctrl_panel.headingLabel,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.headingTB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.headingTB,'Style','Edit');
set(hui.fig.ctrl_panel.headingTB,'String','5');
set(hui.fig.ctrl_panel.headingTB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.headingTB,'Enable','off');

hui.fig.ctrl_panel.lengthLabel = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.lengthLabel,'Style','Text');
set(hui.fig.ctrl_panel.lengthLabel,'String','Min Length (m)');
set(hui.fig.ctrl_panel.lengthLabel,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.lengthTB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.lengthTB,'Style','Edit');
set(hui.fig.ctrl_panel.lengthTB,'String','5000');
set(hui.fig.ctrl_panel.lengthTB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.lengthTB,'Enable','off');

hui.fig.ctrl_panel.maxLengthLabel = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.maxLengthLabel,'Style','Text');
set(hui.fig.ctrl_panel.maxLengthLabel,'String','Max Length (m)');
set(hui.fig.ctrl_panel.maxLengthLabel,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.maxLengthTB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.maxLengthTB,'Style','Edit');
set(hui.fig.ctrl_panel.maxLengthTB,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.autoTable.ui = []; % Parent is a table container
hui.fig.ctrl_panel.autoTable.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.fig.ctrl_panel.autoTable.height_margin = NaN*zeros(30,30);
hui.fig.ctrl_panel.autoTable.false_width = NaN*zeros(30,30);
hui.fig.ctrl_panel.autoTable.false_height = NaN*zeros(30,30);
row = 1; col = 1;
hui.fig.ctrl_panel.autoTable.handles{row,col}   = hui.fig.ctrl_panel.headingLabel;
hui.fig.ctrl_panel.autoTable.width(row,col)     = inf;
hui.fig.ctrl_panel.autoTable.height(row,col)    = 25;
hui.fig.ctrl_panel.autoTable.false_width(row,col) = 5;
col = col + 1;
hui.fig.ctrl_panel.autoTable.handles{row,col}   = hui.fig.ctrl_panel.headingTB;
hui.fig.ctrl_panel.autoTable.width(row,col)     = inf;
hui.fig.ctrl_panel.autoTable.height(row,col)    = 25;
hui.fig.ctrl_panel.autoTable.false_width(row,col) = 5;
col = 1; row = row + 1;
hui.fig.ctrl_panel.autoTable.handles{row,col}   = hui.fig.ctrl_panel.lengthLabel;
hui.fig.ctrl_panel.autoTable.width(row,col)     = inf;
hui.fig.ctrl_panel.autoTable.height(row,col)    = 25;
hui.fig.ctrl_panel.autoTable.false_width(row,col) = 5;
col = col + 1;
hui.fig.ctrl_panel.autoTable.handles{row,col}   = hui.fig.ctrl_panel.lengthTB;
hui.fig.ctrl_panel.autoTable.width(row,col)     = inf;
hui.fig.ctrl_panel.autoTable.height(row,col)    = 25;
hui.fig.ctrl_panel.autoTable.false_width(row,col) = 5;
col = 1; row = row + 1;
hui.fig.ctrl_panel.autoTable.handles{row,col}   = hui.fig.ctrl_panel.maxLengthLabel;
hui.fig.ctrl_panel.autoTable.width(row,col)     = inf;
hui.fig.ctrl_panel.autoTable.height(row,col)    = 25;
hui.fig.ctrl_panel.autoTable.false_width(row,col) = 5;
col = col + 1;
hui.fig.ctrl_panel.autoTable.handles{row,col}   = hui.fig.ctrl_panel.maxLengthTB;
hui.fig.ctrl_panel.autoTable.width(row,col)     = inf;
hui.fig.ctrl_panel.autoTable.height(row,col)    = 25;
hui.fig.ctrl_panel.autoTable.false_width(row,col) = 5;

hui.fig.ctrl_panel.autoGenPB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.autoGenPB,'Style','PushButton');
set(hui.fig.ctrl_panel.autoGenPB,'String','Autogenerate');
set(hui.fig.ctrl_panel.autoGenPB,'Callback',@autoGenPB_callback);

hui.fig.ctrl_panel.saveSpacer = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.saveSpacer,'Style','Text');
set(hui.fig.ctrl_panel.saveSpacer,'String','__________________');
set(hui.fig.ctrl_panel.saveSpacer,'HorizontalAlignment','Center');

hui.fig.ctrl_panel.savePB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.savePB,'Style','PushButton');
set(hui.fig.ctrl_panel.savePB,'String','Save Frames');
set(hui.fig.ctrl_panel.savePB,'Callback',@savePB_callback);

hui.fig.ctrl_panel.table.ui = hui.fig.ctrl_panel.handle;
hui.fig.ctrl_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.fig.ctrl_panel.table.height_margin = NaN*zeros(30,30);
hui.fig.ctrl_panel.table.false_width = NaN*zeros(30,30);
hui.fig.ctrl_panel.table.false_height = NaN*zeros(30,30);
hui.fig.ctrl_panel.table.offset = [0 10];
row = 1; col = 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.framesLB;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = inf;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.frmTable;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 100;
hui.fig.ctrl_panel.table.false_width(row,col) = 0;
hui.fig.ctrl_panel.table.width_margin(row,col) = 0;
hui.fig.ctrl_panel.table.height_margin(row,col) = 0;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.autoGenSpacer;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 20;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.autoTable;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 75;
hui.fig.ctrl_panel.table.false_width(row,col) = 0;
hui.fig.ctrl_panel.table.width_margin(row,col) = 0;
hui.fig.ctrl_panel.table.height_margin(row,col) = 0;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.autoGenPB;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 25;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.saveSpacer;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 20;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.savePB;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 25;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
clear row col
table_draw(hui.fig.ctrl_panel.table);

%% Read in records and frames files and populate GUI
% =====================================================================
% =====================================================================

GB = struct('param',param);

GB.records = records_file.records;
GB.records.along_track = geodetic_to_along_track(GB.records.lat,GB.records.lon);

if any(strcmpi(output_dir,'rds'))
  GB.default_frame_len = 50000;
elseif any(strcmpi(output_dir,'accum'))
  GB.default_frame_len = 20000;
elseif any(strcmpi(output_dir,{'kaband','kuband','snow'}))
  GB.default_frame_len = 5000;
end
if isfield(param.records,'frame_length') & ~isempty(param.records.frame_length)
  GB.default_frame_len = param.records.frame_length;
end
set(hui.fig.ctrl_panel.maxLengthTB,'String',sprintf('%i',GB.default_frame_len));

% Load frames file or create default frame list if it does not exist
if exist(ct_filename_support(GB.param,'','frames'),'file')
  frames_fn = ct_filename_support(GB.param,'','frames');
  frames_file = load(frames_fn,'frames');
  if find(frames_file.frames.frame_idxs > length(GB.records.lat))
    warning('Frames file %s\ncontains indices past the end of the records file (delete old file?)', frames_fn);
    keyboard
    frames_file.frames.frame_idxs = frames_file.frames.frame_idxs(frames_file.frames.frame_idxs <= length(GB.records.lat));
  end
  GB.frames = frames_file.frames;
else
  GB.frames.frame_idxs = 1;
  GB.frames.nyquist_zone = NaN;
  GB.frames.proc_mode = 0;
end
GB.cur_frame = 1;

if isempty(GB.param.records.geotiff_fn)
  GB.map_mode = false;
else
  GB.param.records.geotiff_fn = ct_filename_gis(param, GB.param.records.geotiff_fn);
  % Using Matlab's mapping mode with the geotiff provided
  GB.map_mode = true;
  % Using Matlab's mapping mode with the geotiff provided
  GB.map_mode = true;
  % Obtain the projection structure.
  GB.proj = geotiffinfo(GB.param.records.geotiff_fn);
  % Read the image
  [RGB, R, tmp] = geotiffread(GB.param.records.geotiff_fn);
  if size(RGB,3) == 3 && strcmp(class(RGB),'uint16') && max(RGB(:)) <= 255
    RGB = uint8(RGB);
  end
  if strcmpi(class(RGB),'int16') || strcmpi(class(RGB),'single')
    RGB = double(RGB);
  end
  R = R/1e3;
  
  [GB.X,GB.Y] = projfwd(GB.proj,double(GB.records.lat),double(GB.records.lon));
  GB.X = GB.X/1e3;
  GB.Y = GB.Y/1e3;
  
  mapshow(RGB, R);
end

% Update GUI with new frame list
update_frame_gui;

% Set callback properties
zoom on; zoom off;
set(hui.fig.handle,'Interruptible','off');
set(hui.fig.handle,'WindowButtonUpFcn',@gui_windowbuttonupfcn)
set(hui.fig.handle,'WindowButtonMotionFcn',[])
set(hui.fig.handle,'WindowKeyPressFcn',@gui_windowkeypressfcn);

end

%% framesLB_callback
function framesLB_callback(hObj,event)

global hui;
global GB;

GB.cur_frame = get(hui.fig.ctrl_panel.framesLB,'Value');
update_frame_gui;

end

%% procTB_callback
function procTB_callback(hObj,event)

global hui;
global GB;

GB.frames.proc_mode(GB.cur_frame) = str2double(get(hObj,'String'));
    
end

%% createPB_callback
function createPB_callback(hObj,event)

global hui;

end

%% deletePB_callback
function deletePB_callback(hObj,event)

global hui;
global GB;

if GB.cur_frame == 1
  GB.frames.frame_idxs = GB.frames.frame_idxs([1 3:end]);
  GB.frames.nyquist_zone = GB.frames.nyquist_zone([1 3:end]);
  GB.frames.proc_mode = GB.frames.proc_mode([1 3:end]);
else
  GB.frames.frame_idxs = GB.frames.frame_idxs([1:GB.cur_frame-1 GB.cur_frame+1:end]);
  GB.frames.nyquist_zone = GB.frames.nyquist_zone([1:GB.cur_frame-1 GB.cur_frame+1:end]);
  GB.frames.proc_mode = GB.frames.proc_mode([1:GB.cur_frame-1 GB.cur_frame+1:end]);
  GB.cur_frame = GB.cur_frame - 1;
end
update_frame_gui;

end

%% autoGenPB_callback
% Autogenerate
% =====================================================================
function autoGenPB_callback(hObj,event)

global hui; % hui: user interface handles
global GB; % Geobase: contains all the geographic info

cur_along_track = GB.records.along_track(GB.frames.frame_idxs(GB.cur_frame));
if GB.cur_frame < length(GB.frames.frame_idxs)
  next_idx = GB.frames.frame_idxs(GB.cur_frame+1);
else
  next_idx = length(GB.records.along_track);
end
max_frm_length = str2double(get(hui.fig.ctrl_panel.maxLengthTB,'String'));
for idx = GB.frames.frame_idxs(GB.cur_frame):next_idx-1
  if GB.records.along_track(idx)-cur_along_track > max_frm_length
    cur_along_track = GB.records.along_track(idx);
    
    % Create a new frame
    first_ind = GB.cur_frame+1;
    GB.frames.frame_idxs(first_ind+1:end+1) = GB.frames.frame_idxs(first_ind:end);
    GB.frames.frame_idxs(first_ind) = idx;
    GB.frames.nyquist_zone(first_ind+1:end+1) = GB.frames.nyquist_zone(first_ind:end);
    GB.frames.nyquist_zone(first_ind) = NaN;
    GB.frames.proc_mode(first_ind+1:end+1) = GB.frames.proc_mode(first_ind:end);
    GB.frames.proc_mode(first_ind) = str2double(get(hui.fig.ctrl_panel.procTB,'String'));

    GB.cur_frame = first_ind;
  end
end
update_frame_gui;

end

%% savePB_callback
% Save Frames
% =====================================================================
function savePB_callback(hObj,event)

global hui; % hui: user interface handles
global GB; % Geobase: contains all the geographic info

fn = ct_filename_support(GB.param,'','frames');
[out_path] = fileparts(fn);
if ~exist(out_path,'dir')
  fprintf('Making directory %s\n', out_path);
  mkdir(out_path);
end
save(fn,'-STRUCT','GB','frames');
fprintf('Frames saved in %s\n', fn);

end

%% find_closest_point
% Support function for finding closest point
% =====================================================================
function [min_dist min_ind] = find_closest_point(x,y)

global hui; % hui: user interface handles
global GB; % Geobase: contains all the geographic info

if GB.map_mode
  [min_dist min_ind] = min(sqrt((x - GB.X).^2 + (y - GB.Y).^2));
else
  [min_dist min_ind] = min(sqrt((x - GB.records.lon).^2 + (y - GB.records.lat).^2));
end

end

%% update_frame_gui
% Support function for updating frame listbox, map, and frame
% controls
% =====================================================================
function update_frame_gui

global hui; % hui: user interface handles
global GB; % Geobase: contains all the geographic info

% Set strings in GUI according to frames list
for idx = 1:length(GB.frames.frame_idxs)
  if idx == length(GB.frames.frame_idxs)
    frm = [GB.frames.frame_idxs(idx) length(GB.records.lat)];
  else
    frm = [GB.frames.frame_idxs(idx); GB.frames.frame_idxs(idx+1)-1];
  end
  frm_text{idx} = sprintf('Frm %2d: %8d-%8d %4.1f km', idx, frm(1), frm(2), diff(GB.records.along_track(frm(1:2))/1e3));
end
if GB.cur_frame == length(GB.frames.frame_idxs)
  frm = [GB.frames.frame_idxs(GB.cur_frame) length(GB.records.lat)];
else
  frm = [GB.frames.frame_idxs(GB.cur_frame); GB.frames.frame_idxs(GB.cur_frame+1)-1];
end
set(hui.fig.ctrl_panel.framesLB,'String',frm_text);
set(hui.fig.ctrl_panel.startTB,'String',sprintf('%d',frm(1)));
set(hui.fig.ctrl_panel.stopTB,'String',sprintf('%d',frm(2)));
set(hui.fig.ctrl_panel.procTB,'String',sprintf('%d',GB.frames.proc_mode(GB.cur_frame)));

% =====================================================================
% Plot map
axes(hui.fig.map_axes.handle);

% Delete the old handles
if isfield(hui.fig.map_axes,'frm_starts')
  for idx = 1:length(hui.fig.map_axes.frm_starts)
    delete(hui.fig.map_axes.frm_starts(idx));
    delete(hui.fig.map_axes.frm_stops(idx));
  end
  hui.fig.map_axes.frm_starts = [];
  hui.fig.map_axes.frm_stops = [];
end
if isfield(hui.fig.map_axes,'line')
  delete(hui.fig.map_axes.line);
end
if isfield(hui.fig.map_axes,'cur_frame')
  delete(hui.fig.map_axes.cur_frame);
end

if GB.map_mode
  hold on;
  % Plot the whole flight line
  hui.fig.map_axes.line = plot(GB.X(1:50:end),GB.Y(1:50:end));
  xlabel('X (km)');
  ylabel('Y (km)');
  title(GB.proj.Projection);
  
  % Plot each frame with special for the current frame
  for idx = 1:length(GB.frames.frame_idxs)
    if idx == length(GB.frames.frame_idxs)
      frm = [GB.frames.frame_idxs(idx) length(GB.records.lat)];
    else
      frm = [GB.frames.frame_idxs(idx); GB.frames.frame_idxs(idx+1)-1];
    end
    hui.fig.map_axes.frm_starts(idx) = plot(GB.X(frm(1)),GB.Y(frm(1)),'ro');
    hui.fig.map_axes.frm_stops(idx) = plot(GB.X(frm(2)),GB.Y(frm(2)),'rx');
    if idx == GB.cur_frame
      hui.fig.map_axes.cur_frame = plot(GB.X(frm(1):50:frm(2)),GB.Y(frm(1):50:frm(2)),'r-','LineWidth',2);
    end
  end
  hold off;
  
else
  hold on;
  % Plot the whole flight line
  hui.fig.map_axes.line = plot(GB.records.lon(1:50:end),GB.records.lat(1:50:end));
  xlabel('Longitude (E)');
  ylabel('Latitude (N)');
  
  % Plot each frame with special for the current frame
  for idx = 1:length(GB.frames.frame_idxs)
    if idx == length(GB.frames.frame_idxs)
      frm = [GB.frames.frame_idxs(idx) length(GB.records.lat)];
    else
      frm = [GB.frames.frame_idxs(idx); GB.frames.frame_idxs(idx+1)-1];
    end
    hui.fig.map_axes.frm_starts(idx) = plot(GB.records.lon(frm(1)),GB.records.lat(frm(1)),'ro');
    hui.fig.map_axes.frm_stops(idx) = plot(GB.records.lon(frm(2)),GB.records.lat(frm(2)),'rx');
    if idx == GB.cur_frame
      hui.fig.map_axes.cur_frame = plot(GB.records.lon(frm(1):50:frm(2)),GB.records.lat(frm(1):50:frm(2)),'r-');
    end
  end
  hold off;
end
grid on;

set(hui.fig.ctrl_panel.framesLB,'Value',GB.cur_frame);

end

%% gui_windowbuttonupfcn
% WindowButtonUpFcn call back function
% =====================================================================
function gui_windowbuttonupfcn(src,event)

global hui; % hui: user interface handles
global GB; % Geobase: contains all the geographic info

[x,y,but] = get_mouse_info(hui.fig.handle,hui.fig.map_axes.handle);
%fprintf('1: x = %.2f, y = %.2f, but = %d\n', x, y, but);

xlimits = xlim;
ylimits = ylim;
if x < xlimits(1) || x > xlimits(2) || y < ylimits(1) || y > ylimits(2)
  return;
end

if but == 1
  % ===================================================================
  % Left mouse button: add a frame
  
  % Find the closest point
  [min_dist min_ind] = find_closest_point(x,y);

  % Create a new frame
  first_ind = find(GB.frames.frame_idxs >= min_ind,1);
  if isempty(first_ind)
    first_ind = length(GB.frames.frame_idxs) + 1;
  elseif min_ind == GB.frames.frame_idxs(first_ind)
    return;
  end
  GB.frames.frame_idxs(first_ind+1:end+1) = GB.frames.frame_idxs(first_ind:end);
  GB.frames.frame_idxs(first_ind) = min_ind;
  GB.frames.nyquist_zone(first_ind+1:end+1) = GB.frames.nyquist_zone(first_ind:end);
  GB.frames.nyquist_zone(first_ind) = NaN;
  GB.frames.proc_mode(first_ind+1:end+1) = GB.frames.proc_mode(first_ind:end);
  GB.frames.proc_mode(first_ind) = str2double(get(hui.fig.ctrl_panel.procTB,'String'));

  GB.cur_frame = first_ind;
  update_frame_gui;
else
  % ===================================================================
  % Right mouse button: select a frame
  
  % Find the closest point
  [min_dist min_ind] = find_closest_point(x,y);
  
  % Determine which frame it is in
  GB.cur_frame = find(min_ind >= GB.frames.frame_idxs,1,'last');
  update_frame_gui;
end

end

%% WindowKeyPressFcn call back function
% =====================================================================
function gui_windowkeypressfcn(src,event)

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)

% Check to make sure that a key was pressed and not
% just a modifier (e.g. shift, ctrl, alt)
if isempty(event.Key) || strcmpi(event.Key,'shift') || strcmpi(event.Key,'alt') || strcmpi(event.Key,'ctrl')
  return;
end

global hui; % hui: user interface handles
global GB; % Geobase: contains all the geographic info

[x,y,but] = get_mouse_info(hui.fig.handle,hui.fig.map_axes.handle);
%fprintf('1: x = %.2f, y = %.2f, but = %d\n', x, y, but);

xlimits = xlim;
ylimits = ylim;
if x < xlimits(1) || x > xlimits(2) || y < ylimits(1) || y > ylimits(2)
  return;
end

if 0 % This code for debugging
  if ischar(event.Key)
    fprintf('x = %.2f, y = %.2f, key = %s %d\n', x, y, event.Key, event.Character);
  else
    fprintf('x = %.2f, y = %.2f, key = %d %d\n', x, y, event.Key, event.Character);
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
switch event.Key
  case 'backspace' % Backspace
    deletePB_callback;
  case 'delete' % Delete
    if GB.cur_frame < length(GB.frames.frame_idxs)
      GB.cur_frame = GB.cur_frame + 1;
    end
    deletePB_callback;
  case 'downarrow' % Down-arrow
    if GB.cur_frame < length(GB.frames.frame_idxs)
      GB.cur_frame = GB.cur_frame + 1;
    end
    set(hui.fig.ctrl_panel.framesLB,'Value',GB.cur_frame);
    update_frame_gui;
  case 'uparrow' % Up-arrow
    if GB.cur_frame > 1
      GB.cur_frame = GB.cur_frame - 1;
    end
    set(hui.fig.ctrl_panel.framesLB,'Value',GB.cur_frame);
    update_frame_gui;
  case 'z'
    if any(strcmp('shift',event.Modifier))
      % reset zoom/map axis
      axis tight
    else
      % Turn zoom mode on
      zoom on
    end
  case 'f1' % Help information
    fprintf('--------------------create_frames help--------------------\n\n');
    fprintf('F1: Help\n');
    fprintf('left click: Create a new frame break at click\n');
    fprintf('right click: Select the closest frame\n');
    fprintf('Backspace: Delete previous frame break\n');
    fprintf('Delete: Delete next frame break\n');
    fprintf('Up: Select the previous frame\n');
    fprintf('Down: Select the next frame\n');
    fprintf('z: Switch to Matlab zoom mode (click on magnifying glass to switch back)\n');
    fprintf('Z: Reset the zoom/map axis\n');
end

end
