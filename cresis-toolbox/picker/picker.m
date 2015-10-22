function picker(source_data, layer_data, geotiff_fns, param)
% picker(source_data, layer_data, geotiff_fns, param)
%
% source_data = cell column vector of paths to source data directories
%   First path listed is master path. Files in each source directory
%   must have the same name otherwise they will not link properly.
% layer_data = path to layer outputs directory (name of layer files
%   must match master source filename)
% geotiff_fns = cell column vector of geotiff's that can be loaded into map view
%
% Figure 3 and 4 Hotkeys (figure 4 does not have tools enabled)
%  SPACEBAR: toggle layers on/off
%  SHIFT-SPACEBAR: toggle between layer data and quality
%  zZ: zoom to full view
%  bemdiqsM: Choose corresponding tool
%  123: Choose surface/bottom/both
%  S: save layer
%  f: Flip display left/right
%  Q: Cycle through quality levels
%  u: undo last tool command
%  np: load (n)ext or (p)revious frame
%  ARROW-KEYS: Move axis around
%  CTRL-ARROW-KEYS: For up/down moves the current layer up/down
%  SHIFT-ARROW-KEYS: For left/right arrow keys goes to prev/next frame, but
%    keeps axes
%    For up/down arrow keys moves through source list
%  ESCAPE: escape last click
%  HOME/END/PAGE-UP/PAGE-DOWN: increment/decrement caxis bounds
%  Left-click: apply tool
%  Right-click: move cursor
%
% Figure 2 Hotkeys
%  Left-click: Select closest frame
%  Double-click: Load currently selected frame
%  Right-click: Move cursors in pick and view
%
% Example: See run_picker.m
%
% Author: Anthony Hoch, John Paden, Aric Beaver

% Figure 1: Control Window
% Figure 2: Map Window
% Figure 3: Pick Window
% Figure 4: View Window

% Figure 1:
% - Select map to display (default is none)
% - Select data segment (default is first)
% - Select source folder (default is first)
% - Image domain processing (Raw, Average, Detrend, default is Raw)
% - Load to pick
% - Load to view
% - Y-limits to use when loading new figures
%
% - Active picking tool (browse, enter pnt, max pnt, delete region,
%   quality region, interp region, snake region, max region, default is browse)
% - Active layer (surface, bottom, both, default surface)
% - Active quality (good, moderate, derived, default good)
% - Save layers of current segment loaded into picking window
% - Tool parameters

% ---------------------------------------------------------------------
% Layer Files:
%
% Two Layers, 1: surface, 2: bottom
% Layer Data Value 1: Manually picked points
%  inf: no pick
%  finite: propagation time to target
% Layer Data Value 2: Automatically generated points
%  inf: no pick
%  finite: propagation time to target
% Quality control level
%  0: good
%  1: moderate
%  2: derived

% ===================================================================
% Setup global variables
% ===================================================================

% hui: contains all the figure and user-interface/widget handles
% gCtrl: contains all the control variables
global gCtrl;
global hui;
gCtrl = [];
hui = [];

gCtrl.tic = tic;
gCtrl.source.paths = source_data;
gCtrl.layer_data = layer_data;
gCtrl.geotiff.fns = geotiff_fns;
gCtrl.geotiff.fns(2:end+1) = gCtrl.geotiff.fns(1:end);
gCtrl.geotiff.fns{1} = 'None';

gCtrl.radar_name = param.radar_name;
gCtrl.season_name = param.season_name;

if ~exist('param','var')
  param.fast_load.en = false;
end
if ~isfield(param,'fast_load')
  param.fast_load.en = false;
end
if ~isfield(param.fast_load,'en')
  param.fast_load.en = false;
end
if ~isfield(param,'landmarks')
  param.landmarks = '';
end
gCtrl.source.landmarks.fn = param.landmarks;
if ~isfield(param,'day_seg_mask')
  param.day_seg_mask = '';
end

gCtrl.geotiff.cur = 1;

hui.fig.handle = 1;
try; delete(hui.fig.handle); end;
hui.mapfig.handle = 2;
try; delete(hui.mapfig.handle); end;
hui.pickfig.handle = 3;
try; delete(hui.pickfig.handle); end;
hui.viewfig.handle = 4;
try; delete(hui.viewfig.handle); end;
hui.landfig.handle = 7;
try; delete(hui.landfig.handle); end;

% Initialize current source path
gCtrl.source.cur_src = 1;

% Initialize current map selection, pick, and view frame
gCtrl.source.cur_sel = 1;
gCtrl.source.cur_pick = 1;
gCtrl.source.cur_view = 1;

% Initialize cursor indices
gCtrl.pick.cur_idx = 1;
gCtrl.view.cur_idx = 1;

% Initialize tool settings
gCtrl.undo_stack_ptr = 0;
gCtrl.undo_stack = [];
gCtrl.tool.layer_multiple = 1;
gCtrl.tool.selection.x = NaN;
gCtrl.tool.layer_switch = false;
gCtrl.tool.layer_visible_pick = [1 1 1 1 1]; % [layer_or_quality master-visible 1-visible 2-visible 3-visible] as [1or2 boolean boolean boolean boolean]
gCtrl.tool.layer_visible_view = [1 1 1 1 1]; % [layer_or_quality master-visible 1-visible 2-visible 3-visible] as [1or2 boolean boolean boolean boolean]

% ===================================================================
% ===================================================================
% Prepare source directory list
%  - Load in layer files
% ===================================================================
% ===================================================================

% Remove trailing file separation character from paths if it exists
for path_idx = 1:length(gCtrl.source.paths)
  if gCtrl.source.paths{path_idx}(end) == filesep
    gCtrl.source.paths{path_idx} = gCtrl.source.paths{path_idx}(1:end-1);
  end
end

%% Search for all layer files which match param.day_seg_mask
source_dirs = get_filenames(gCtrl.layer_data,param.day_seg_mask,'','',struct('type','d'));
if isempty(source_dirs)
  error('No layer directories in %s match param.day_seg_mask = %s', gCtrl.layer_data, param.day_seg_mask);
end
gCtrl.source.fns = get_filenames(source_dirs{1},'Data','','.mat');
for source_idx = 2:length(source_dirs)
  gCtrl.source.fns = cat(1,gCtrl.source.fns, ...
    get_filenames(source_dirs{source_idx},'Data','','.mat'));
end

if isempty(gCtrl.source.fns)
  error('No layer files found in %s', gCtrl.layer_data);
end
% Sort the filenames by frame ID rather than full path
clear names;
for idx = 1:length(gCtrl.source.fns)
  % Grab the filename from the full path
  [tmp names{idx}] = fileparts(gCtrl.source.fns{idx});
  % Grab just the frame ID
  names{idx} = names{idx}(6:end);
end
[names new_ordering] = sort(names);
gCtrl.source.fns = gCtrl.source.fns(new_ordering);

% If a source data mask is supplied, only keep frame IDs which
% satisfy the mask
if isfield(param,'source_data_mask')
  good_mask = logical(zeros(size(names)));
  for idx = 1:length(param.source_data_mask)
    good_idxs = strmatch(param.source_data_mask{idx},names);
    good_mask(good_idxs) = 1;
  end
  gCtrl.source.fns = gCtrl.source.fns(good_mask);
  names = names(good_mask);
end
clear tmp;
if isempty(names)
  error('No layer files found which match source_data_mask');
end

tmp_file = ct_filename_tmp(param,param.fast_load.tmp_file,'picker','picker_fast_load.mat');

if param.fast_load.en && exist(tmp_file,'file') ...
    && ~param.fast_load.recreate
  fprintf('Fast load\n');
  warning off;
  for file_idx = 1:length(gCtrl.source.fns)
    fn = gCtrl.source.fns{file_idx};
    gCtrl.source.layers{file_idx} = load(fn,'layerData');
    gCtrl.source.geo{file_idx} = load(fn,'Latitude','Longitude','GPS_time','Elevation');
  end
  warning on;
  load(tmp_file);
  gCtrl.source.frm_id = fast_load.frm_id;
  gCtrl.source.src_fns = fast_load.src_fns;
  gCtrl.source.src_disp = fast_load.src_disp;
  
  fprintf('  Done (%.1f sec)\n', toc(gCtrl.tic));
else
  fprintf('Loading layer and GPS from the layer files\n');
  warning off;
  for file_idx = 1:length(gCtrl.source.fns)
    fn = gCtrl.source.fns{file_idx};
    gCtrl.source.layers{file_idx} = load(fn,'layerData');
    gCtrl.source.geo{file_idx} = load(fn,'Latitude','Longitude','GPS_time','Elevation');
    [fpath name] = fileparts(fn);
    gCtrl.source.frm_id(file_idx,:) = name(end-14:end);
    gCtrl.source.src_fns{file_idx} = {};
    gCtrl.source.src_disp{file_idx} = {};
  end
  warning on;
  fprintf('  Done (%.1f sec)\n', toc(gCtrl.tic));
  
  fprintf('Cataloging source directories (could take several minutes)\n');
  for path_idx = 1:length(gCtrl.source.paths)
    % We now have a list of layer files in gCtrl.source.fns
    % Our next task is to find all the echogram files associated with
    % each layer file. We do this by searching through every path in
    % gCtrl.source.paths and matching the frame ID (YYYYMMDD_segS_FF)
    fpath = gCtrl.source.paths{path_idx};
    if exist(fpath,'dir')
      
      %% Search for all echogram files which match param.day_seg_mask
      source_dirs = get_filenames(fpath,param.day_seg_mask,'','',struct('type','d'));
      if isempty(source_dirs)
        fns_list = {};
      else
        fns_list = get_filenames(source_dirs{1},'Data','','.mat');
        for source_idx = 2:length(source_dirs)
          fns_list = cat(1,fns_list, ...
            get_filenames(source_dirs{source_idx},'Data','','.mat'));
        end
      end
      
      % Get directory name of this echogram (this will be used in the
      % listmenu GUI to indicate echograms from this directory)
      [tmp dir_name] = fileparts(fpath);
      if length(dir_name) > 6 && strcmp(dir_name(1:6),'CSARP_')
        dir_name = dir_name(7:end);
      end
      
      % Pull out just the file names to make the strfind in the
      % next loop a little faster
      fns_name = {};
      for file_idx = 1:length(fns_list)
        [tmp fns_name{file_idx}] = fileparts(fns_list{file_idx});
        if fns_name{file_idx}(6) == 'i'
          fns_name{file_idx} = fns_name{file_idx}(13:end);
        else
          fns_name{file_idx} = fns_name{file_idx}(6:end);
        end
      end
      
      % For each echogram file, search through the list of layer files
      % for matching frame IDs
      for echo_idx = 1:length(fns_name)
        % Check if current file name is before the first layer item
        if frame_id_comp(fns_name{echo_idx},names{1}) < 0
          continue;
        end
        % Check if current file name is after the last layer item
        if frame_id_comp(fns_name{echo_idx},names{end}) > 0
          continue;
        end
        % Binary search through sorted layer list
        %  10 --> 5/6, 11 --> 6
        bot_idx = 1;
        top_idx = length(gCtrl.source.fns);
        comp_result = ~0;
        while comp_result
          file_idx = floor((bot_idx+top_idx)/2);
          comp_result = frame_id_comp(fns_name{echo_idx},names{file_idx});
          if bot_idx >= top_idx
            break;
          end
          if comp_result < 0
            top_idx = file_idx - 1;
          elseif comp_result > 0
            bot_idx = file_idx + 1;
          end
        end
        if comp_result == 0
          % A frame ID matched, so add this echogram to this layer files
          % list.
          
          % For backwards compatibility, add the geometry if the layer
          % file did not have it
          if ~isfield(gCtrl.source.geo{file_idx},'Latitude')
            gCtrl.source.geo{file_idx} ...
              = load(fns_list{echo_idx},'Latitude','Longitude','GPS_time','Elevation');
          end
          
          % Now add the full path to the echogram
          gCtrl.source.src_fns{file_idx} ...
            = cat(1,gCtrl.source.src_fns{file_idx}, fns_list{echo_idx});
          
          % Now add the GUI name for the echogram
          [fpath name] = fileparts(fns_list{echo_idx});
          if name(6:8) == 'img'
            img_name = ['/' name(6:11)];
          else
            img_name = [];
          end
          gCtrl.source.src_disp{file_idx}{end+1} = [dir_name img_name];
        end
      end
    end
    
  end
  
  fprintf('  Done (%.1f sec)\n', toc(gCtrl.tic));
  if param.fast_load.en
    fast_load.layers = gCtrl.source.layers;
    fast_load.geo = gCtrl.source.geo;
    fast_load.frm_id = gCtrl.source.frm_id;
    fast_load.src_fns = gCtrl.source.src_fns;
    fast_load.src_disp = gCtrl.source.src_disp;
    tmp_path = fileparts(tmp_file);
    if ~exist(tmp_path,'dir')
      mkdir(tmp_path);
    end
    save('-v6',tmp_file,'fast_load');
  end
end

% Create global modified flag array, we will use ' ' to represent that the
% frame has not been modified since last saved and '*' to represent that it
% has been modified
gCtrl.source.modified = ' '*ones(length(gCtrl.source.fns),1);

% Initialize cursor positions
gCtrl.pick.cursor = gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time(gCtrl.pick.cur_idx);
gCtrl.view.cursor = gCtrl.source.geo{gCtrl.source.cur_view}.GPS_time(gCtrl.view.cur_idx);

% ===================================================================
% ===================================================================
% Create user interface figure
% ===================================================================
% ===================================================================
figure(hui.fig.handle);
set(hui.fig.handle,'Position',[10 200 240 640]);
set(hui.fig.handle,'DockControls','off')
set(hui.fig.handle,'NumberTitle','off');
set(hui.fig.handle,'Name','1: Control');
set(hui.fig.handle,'ToolBar','none');
set(hui.fig.handle,'MenuBar','none');
set(hui.fig.handle,'CloseRequestFcn',@windowClose_callback);

% ===================================================================
% Create user interface figure table contents
% ===================================================================

% -------------------------------------------------------------------
% Control Panel
hui.fig.ctrl_panel.handle = uipanel('Parent',hui.fig.handle);
set(hui.fig.ctrl_panel.handle,'Title','Pick Controls');
set(hui.fig.ctrl_panel.handle,'TitlePosition','CenterTop');
set(hui.fig.ctrl_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(hui.fig.ctrl_panel.handle,'ShadowColor',[0.6 0.6 0.6]);

% ===================================================================
% Setup user interface figure table
% ===================================================================
hui.fig.table.ui = hui.fig.handle;
row = 1; col = 1;
hui.fig.table.handles{row,col}   = hui.fig.ctrl_panel.handle;
hui.fig.table.width(row,col)     = inf;
hui.fig.table.height(row,col)    = inf;
hui.fig.table.false_height(row,col) = 0;
clear row col
table_draw(hui.fig.table);

% ===================================================================
% Create ctrl_panel table contents
% ===================================================================
% -------------------------------------------------------------------
% Map Popup Menu
hui.fig.ctrl_panel.mapsLB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
menuString = {};
for idx_menu = 1:length(gCtrl.geotiff.fns)
  [path name] = fileparts(gCtrl.geotiff.fns{idx_menu});
  menuString{end+1} = name;
end
set(hui.fig.ctrl_panel.mapsLB,'String',menuString);
clear menuString;
set(hui.fig.ctrl_panel.mapsLB,'Value',gCtrl.geotiff.cur);
set(hui.fig.ctrl_panel.mapsLB,'Style','popupmenu');
set(hui.fig.ctrl_panel.mapsLB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.mapsLB,'FontName','fixed');
set(hui.fig.ctrl_panel.mapsLB,'Callback',@mapsLB_callback);

% -------------------------------------------------------------------
% Frames Listbox
hui.fig.ctrl_panel.framesLB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
menuString = {};
for idx_menu = 1:length(gCtrl.source.fns)
  menuString{end+1} = [gCtrl.source.frm_id(idx_menu,:) ' ' gCtrl.source.modified(idx_menu)];
end
set(hui.fig.ctrl_panel.framesLB,'String',menuString);
set(hui.fig.ctrl_panel.framesLB,'Value',gCtrl.source.cur_sel);
set(hui.fig.ctrl_panel.framesLB,'Style','listbox');
set(hui.fig.ctrl_panel.framesLB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.framesLB,'FontName','fixed');
set(hui.fig.ctrl_panel.framesLB,'Callback',@framesLB_callback);

% -------------------------------------------------------------------
% Source data Listbox
hui.fig.ctrl_panel.sourceLB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
menuString = gCtrl.source.src_disp{gCtrl.source.cur_sel};
set(hui.fig.ctrl_panel.sourceLB,'String',menuString);
set(hui.fig.ctrl_panel.sourceLB,'Value',gCtrl.source.cur_src);
set(hui.fig.ctrl_panel.sourceLB,'Style','listbox');
set(hui.fig.ctrl_panel.sourceLB,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.sourceLB,'FontName','fixed');
set(hui.fig.ctrl_panel.sourceLB,'Callback',@sourceLB_callback);

% -------------------------------------------------------------------
% Display Mode and Colormap Popup Menus
hui.fig.ctrl_panel.display_modePM = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.display_modePM,'Style','PopupMenu');
set(hui.fig.ctrl_panel.display_modePM,'String',{'Raw','Averaged','Detrend','Local Detrend','New_Detrend'});
set(hui.fig.ctrl_panel.display_modePM,'Callback',@display_modePM_callback);

hui.fig.ctrl_panel.colormapPM = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.colormapPM,'Style','PopupMenu');
set(hui.fig.ctrl_panel.colormapPM,'String',{'gray','jet','copper'});
set(hui.fig.ctrl_panel.colormapPM,'Callback',@colormapPM_callback);

% -------------------------------------------------------------------
% Load to View/Pick Pushbuttons and Y-limits
hui.fig.ctrl_panel.load_pickPB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.load_pickPB,'Style','PushButton');
set(hui.fig.ctrl_panel.load_pickPB,'String','Load Pick');
set(hui.fig.ctrl_panel.load_pickPB,'Callback',@load_pickPB_callback);

hui.fig.ctrl_panel.load_viewPB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.load_viewPB,'Style','PushButton');
set(hui.fig.ctrl_panel.load_viewPB,'String','Load View');
set(hui.fig.ctrl_panel.load_viewPB,'Callback',@load_viewPB_callback);

hui.fig.ctrl_panel.ylim_label= uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.ylim_label,'Style','text');
set(hui.fig.ctrl_panel.ylim_label,'String','Y-Limits (us)');

hui.fig.ctrl_panel.ylim_TE = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.ylim_TE,'Style','edit');
set(hui.fig.ctrl_panel.ylim_TE,'String','');
set(hui.fig.ctrl_panel.ylim_TE,'Callback',@ylim_TE_callback);

hui.fig.ctrl_panel.toolPM = uicontrol('Parent',hui.fig.ctrl_panel.handle);
menuString = {};
menuString{1} = '(b)rowse';
menuString{2} = '(e)nter pnt';
menuString{3} = '(m)ax pnt';
menuString{4} = '(d)elete';
menuString{5} = '(q)uality';
menuString{6} = '(i)nterp (w)';
menuString{7} = '(s)nake';
menuString{8} = '(M)ax';
menuString{9} = '(l)eading edge';
menuString{10} = 'landmar(k)';
menuString{11} = '(c)onvert layers';
set(hui.fig.ctrl_panel.toolPM,'String',menuString);
set(hui.fig.ctrl_panel.toolPM,'Value',1);
clear menuString;
set(hui.fig.ctrl_panel.toolPM,'Style','popupmenu');
set(hui.fig.ctrl_panel.toolPM,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.toolPM,'FontName','fixed');
set(hui.fig.ctrl_panel.toolPM,'Callback',@toolPM_callback);

hui.fig.ctrl_panel.layerPM = uicontrol('Parent',hui.fig.ctrl_panel.handle);
menuString = {};
menuString{1} = '(1) surface';
menuString{2} = '(2) bottom';
menuString{3} = '(3) both';
set(hui.fig.ctrl_panel.layerPM,'String',menuString);
set(hui.fig.ctrl_panel.layerPM,'Value',1);
clear menuString;
set(hui.fig.ctrl_panel.layerPM,'Style','popupmenu');
set(hui.fig.ctrl_panel.layerPM,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.layerPM,'FontName','fixed');

hui.fig.ctrl_panel.qualityPM = uicontrol('Parent',hui.fig.ctrl_panel.handle);
menuString = {};
menuString{1} = 'good';
menuString{2} = 'moderate';
menuString{3} = 'derived';
set(hui.fig.ctrl_panel.qualityPM,'String',menuString);
set(hui.fig.ctrl_panel.qualityPM,'Value',1);
clear menuString;
set(hui.fig.ctrl_panel.qualityPM,'Style','popupmenu');
set(hui.fig.ctrl_panel.qualityPM,'HorizontalAlignment','Center');
set(hui.fig.ctrl_panel.qualityPM,'FontName','fixed');

hui.fig.ctrl_panel.savePB = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.savePB,'Style','PushButton');
set(hui.fig.ctrl_panel.savePB,'String','Save Layer');
set(hui.fig.ctrl_panel.savePB,'Callback',@picker_save_layers);

hui.fig.ctrl_panel.tool_param1_TE = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.tool_param1_TE,'Style','edit');
set(hui.fig.ctrl_panel.tool_param1_TE,'String','');

hui.fig.ctrl_panel.tool_param2_TE = uicontrol('Parent',hui.fig.ctrl_panel.handle);
set(hui.fig.ctrl_panel.tool_param2_TE,'Style','edit');
set(hui.fig.ctrl_panel.tool_param2_TE,'String','');

% -------------------------------------------------------------------
% Load to View/Pick, Y-lim, Tools Table
hui.fig.ctrl_panel.butTable.ui = []; % Parent is a table container
hui.fig.ctrl_panel.butTable.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.fig.ctrl_panel.butTable.height_margin = NaN*zeros(30,30);
hui.fig.ctrl_panel.butTable.false_width = NaN*zeros(30,30);
hui.fig.ctrl_panel.butTable.false_height = NaN*zeros(30,30);
row = 1; col = 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.display_modePM;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = col + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.colormapPM;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = 1; row = row + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.load_pickPB;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = col + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.load_viewPB;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = 1; row = row + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.ylim_label;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height_margin(row,col) = 7;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = col + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.ylim_TE;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = 1; row = row + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.toolPM;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = col + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.layerPM;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = 1; row = row + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.qualityPM;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = col + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.savePB;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = 1; row = row + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.tool_param1_TE;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
col = col + 1;
hui.fig.ctrl_panel.butTable.handles{row,col}   = hui.fig.ctrl_panel.tool_param2_TE;
hui.fig.ctrl_panel.butTable.width(row,col)     = inf;
hui.fig.ctrl_panel.butTable.height(row,col)    = 25;
hui.fig.ctrl_panel.butTable.false_width(row,col) = 3;
clear row col

% ===================================================================
% Setup ctrl_panel table
% ===================================================================
hui.fig.ctrl_panel.table.ui = hui.fig.ctrl_panel.handle;
hui.fig.ctrl_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.fig.ctrl_panel.table.height_margin = NaN*zeros(30,30);
hui.fig.ctrl_panel.table.false_width = NaN*zeros(30,30);
hui.fig.ctrl_panel.table.false_height = NaN*zeros(30,30);
hui.fig.ctrl_panel.table.offset = [0 10];
row = 1; col = 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.mapsLB;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 25;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.framesLB;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = inf;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.sourceLB;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 80;
hui.fig.ctrl_panel.table.false_height(row,col) = 5;
hui.fig.ctrl_panel.table.false_width(row,col) = 5;
row = row + 1;
hui.fig.ctrl_panel.table.handles{row,col}   = hui.fig.ctrl_panel.butTable;
hui.fig.ctrl_panel.table.width(row,col)     = inf;
hui.fig.ctrl_panel.table.height(row,col)    = 160;
hui.fig.ctrl_panel.table.false_height(row,col) = 5;
hui.fig.ctrl_panel.table.false_width(row,col) = 0;
hui.fig.ctrl_panel.table.width_margin(row,col) = 0;
hui.fig.ctrl_panel.table.height_margin(row,col) = 0;
clear row col
table_draw(hui.fig.ctrl_panel.table);

% ===================================================================
% ===================================================================
% Create map figure
% ===================================================================
% ===================================================================
fprintf('Creating map figure (may take a minute)\n');
figure(hui.mapfig.handle);
set(hui.mapfig.handle,'NumberTitle','off');
set(hui.mapfig.handle,'Name','2: Map');
set(hui.mapfig.handle,'CloseRequestFcn',@windowClose_callback);
picker_map(1);
fprintf('  Done (%.1f sec)\n', toc(gCtrl.tic));

% ===================================================================
% ===================================================================
% Create landmark figure
% ===================================================================
% ===================================================================

fprintf('Creating landmark figure \n');

figure(hui.landfig.handle);
% [10 200 240 640]
set(hui.landfig.handle,'Position',[280 200 240 500]);
set(hui.landfig.handle,'DockControls','off')
set(hui.landfig.handle,'NumberTitle','off');
set(hui.landfig.handle,'Name','7: Landmark');
set(hui.landfig.handle,'ToolBar','none');
set(hui.landfig.handle,'MenuBar','none');
set(hui.landfig.handle,'CloseRequestFcn',@landmark_CloseRequestFcn)

% ===================================================================
% Create user interface figure table contents
% ===================================================================

% -------------------------------------------------------------------
% Control Panel
hui.landfig.ctrl_panel.handle = uipanel('Parent',hui.landfig.handle);
set(hui.landfig.ctrl_panel.handle,'Title','Landmark Aquisition');
set(hui.landfig.ctrl_panel.handle,'TitlePosition','CenterTop');
set(hui.landfig.ctrl_panel.handle,'HighlightColor',[0.8 0.8 0.8]);
set(hui.landfig.ctrl_panel.handle,'ShadowColor',[0.6 0.6 0.6]);

% ===================================================================
% Create ctrl_panel table contents
% ===================================================================

% -------------------------------------------------------------------
% Radio buttons

hui.landfig.ctrl_panel.radioBG = uibuttongroup('Parent',hui.landfig.ctrl_panel.handle);
% Create two radio buttons in the button group.
hui.landfig.ctrl_panel.radio0 = uicontrol('Style','Radio','String','Current Frame',...
  'parent',hui.landfig.ctrl_panel.radioBG,'HandleVisibility','on');
hui.landfig.ctrl_panel.radio2 = uicontrol('Style','Radio','String','None',...
  'parent',hui.landfig.ctrl_panel.radioBG,'HandleVisibility','on');

% Initialize some button group properties.
set(hui.landfig.ctrl_panel.radioBG,'SelectionChangeFcn',@landmark_frameRB);
set(hui.landfig.ctrl_panel.radioBG,'SelectedObject',hui.landfig.ctrl_panel.radio0);  % Selection
set(hui.landfig.ctrl_panel.radioBG,'BorderType','none');

% -------------------------------------------------------------------
% Landmark Control Buttons Table
hui.landfig.ctrl_panel.rbTable.ui = []; % Parent is a table container
hui.landfig.ctrl_panel.rbTable.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.landfig.ctrl_panel.rbTable.height_margin = NaN*zeros(30,30);
hui.landfig.ctrl_panel.rbTable.false_width = NaN*zeros(30,30);
hui.landfig.ctrl_panel.rbTable.false_height = NaN*zeros(30,30);

row = 1; col = 1;
hui.landfig.ctrl_panel.rbTable.handles{row,col}   = hui.landfig.ctrl_panel.radio0;
hui.landfig.ctrl_panel.rbTable.width(row,col)     = inf;
hui.landfig.ctrl_panel.rbTable.height(row,col)    = 25;
hui.landfig.ctrl_panel.rbTable.false_width(row,col) = 3;

col = col + 1;
hui.landfig.ctrl_panel.rbTable.handles{row,col}   = hui.landfig.ctrl_panel.radio2;
hui.landfig.ctrl_panel.rbTable.width(row,col)     = inf;
hui.landfig.ctrl_panel.rbTable.height(row,col)    = 25;
hui.landfig.ctrl_panel.rbTable.false_width(row,col) = 3;

clear row col

% -------------------------------------------------------------------
% Landmarks Listbox
hui.landfig.ctrl_panel.landmarksLB = uicontrol('Parent',hui.landfig.ctrl_panel.handle);
menuString = {};

if ~isempty(gCtrl.source.landmarks.fn) && exist(gCtrl.source.landmarks.fn,'file')
  load(gCtrl.source.landmarks.fn);
  if ~isfield(landmarks,{'frm_id','rbin_start','rbin_stop','rline_start',...
      'rline_stop','gpstime_start','gpstime_stop','radar_name','season_name'...
      'notes','type'})
    landmarks = '';
    save('-v6',gCtrl.source.landmarks.fn,'landmarks')
  end
  
  if ~isempty(landmarks)
    for landmark_idx = 1:length(landmarks)
      menuString{landmark_idx} = horzcat(landmarks(landmark_idx).frm_id,' - ',landmarks(landmark_idx).type);
    end
  end
  if ~isempty(landmarks)
    set(hui.landfig.ctrl_panel.landmarksLB,'String',menuString);
    set(hui.landfig.ctrl_panel.landmarksLB,'Value',1);
  end
else
  warning('Landmark file (%s) does not exist.\n', gCtrl.source.landmarks.fn)
  set(hui.landfig.ctrl_panel.landmarksLB,'Value',NaN);
  landmarks = [];
end

set(hui.landfig.ctrl_panel.landmarksLB,'Style','listbox');
set(hui.landfig.ctrl_panel.landmarksLB,'HorizontalAlignment','Center');
set(hui.landfig.ctrl_panel.landmarksLB,'FontName','fixed');
set(hui.landfig.ctrl_panel.landmarksLB,'Callback',@landmarkLB_callback);
clear menuString;

% -------------------------------------------------------------------
% Notes box
hui.landfig.ctrl_panel.notesBox = uicontrol('Style','text');%,'Parent',hui.landfig.ctrl_panel.handle);
set(hui.landfig.ctrl_panel.notesBox,'Value',1);
set(hui.landfig.ctrl_panel.notesBox,'String','Notes (save automatically):');
set(hui.landfig.ctrl_panel.notesBox,'HorizontalAlignment','Left');

hui.landfig.ctrl_panel.notesBox_TE = uicontrol('Parent',hui.landfig.ctrl_panel.handle);
set(hui.landfig.ctrl_panel.notesBox_TE,'Style','edit','Max',2);
set(hui.landfig.ctrl_panel.notesBox_TE,'HorizontalAlignment','Left');
set(hui.landfig.ctrl_panel.notesBox_TE,'Callback',@notesBox_TE_callback);
set(hui.landfig.ctrl_panel.notesBox_TE,'ForegroundColor',[0 0 0])
if isempty(landmarks)
  set(hui.landfig.ctrl_panel.notesBox_TE,'Enable','off');
else
  set(hui.landfig.ctrl_panel.notesBox_TE,'String',landmarks(1).notes);
end


% -------------------------------------------------------------------
% Tool buttons
hui.landfig.ctrl_panel.delete_landPB = uicontrol('Parent',hui.landfig.ctrl_panel.handle);
set(hui.landfig.ctrl_panel.delete_landPB,'Style','PushButton');
set(hui.landfig.ctrl_panel.delete_landPB,'String','Delete Landmark');
set(hui.landfig.ctrl_panel.delete_landPB,'Callback',@landmark_deletePB);

hui.landfig.ctrl_panel.typePM = uicontrol('Parent',hui.landfig.ctrl_panel.handle);
menuString = {};
menuString{1} = 'lake';
menuString{2} = 'mountain';
menuString{3} = 'ocean';
menuString{4} = 'river';
menuString{5} = 'ice';
menuString{6} = 'fjord';
menuString{7} = 'calving front';
menuString{8} = 'surface rocks';
menuString{9} = 'bottom';
menuString{10} = 'clutter';
menuString{11} = 'noise';
set(hui.landfig.ctrl_panel.typePM,'String',menuString);
set(hui.landfig.ctrl_panel.typePM,'Value',1);
clear menuString;
set(hui.landfig.ctrl_panel.typePM,'Style','popupmenu');
set(hui.landfig.ctrl_panel.typePM,'HorizontalAlignment','Center');
set(hui.landfig.ctrl_panel.typePM,'FontName','fixed');
set(hui.landfig.ctrl_panel.typePM,'Callback',@landmark_typePM_callback);

% -------------------------------------------------------------------
% Landmark Control Buttons Table
hui.landfig.ctrl_panel.butTable.ui = []; % Parent is a table container
hui.landfig.ctrl_panel.butTable.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.landfig.ctrl_panel.butTable.height_margin = NaN*zeros(30,30);
hui.landfig.ctrl_panel.butTable.false_width = NaN*zeros(30,30);
hui.landfig.ctrl_panel.butTable.false_height = NaN*zeros(30,30);


row = 1; col = 1;
hui.landfig.ctrl_panel.butTable.handles{row,col}   = hui.landfig.ctrl_panel.delete_landPB;
hui.landfig.ctrl_panel.butTable.width(row,col)     = inf;
hui.landfig.ctrl_panel.butTable.height(row,col)    = 25;
hui.landfig.ctrl_panel.butTable.false_width(row,col) = 3;

col = col + 1;
hui.landfig.ctrl_panel.butTable.handles{row,col}   = hui.landfig.ctrl_panel.typePM;
hui.landfig.ctrl_panel.butTable.width(row,col)     = inf;
hui.landfig.ctrl_panel.butTable.height(row,col)    = 25;
hui.landfig.ctrl_panel.butTable.false_width(row,col) = 3;

clear row col

% ===================================================================
% Setup ctrl_panel table
% ===================================================================
hui.landfig.ctrl_panel.table.ui = hui.landfig.ctrl_panel.handle;
hui.landfig.ctrl_panel.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
hui.landfig.ctrl_panel.table.height_margin = NaN*zeros(30,30);
hui.landfig.ctrl_panel.table.false_width = NaN*zeros(30,30);
hui.landfig.ctrl_panel.table.false_height = NaN*zeros(30,30);
hui.landfig.ctrl_panel.table.offset = [0 10];

row = 1; col = 1;
hui.landfig.ctrl_panel.table.handles{row,col}   = hui.landfig.ctrl_panel.rbTable;
hui.landfig.ctrl_panel.table.width(row,col)     = inf;
hui.landfig.ctrl_panel.table.height(row,col)    = 30;
hui.landfig.ctrl_panel.table.false_width(row,col) = 5;

row = row + 1;
hui.landfig.ctrl_panel.table.handles{row,col}   = hui.landfig.ctrl_panel.landmarksLB;
hui.landfig.ctrl_panel.table.width(row,col)     = inf;
hui.landfig.ctrl_panel.table.height(row,col)    = inf;
hui.landfig.ctrl_panel.table.false_width(row,col) = 3;

row = row + 1;
hui.landfig.ctrl_panel.table.handles{row,col}   = hui.landfig.ctrl_panel.butTable;
hui.landfig.ctrl_panel.table.width(row,col)     = inf;
hui.landfig.ctrl_panel.table.height(row,col)    = 30;
hui.landfig.ctrl_panel.table.false_height(row,col) = 0;
hui.landfig.ctrl_panel.table.false_width(row,col) = 0;
hui.landfig.ctrl_panel.table.width_margin(row,col) = 0;
hui.landfig.ctrl_panel.table.height_margin(row,col) = 3;

row = row + 1;
hui.landfig.ctrl_panel.table.handles{row,col}   = hui.landfig.ctrl_panel.notesBox;
hui.landfig.ctrl_panel.table.width(row,col)     = inf;
hui.landfig.ctrl_panel.table.height(row,col)    = 25;
hui.landfig.ctrl_panel.table.false_width(row,col) = 0;
hui.landfig.ctrl_panel.table.false_height(row,col) = 5;

row = row + 1;
hui.landfig.ctrl_panel.table.handles{row,col}   = hui.landfig.ctrl_panel.notesBox_TE;
hui.landfig.ctrl_panel.table.width(row,col)     = inf;
hui.landfig.ctrl_panel.table.height(row,col)    = 50;
hui.landfig.ctrl_panel.table.false_width(row,col) = 0;
hui.landfig.ctrl_panel.table.false_height(row,col) = 3;

clear row col

% ===================================================================
% Create
% ===================================================================
table_draw(hui.landfig.ctrl_panel.table);
set(hui.landfig.handle,'Visible','off'); % Debuggin'
fprintf('  Done (%.1f sec)\n', toc(gCtrl.tic));

% % ===================================================================
% % Create landmark handles
% % ===================================================================
hui.pickfig.landmarks_h = [];

return;

% ===================================================================
% ===================================================================
% Callback Functions
% ===================================================================
% ===================================================================

function windowClose_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

% Check to see if all layers saved
check_var = gCtrl.source.modified(1);
check_idx = 1;
while (gCtrl.source.modified(check_idx) == check_var) && (check_idx < length(gCtrl.source.modified))
  check_idx = check_idx + 1;
end

if (check_idx == length(gCtrl.source.modified))
  saved = true;
else
  saved = false;
end

if ~saved % notsaved
  if ishandle(hui.fig.handle) && ishandle(hui.mapfig.handle) % first prompt
    prompt = questdlg('There are unsaved layers. Do you want to close?','CReSIS picker toolbox','Yes','No','Yes');
  else % already prompted answer yes, close all
    closereq();
    try; delete(hui.fig.handle); end;
    try; delete(hui.mapfig.handle); end;
    try; delete(hui.pickfig.handle); end;
    try; delete(hui.viewfig.handle); end;
    try; delete(hui.landfig.handle); end;
  end
  switch prompt
    case 'Yes'
      % Close all windows
      closereq();
      try; delete(hui.fig.handle); end;
      try; delete(hui.mapfig.handle); end;
      try; delete(hui.pickfig.handle); end;
      try; delete(hui.viewfig.handle); end;
      try; delete(hui.landfig.handle); end;
    case 'No'
      % Do nothing, close the questdlg
      saved = 'false';
    case ''
      % Do nothing, close the questdlg
      saved = 'false';
  end
  
else % allsaved
  % close all
  closereq();
  try; delete(hui.fig.handle); end;
  try; delete(hui.mapfig.handle); end;
  try; delete(hui.pickfig.handle); end;
  try; delete(hui.viewfig.handle); end;
  try; delete(hui.landfig.handle); end;
end

return

function mapsLB_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

picker_map(1);

return

function framesLB_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

if strcmpi(get(1,'SelectionType'),'open')
  % Double click
  %  - On a double click, the callback is called twice: once with
  %    'normal' SelectionType and then with 'open'.
  %  - On double click, load the selection to the picking window
  picker_pick(1);
else
  % Regular click
  picker_map(2,get(hui.fig.ctrl_panel.framesLB,'Value'));
end

return

function sourceLB_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

if strcmpi(get(1,'SelectionType'),'open')
  % Double click
  %  - On a double click, the callback is called twice: once with
  %    'normal' SelectionType and then with 'open'.
  %  - On double click, load the selection to the picking window
  picker_pick(1);
else
  % Regular click
  gCtrl.source.cur_src = get(hui.fig.ctrl_panel.sourceLB,'Value');
end

return

function display_modePM_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

picker_pick(2);
picker_view(2);

return

function colormapPM_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

picker_colormap(hui.pickfig.handle);
picker_colormap(hui.viewfig.handle);

return

function load_pickPB_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

picker_pick(1);

return

function load_viewPB_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

picker_view(1);

return

function ylim_TE_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

picker_pick(3);
picker_view(3);

return

function toolPM_callback(hObj,event)

global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

% Reset tool settings
gCtrl.tool.selection.x = NaN;

if get(hui.fig.ctrl_panel.toolPM,'Value') == 10
  set(hui.landfig.handle,'Visible','on')
end


return

% ===================================================================
% Landmark Callback Functions
% ===================================================================

function landmark_frameRB(hObj,event)
global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

% Radio button mapping
current = get(hui.landfig.ctrl_panel.radio0,'Value');     % Current Frames
none = get(hui.landfig.ctrl_panel.radio2,'Value');    % None
menuString = {};

% Current or None settings
if ishandle(hui.pickfig.handle)
  if current
    % Shown all frame landmark highligthed boxes
    for frm_idx = 1:length(hui.pickfig.landmarks_h)
      set(hui.pickfig.landmarks_h(frm_idx),'Visible','on');
    end
    set(hui.fig.ctrl_panel.toolPM,'Value',10); % set tool to landmark
  elseif none
    % Show no highlight boxes
    for frm_idx = 1:length(hui.pickfig.landmarks_h)
      set(hui.pickfig.landmarks_h(frm_idx),'Visible','off');
    end
    set(hui.fig.ctrl_panel.toolPM,'Value',1); % set tool to browse
  end
end
return

function landmark_CloseRequestFcn(hObj,event)
global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info


if ishandle(hui.fig.ctrl_panel.toolPM)
  if ishandle(hui.fig.handle)
    if get(hui.fig.ctrl_panel.toolPM,'Value') == 10
      set(hui.fig.ctrl_panel.toolPM,'Value',1);
    end
  end
  set(hui.landfig.handle,'Visible','off');
  if ishandle(hui.pickfig.landmarks_h)
    set(hui.pickfig.landmarks_h,'Visible','off');
  end
else
  % do nothing, closed by map or cmd window
end

return

function landmarkLB_callback(hObj,event)
% Listbox of landmarks callback
global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

if isempty(gCtrl.source.landmarks.fn) || ~exist(gCtrl.source.landmarks.fn,'file')
  return;
end

load(gCtrl.source.landmarks.fn);

if strcmpi(get(hui.landfig.handle,'SelectionType'),'open')
  % Double click
  %  - On a double click, the callback is called twice: once with
  %    'normal' SelectionType and then with 'open'.
  %  - On double click, call picker_pick_landmarks
 load_sel = 1;
  current_frame = landmarks(get(hui.landfig.ctrl_panel.landmarksLB,'Value')).frm_id;
  while ~strcmp(current_frame,gCtrl.source.frm_id(load_sel,:))
    load_sel = load_sel + 1;
  end
  gCtrl.source.cur_sel = load_sel;
  picker_pick(1);
else
  % On regular click, update notes
  current_selection = get(hui.landfig.ctrl_panel.landmarksLB,'Value');
  load(gCtrl.source.landmarks.fn);
  
  if ~isempty(landmarks)
    set(hui.landfig.ctrl_panel.notesBox_TE,'String',landmarks(current_selection).notes);
  end
  % And make current landmark solid line style for distinction
  if ~isempty(hui.pickfig.landmarks_h)
    for landmark_idx = 1:length(landmarks)
      if strcmp(landmarks(landmark_idx).frm_id,gCtrl.source.frm_id(gCtrl.source.cur_pick,:))
        if strcmp(landmarks(landmark_idx).season_name,gCtrl.season_name)
          set(hui.pickfig.landmarks_h(landmark_idx),'LineStyle','- -');
        end
      end
    end
    if strcmp(landmarks(current_selection).frm_id,gCtrl.source.frm_id(gCtrl.source.cur_pick,:))
      if strcmp(landmarks(current_selection).season_name,gCtrl.season_name)
        set(hui.pickfig.landmarks_h(current_selection),'LineStyle','-');
      end
    end
  end
end

return

function landmark_deletePB(hObj,event)
% Landmark removal tool
global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

if isempty(gCtrl.source.landmarks.fn) || ~exist(gCtrl.source.landmarks.fn,'file')
  return;
end

load(gCtrl.source.landmarks.fn);
if isempty(landmarks)
  return;
end

current_selection = get(hui.landfig.ctrl_panel.landmarksLB,'Value');

% Delete specific field of landmark mat file
if ~isempty(landmarks)
  current_frame = landmarks(current_selection).frm_id;
  landmarks(current_selection) = '';
end
save('-v6',gCtrl.source.landmarks.fn,'landmarks');

menuString = {};
for landmark_idx = 1:length(landmarks)
  menuString{landmark_idx} = horzcat(landmarks(landmark_idx).frm_id,' - ',landmarks(landmark_idx).type);
end

if isempty(landmarks)
  set(hui.landfig.ctrl_panel.notesBox_TE,'Enable','off')
end
if isempty(menuString)
  set(hui.landfig.ctrl_panel.landmarksLB,'Value',NaN);
  % Disable notes box
  set(hui.landfig.ctrl_panel.notesBox_TE,'Enable','off');
  set(hui.landfig.ctrl_panel.notesBox_TE,'String','');
elseif  current_selection > length(menuString)
  set(hui.landfig.ctrl_panel.landmarksLB,'Value',length(menuString));
end
set(hui.landfig.ctrl_panel.landmarksLB,'String',menuString);

if ishandle(hui.pickfig.handle)
  if strcmp(gCtrl.source.frm_id(gCtrl.source.cur_pick,:),current_frame)
    % Current frame contains the landmark that we just deleted, so
    % we need to delete the handle from the handle list
    if length(menuString) >= 1
      delete(hui.pickfig.landmarks_h(current_selection));
    else
      for landmark_idx = 1:length(hui.pickfig.landmarks_h)
        if ishandle(hui.pickfig.landmarks_h(landmark_idx))
          delete(hui.pickfig.landmarks_h(landmark_idx)); 
        end
      end
    end
    
    hui.pickfig.landmarks_h(current_selection) = [];
  end
  figure(hui.pickfig.handle)
  hold on;
  plot(hui.pickfig.landmarks_h);
  hold off;
end

current_selection = get(hui.landfig.ctrl_panel.landmarksLB,'Value');
if ~isnan(current_selection)
  set(hui.landfig.ctrl_panel.notesBox_TE,'String',landmarks(current_selection).notes);
end

return

function landmark_typePM_callback(hObj,event)
% Landmark type tool selection
global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

if isempty(gCtrl.source.landmarks.fn) || ~exist(gCtrl.source.landmarks.fn,'file')
  return;
end

load(gCtrl.source.landmarks.fn)

new_types = get(hui.landfig.ctrl_panel.typePM,'String');
new_type = char(new_types(get(hui.landfig.ctrl_panel.typePM,'Value')));
current_selection = get(hui.landfig.ctrl_panel.landmarksLB,'Value');
landmarks(current_selection).type = new_type;

menuString = {};
for landmark_idx = 1:length(landmarks)
  menuString{landmark_idx} = horzcat(landmarks(landmark_idx).frm_id,' - ',landmarks(landmark_idx).type);
end
set(hui.landfig.ctrl_panel.landmarksLB,'String',menuString);

save('-v6',gCtrl.source.landmarks.fn,'landmarks');

return

function notesBox_TE_callback(hObj,event)
% Landmark note description text edit field
global hui; % hui: user interface handles
global gCtrl; % gCtrl: contains all the control info

if isempty(gCtrl.source.landmarks.fn) || ~exist(gCtrl.source.landmarks.fn,'file')
  return;
end

load(gCtrl.source.landmarks.fn)

current_selection = get(hui.landfig.ctrl_panel.landmarksLB,'Value');
notes = get(hui.landfig.ctrl_panel.notesBox_TE,'String');
load(gCtrl.source.landmarks.fn);
landmarks(current_selection).notes = notes;

save('-v6',gCtrl.source.landmarks.fn,'landmarks');

return
