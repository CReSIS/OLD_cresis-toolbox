function create_ui(obj)

%% Get the system info from the database (systems,seasons,locations)
fprintf('  Querying database to get data list (%s)\n', datestr(now,'HH:MM:SS'));
[status,data] = opsGetSystemInfo();
unique_system_list = unique(data.properties.systems);
obj.profile = data.opsProfile;
obj.systems = data.properties.systems;
obj.seasons = data.properties.seasons;
obj.locations = data.properties.locations;

%% Connect to CReSIS WMS and set obj.full_maplist
fprintf('  Querying WMS to get map list (%s)\n', datestr(now,'HH:MM:SS'));
opsCmd; % Populate gOPS variable
wms = WebMapServer(sprintf('%swms/',gOps.geoServerUrl));
try
  wms_capabilities = wms.getCapabilities();
catch ME
  rethrow(ME);
end
wms_layers = wms_capabilities.Layer;

% Create a list of all layers excluding "*line_paths*",
% "*data_quality*","*data_coverage*", "*crossover_errors*","*data_elevation*"
obj.wms_maps = {};
obj.wms_maps{end+1} = 'arctic:Google';
obj.wms_maps{end+1} = 'antarctic:Google';
for idx = 1:length(wms_layers)
  match1_idx = strfind(wms_layers(idx).LayerName,'line_paths');
  match2_idx = strfind(wms_layers(idx).LayerName,'crossover_errors');
  match3_idx = strfind(wms_layers(idx).LayerName,'data_quality');
  match4_idx = strfind(wms_layers(idx).LayerName,'data_coverage');
  match5_idx = strfind(wms_layers(idx).LayerName,'data_elevation');
  if isempty(match1_idx) && isempty(match2_idx) && isempty(match3_idx) && isempty(match4_idx) && isempty(match5_idx)
    obj.wms_maps{end+1} = wms_layers(idx).LayerName;
  end
end

%%====================================================================
%% Create User Interface for prefwin class
%%====================================================================

%% Set prefwin figure properties
set(obj.h_fig,'Position',[obj.default_params.x obj.default_params.y obj.default_params.w obj.default_params.h]);
set(obj.h_fig,'DockControls','off')
set(obj.h_fig,'NumberTitle','off');
if strcmpi(class(obj.h_fig),'double')
  set(obj.h_fig,'Name',sprintf('%d: preference',obj.h_fig));
else
  set(obj.h_fig,'Name',sprintf('%d: preference',obj.h_fig.Number));
end
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none');
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);

%% Create the widgets

% Layer selection class (populate later from preference file)
obj.h_gui.layers = selectionbox(obj.h_fig,'Layers',[],1);
set(obj.h_gui.layers.h_list_available,'TooltipString','Available layers (double or right click to select).');
set(obj.h_gui.layers.h_list_selected,'TooltipString','Selected layers (double or right click to remove).');

uimenu(obj.h_gui.layers.h_list_availableCM, 'Label', 'New', 'Callback', @obj.layers_callback);
uimenu(obj.h_gui.layers.h_list_availableCM, 'Label', 'Delete', 'Callback', @obj.layers_callback);
uimenu(obj.h_gui.layers.h_list_availableCM, 'Label', 'Rename', 'Callback', @obj.layers_callback);
uimenu(obj.h_gui.layers.h_list_availableCM, 'Label', 'Refresh', 'Callback', @obj.layers_callback);

% Season selection class (populate later from preference file)
obj.h_gui.seasons = selectionbox(obj.h_fig,'Seasons',[],1);
set(obj.h_gui.seasons.h_list_available,'TooltipString','Available seasons (double click or right click to select).');
set(obj.h_gui.seasons.h_list_selected,'TooltipString','Selected seasons(double click or right click to remove).');

% System list box label
obj.h_gui.systemsText = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.systemsText,'Style','Text');
set(obj.h_gui.systemsText,'String','Systems');

% Source list box label
obj.h_gui.sourceText = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.sourceText,'Style','Text');
set(obj.h_gui.sourceText,'String','Echogram Sources');

% System list box
obj.h_gui.systemsLB = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.systemsLB,'String',unique_system_list);
set(obj.h_gui.systemsLB,'Style','listbox');
set(obj.h_gui.systemsLB,'HorizontalAlignment','Center');
set(obj.h_gui.systemsLB,'FontName','fixed');
set(obj.h_gui.systemsLB,'Callback',@obj.systemsLB_callback);
set(obj.h_gui.systemsLB,'Min',1); % One must always be selected
set(obj.h_gui.systemsLB,'TooltipString','Systems (choose one)');

% Source list box (populate later from preference file)
obj.h_gui.sourceLB = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.sourceLB,'String','');
set(obj.h_gui.sourceLB,'Style','listbox');
set(obj.h_gui.sourceLB,'HorizontalAlignment','Center');
set(obj.h_gui.sourceLB,'FontName','fixed');
set(obj.h_gui.sourceLB,'Callback',@obj.sourceLB_callback);
set(obj.h_gui.sourceLB,'TooltipString',...
  sprintf('List of echogram sources to load\n Left click to select\n Right click to add or remove'));
set(obj.h_gui.sourceLB,'Min',1); % One must always be selected
set(obj.h_gui.sourceLB,'Max',1e9); % Allow multiple selections

% Source list box context menu
obj.h_gui.sourceCM = uicontextmenu;
% Define the context menu items and install their callbacks
obj.h_gui.sourceCM_item1 = uimenu(obj.h_gui.sourceCM, 'Label', 'Add', 'Callback', @obj.sourceLB_callback);
obj.h_gui.sourceCM_item2 = uimenu(obj.h_gui.sourceCM, 'Label', 'Remove', 'Callback', @obj.sourceLB_callback);
set(obj.h_gui.sourceLB,'uicontextmenu',obj.h_gui.sourceCM)

% Map Popup Menu (populate later from preference file)
obj.h_gui.mapsPM = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.mapsPM,'String',obj.wms_maps);
set(obj.h_gui.mapsPM,'Value',1);
set(obj.h_gui.mapsPM,'Style','popupmenu');
set(obj.h_gui.mapsPM,'HorizontalAlignment','Center');
set(obj.h_gui.mapsPM,'FontName','fixed');
set(obj.h_gui.mapsPM,'TooltipString','Available maps (select one which matches seasons'' location).');

% Map flightline/vectors Popup Menu (populate later from preference file)
obj.h_gui.flightlinesPM = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.flightlinesPM,'String',{'Regular Flightlines','Quality Flightlines','Coverage Flightlines', 'Crossover Errors','Bed Elevation'});
set(obj.h_gui.flightlinesPM,'Value',1);
set(obj.h_gui.flightlinesPM,'Style','popupmenu');
set(obj.h_gui.flightlinesPM,'HorizontalAlignment','Center');
set(obj.h_gui.flightlinesPM,'FontName','fixed');
set(obj.h_gui.flightlinesPM,'TooltipString','Available flightlines.');

% LayerSource pop up menu (populate later from preference window)%%
obj.h_gui.LayerSourcePM = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.LayerSourcePM,'String',{'layerdata','OPS'});
set(obj.h_gui.LayerSourcePM,'Value',1);
set(obj.h_gui.LayerSourcePM,'Style','popupmenu');
set(obj.h_gui.LayerSourcePM,'HorizontalAlignment','Center');
set(obj.h_gui.LayerSourcePM,'FontName','fixed');
set(obj.h_gui.LayerSourcePM,'TooltipString','Available layer sources (select one)');
set(obj.h_gui.LayerSourcePM,'Callback',@obj.LayerSourcePM_callback);
%V = get(obj.h_gui.LayerSourcePM,'Value');

% layerdata sources pop up menu (populate later from preference window)%%
obj.h_gui.layerDataSourcePM = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.layerDataSourcePM,'String',{'layerData','CSARP_post/layerData','CSARP_layerdata'});
set(obj.h_gui.layerDataSourcePM,'Value',1);
set(obj.h_gui.layerDataSourcePM,'Style','popupmenu');
set(obj.h_gui.layerDataSourcePM,'HorizontalAlignment','Center');
set(obj.h_gui.layerDataSourcePM,'FontName','fixed');
set(obj.h_gui.layerDataSourcePM,'TooltipString','Available layerdata sources (select one)');

% Okay Button
obj.h_gui.okPB = uicontrol('Parent',obj.h_fig);
set(obj.h_gui.okPB,'Style','PushButton');
set(obj.h_gui.okPB,'String','OK');
set(obj.h_gui.okPB,'Callback',@obj.okPB_callback);

%% Create the table
obj.h_gui.table.ui=obj.h_fig;
obj.h_gui.table.width_margin = NaN*zeros(30,30); % Just make these bigger than they have to be
obj.h_gui.table.height_margin = NaN*zeros(30,30);
obj.h_gui.table.false_width = NaN*zeros(30,30);
obj.h_gui.table.false_height = NaN*zeros(30,30);
obj.h_gui.table.offset = [0 0];

row = 1; col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.layers.h_text;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.width(row,col)     = 0;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;%%
obj.h_gui.table.handles{row,col}   = obj.h_gui.LayerSourcePM;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 25;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col =2;%%
obj.h_gui.table.handles{row,col}   = obj.h_gui.layerDataSourcePM;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 25;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.layers.h_list_available;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 80;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.layers.h_list_selected;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 80;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.seasons.h_text;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.width(row,col)     = 0;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.seasons.h_list_available;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = inf;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.seasons.h_list_selected;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = inf;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.systemsText;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.sourceText;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 20;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.systemsLB;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 80;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.sourceLB;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 80;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;
obj.h_gui.table.handles{row,col}   = obj.h_gui.mapsPM;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 22;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.flightlinesPM;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 22;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

row = row + 1; col = 1;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 25;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

col = 2;
obj.h_gui.table.handles{row,col}   = obj.h_gui.okPB;
obj.h_gui.table.width(row,col)     = inf;
obj.h_gui.table.height(row,col)    = 25;
obj.h_gui.table.width_margin(row,col) = 1;
obj.h_gui.table.height_margin(row,col) = 1;

clear row col
table_draw(obj.h_gui.table);

%% Set settings passed in by obj.default_params
if isfield(obj.default_params,'system') && ischar(obj.default_params.system)
  match_idx = find(strcmp(obj.default_params.system,get(obj.h_gui.systemsLB,'String')));
else
  match_idx = [];
end
if isempty(match_idx)
  set(obj.h_gui.systemsLB,'Value',1);
else
  set(obj.h_gui.systemsLB,'Value',match_idx);
end
obj.systemsLB_callback;

% Select the default seasons
obj.h_gui.seasons.set_selected(obj.default_params.season_names,true);

% Select the default layers
obj.h_gui.layers.set_selected(obj.default_params.layer_names,true);

%% Set default map
if isfield(obj.default_params,'map_name') && ischar(obj.default_params.map_name)
  match_idx = find(strcmp(obj.default_params.map_name,get(obj.h_gui.mapsPM,'String')));
else
  match_idx = [];
end
if isempty(match_idx)
  set(obj.h_gui.mapsPM,'Value',1);
else
  set(obj.h_gui.mapsPM,'Value',match_idx);
end

%% Set default flightlines
if isfield(obj.default_params,'flightlines') && ischar(obj.default_params.flightlines)
  match_idx = find(strcmp(obj.default_params.flightlines,get(obj.h_gui.flightlinesPM,'String')));
else
  match_idx = [];
end
if isempty(match_idx)
  set(obj.h_gui.flightlinesPM,'Value',1);
else
  set(obj.h_gui.flightlinesPM,'Value',match_idx);
end

%% Set default echogram sources
if isfield(obj.default_params,'sources')
  set(obj.h_gui.sourceLB,'String',obj.default_params.sources);
end

%% Set default layer source
if isfield(obj.default_params,'LayerSource') && ischar(obj.default_params.LayerSource)
  match_idx = find(strcmp(obj.default_params.LayerSource,get(obj.h_gui.LayerSourcePM,'String')));
else
  match_idx = [];
end
if isempty(match_idx)
  set(obj.h_gui.LayerSourcePM,'Value',1);
else
  set(obj.h_gui.LayerSourcePM,'Value',match_idx);
end

%% Set default layerdata source
if isfield(obj.default_params,'layerDataSource') && ischar(obj.default_params.layerDataSource)
  match_idx = find(strcmp(obj.default_params.layerDataSource,get(obj.h_gui.layerDataSourcePM,'String')));
else
  match_idx = [];
end
if isempty(match_idx)
  temp = get(obj.h_gui.LayerSourcePM,'String');
  LayerSource = temp{get(obj.h_gui.LayerSourcePM,'Value')};
  if strcmpi(LayerSource,'OPS')
    set(obj.h_gui.layerDataSourcePM,'Enable','off');
    set(obj.h_gui.layerDataSourcePM,'Value',1);
  elseif strcmpi(LayerSource,'layerdata')
    set(obj.h_gui.layerDataSourcePM,'Enable','on');
    set(obj.h_gui.layerDataSourcePM,'Value',1);
  end
else
  temp = get(obj.h_gui.LayerSourcePM,'String');
  LayerSource = temp{get(obj.h_gui.LayerSourcePM,'Value')};
  if strcmpi(LayerSource,'OPS')
    set(obj.h_gui.layerDataSourcePM,'Enable','off');
    set(obj.h_gui.layerDataSourcePM,'Value',match_idx);
  elseif strcmpi(LayerSource,'layerdata')
    set(obj.h_gui.layerDataSourcePM,'Enable','on');
    set(obj.h_gui.layerDataSourcePM,'Value',match_idx);
  end
end

return
