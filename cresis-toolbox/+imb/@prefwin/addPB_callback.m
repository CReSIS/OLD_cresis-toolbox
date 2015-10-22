function addPB_callback(obj, hObj, event)

selectedString = get(obj.top_ctrl_panel.selectedLB,'String');
menuString = get(obj.top_ctrl_panel.seasonLB,'String');
menuString_logic = true(1,length(menuString));
selected_season_idx = get(obj.top_ctrl_panel.seasonLB,'Value');


if strcmp(menuString{1},'none')
  %do nothing
  return;
elseif strcmp(selectedString{1},'none')
  selectedString(1) = menuString(selected_season_idx);
  menuString_logic(selected_season_idx)=0;
  menuString=menuString(menuString_logic);
elseif ~strcmp(selectedString{1},'none') && ~any(strcmp(selectedString,menuString{selected_season_idx}))
  selectedString(end+1)=menuString(selected_season_idx);
  menuString_logic(selected_season_idx)=0;
  menuString=menuString(menuString_logic);
end

if isempty(menuString)
  menuString={'none'};
end

% Filter available maps based on selected system and season
systems = get(obj.bottom_ctrl_panel.sensorsLB,'string');
system = systems{get(obj.bottom_ctrl_panel.sensorsLB,'value')};
selected_idxs = zeros(size(selectedString));
for idx = 1:length(selectedString)
  % For each season selected, find its index in obj.seasons
  found = false;
  for search_idx = 1:length(obj.seasons)
    if strcmp(obj.seasons{search_idx},selectedString{idx}) ...
        && strcmp(obj.systems{search_idx},system)
      selected_idxs(idx) = search_idx;
      found = true;
      break;
    end
  end
  if ~found
    warning('This should never happen, selection was not found in the list');
    keyboard
  end
end
selected_idxs = sort(selected_idxs);
zones = obj.locations(selected_idxs);

% val = 0: nothing selected (not possible in this callback)
% val = 1: only arctic
% val = 2: only antarctic
% val = 3: both arctic and antarctic
val = 0;
if any(strcmp(zones,'arctic'))
  val = val+1;
end
if any(strcmp(zones,'antarctic'))
  val = val+2;
end
% keep the old selected map to match later
old_sel_idx = get(obj.bottom_ctrl_panel.mapsPM,'value');
old_maps = get(obj.bottom_ctrl_panel.mapsPM,'String');
old_sel_map = old_maps{old_sel_idx};

switch val
  case 1
    good_idxs = logical(zeros(1,length(obj.full_maplist)));
    for idx=1:length(obj.full_maplist)
      if strncmp('arctic',obj.full_maplist{idx},6)
        good_idxs(idx) = 1;
      else
        good_idxs(idx) = 0;
      end
    end
    new_mapstring = obj.full_maplist(good_idxs);
    new_mapstring = ['None' new_mapstring];
    sel_idx = find(strcmp(old_sel_map,new_mapstring),'1');
    if isempty(sel_idx)
      sel_idx = 1;
    end
    set(obj.bottom_ctrl_panel.mapsPM,'Value',sel_idx);
    set(obj.bottom_ctrl_panel.mapsPM,'String',new_mapstring);
  case 2
    good_idxs = logical(zeros(1,length(obj.full_maplist)));
    for idx=1:length(obj.full_maplist)
      if strncmp('antarctic',obj.full_maplist{idx},9)
        good_idxs(idx) = 1;
      else
        good_idxs(idx) = 0;
      end
    end
    new_mapstring = obj.full_maplist(good_idxs);
    new_mapstring = ['None' new_mapstring];
    sel_idx = find(strcmp(old_sel_map,new_mapstring),'1');
    if isempty(sel_idx)
      sel_idx = 1;
    end
    set(obj.bottom_ctrl_panel.mapsPM,'Value',sel_idx);
    set(obj.bottom_ctrl_panel.mapsPM,'String',new_mapstring);
  case 3
    sel_idx = find(strcmp(old_sel_map,obj.full_maplist),'1');
    if isempty(sel_idx)
      sel_idx = 1;
    end
    set(obj.bottom_ctrl_panel.mapsPM,'Value',sel_idx);
    set(obj.bottom_ctrl_panel.mapsPM,'String',obj.full_maplist);
  case default
    error('Error in @prefwin\addPB_callback.m');
end

menuString = sort(menuString);
[selectedString sort_idxs] = sort(selectedString);

set(obj.top_ctrl_panel.seasonLB,'String',menuString);
if strcmpi(menuString,'none')
  set(obj.top_ctrl_panel.seasonLB,'Value',1);
elseif selected_season_idx > length(menuString)
  set(obj.top_ctrl_panel.seasonLB,'Value',length(menuString));
end

set(obj.top_ctrl_panel.selectedLB,'Value',find(sort_idxs==length(selectedString)));
set(obj.top_ctrl_panel.selectedLB,'String',selectedString);

return
