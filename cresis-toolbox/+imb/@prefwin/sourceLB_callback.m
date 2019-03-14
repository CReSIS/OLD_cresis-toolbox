function sourceLB_callback(obj,status,event)

% Turn object into struct so we can check existence of fields
h_status = get(status);

% Determine the name of the field that contains the selected menu item
% 9.3 was first??? version using "Text" instead of "Label"
if 0
  matlab_ver = ver('matlab');
  use_text = str2double(matlab_ver.Version) >= 9.3;
elseif 1
  if isfield(h_status,'Text')
    use_text = 1;
  else
    use_text = 0;
  end
end

if use_text
  menu_field_name = 'Text';
else
  menu_field_name = 'Label';
end

if isfield(h_status,menu_field_name) && strcmp(get(status,menu_field_name),'Remove')
  selected_items = get(obj.h_gui.sourceLB,'Value');
  items = get(obj.h_gui.sourceLB,'String');
  mask = logical(zeros(size(items)));
  mask(selected_items) = 1;
  mask = ~mask;
  items = items(mask);
  set(obj.h_gui.sourceLB,'String',items);
  set(obj.h_gui.sourceLB,'Value',[]);
  
elseif isfield(h_status,menu_field_name) && strcmp(get(status,menu_field_name),'Add')
  
  prompt = {'Enter new sources (E.g. standard):'};
  dlg_title = 'Input new sources';
  num_lines = 1;
  def = {''};
  answer = inputdlg(prompt,dlg_title,num_lines,def);
  
  if ~isempty(answer) && ~isempty(answer{1}) && ischar(answer{1})
    answer = answer{1};
    selected_items = get(obj.h_gui.sourceLB,'Value');
    items = get(obj.h_gui.sourceLB,'String');
    items{end+1} = answer;
    [items new_idxs] = unique(items);
    set(obj.h_gui.sourceLB,'String',items);
    selected_items = new_idxs(selected_items);
    set(obj.h_gui.sourceLB,'Value',selected_items);
  end

end
return
