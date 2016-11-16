function sourceLB_callback(obj,status,event)

h_status = get(status);
if isfield(h_status,'Label') && strcmp(get(status,'Label'),'Remove')
  selected_items = get(obj.h_gui.sourceLB,'Value');
  items = get(obj.h_gui.sourceLB,'String');
  mask = logical(zeros(size(items)));
  mask(selected_items) = 1;
  mask = ~mask;
  items = items(mask);
  set(obj.h_gui.sourceLB,'String',items);
  set(obj.h_gui.sourceLB,'Value',[]);
  
elseif isfield(h_status,'Label') && strcmp(get(status,'Label'),'Add')
  
  prompt = {'Enter new sources (E.g. CSARP_standard):'};
  dlg_title = 'Input new sources';
  num_lines = 1;
  def = {'CSARP_'};
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
