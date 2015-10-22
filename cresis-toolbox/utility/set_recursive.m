function set_recursive(h_list,field_name,field_value)
% set_recursive(h_list,field_name,field_value)
%
% Recurses through children of handles in list and changes field name
% to field value. If field value is a function handle, then the function
% handle operates on each of the values.
%
% h_list = list of handles
% field_name = string containing the field name to change
% field_value = one of two:
%   function handle: the current value will be operated on by this function
%   any thing else: the value of the field field_name will be set to field_value
%
% set_recursive(h_axes,'FontSize',8);
% set_recursive(h_axes,'FontWeight','normal');
% dorthset_recursivee_set(h_axes,'String',@dorthe_capitalize_first);
%
% Author: John Paden

%% Apply operation to each handle in h_list
for h_idx = 1:numel(h_list)
  %% Get list of fields for this handle
  h_fields = get(h_list(h_idx));
  
  %% Check for the existence of the field that we want to change, if it
  % exists change it.
  if isfield(h_fields,field_name)
    if strcmp(class(field_value),'function_handle')
      set(h_list(h_idx),field_name,field_value(get(h_list(h_idx),field_name)));
    else
      set(h_list(h_idx),field_name,field_value);
    end
  end
  
  %% Recurse through each of the children of this handle
  h_children = get(h_list(h_idx),'Children');
  set_recursive(h_children,field_name,field_value);
  
  %% Recurse through special children of axes handles
  if isfield(h_fields,'Title')
    h_title = get(h_list(h_idx),'Title');
    set_recursive(h_title,field_name,field_value);
    h_xlabel = get(h_list(h_idx),'XLabel');
    set_recursive(h_xlabel,field_name,field_value);
    h_ylabel = get(h_list(h_idx),'YLabel');
    set_recursive(h_ylabel,field_name,field_value);
  end
end

end
