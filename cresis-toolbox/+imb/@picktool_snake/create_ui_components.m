function create_ui_components(obj)

% create_ui_components(obj)
%
% Creates components for the snake param window's UI when the echogram
% window is first opened.
%

set(obj.h_fig,'DockControls','off')
set(obj.h_fig,'NumberTitle','off');
if strcmpi(class(obj.h_fig),'double')
  set(obj.h_fig,'Name',sprintf('%d: %s tool parameters', obj.h_fig, obj.tool_name_title));
else
  set(obj.h_fig,'Name',sprintf('%d: %s tool parameters', obj.h_fig.Number, obj.tool_name_title));
end
set(obj.h_fig,'ToolBar','none');
set(obj.h_fig,'MenuBar','none'); 
set(obj.h_fig,'KeyPressFcn',@obj.key_press);
set(obj.h_fig,'CloseRequestFcn',@obj.close_win);
%set(obj.h_fig,'Resize','off');

return

