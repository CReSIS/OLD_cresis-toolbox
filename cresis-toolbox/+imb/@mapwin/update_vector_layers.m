function update_vector_layers(obj)
% update_vector_layers(obj)
%
% Update vector graphics for mapwin map
% 1. Draws red line for selected frame
%   obj.cur_sel.X,obj.cur_sel.Y
% 2. Draws green lines associated with each echowin window
%    obj.draw_echowin_flightline()
% 3. Draws cursors
%    obj.draw_cursors()
%
% Only plots flightlines if a selection has been made.

selection_color = [1 0 0];

%------------------------------------------------------------------------
% Draw current selection
set(obj.map_panel.h_cur_sel,{'XData','YData'},{obj.cur_sel.X,obj.cur_sel.Y});

%------------------------------------------------------------------------
% Update echowin flightlines
for echowin_idx = 1:length(obj.echowin_list)
obj.update_echowin_flightlines(obj.echowin_list(echowin_idx));
end

%------------------------------------------------------------------------
% Update echowin cursors
obj.update_cursors();


return
