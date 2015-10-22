function cmds = left_click(obj,param)

% grab parameters that will be used
image_y = param.image_y;
A = param.image_c;
cmds = [];

set(obj.h_fig,'Units','Pixels');
cur_pos = get(obj.h_fig,'Position');
if ~param.keep_tool_pos
  % grab the echowin's position info (passed in)
  tool_pos = [param.tool_x param.tool_y-cur_pos(4) cur_pos(3) cur_pos(4)];  
  set(obj.h_fig,'Position',tool_pos);  
end
set(obj.h_fig,'Visible','on');

% Setup A-scope view
% obtain yaxis choice of the EG
% y_choice = 1;
% if y_choice == 1
%   x_label = 'Time (us)';
% elseif y_choice == 2
%   x_label = 'Elevation (m)';
% elseif y_choice == 3
%   x_label = 'Depth (m)';
% elseif y_choice == 4
%   x_label = 'Range bin';
% end
% cur_xlabel_h = get(obj.h_axes,'XLabel');
% check if the yaxis label is different from the current ascop xaxis label
% if ~strcmp(x_label,get(cur_xlabel_h,'String'))
%   set(obj.h_plot,'XData',[]);
%   set(obj.h_plot,'YData',[]);
% else
%   hold on;
% end

% Get the closest range line
[tmp cursor_idx] = min(abs(param.x - param.image_x));
% Set the A-scope/browse-param window plot
set(obj.h_plot,'XData',param.image_y,'YData',param.image_c(:,cursor_idx));
grid(obj.h_axes,'on');

% xlabel(obj.h_axes,x_label);
% ylabel(obj.h_axes,'Relative power (dB)');

return
