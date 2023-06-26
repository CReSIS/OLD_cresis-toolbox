function cmds = right_click(obj,param)
% cmds = right_click(obj,param)

cmds = [];
if param.echowin.cursor_mode == false
	param.echowin.update_cursor(param.x(1),param.y(1),true);
end
notify(obj,'ascope_memory');


    