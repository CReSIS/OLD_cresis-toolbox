function close_win(obj,varargin)
% imagewin.close_win(obj,varargin)

if obj.hide_only
  set(obj.h_fig,'Visible','off');
else
  delete(obj);
end

return;
