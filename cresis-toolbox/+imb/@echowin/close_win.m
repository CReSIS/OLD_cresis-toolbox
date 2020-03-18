function close_win(obj,varargin)
% close_win(obj,varargin)
%
% We don't delete the object or close the figure here, that is done in the
% parent classes close function which is called when we notify close_window.

cancel_operation = obj.undo_stack_modified_check();

if ~cancel_operation
  % Delete the map preferences window
  notify(obj,'close_window');
end
