function [modified] = ismodified(obj)
% Check to see if there are any commands in the undo stack that have not
% been saved. Returns true if they have, otherwise returns false.

if isempty(obj.stack) || obj.pointer == 0
  modified = false;
else
  modified = true;
end

return;