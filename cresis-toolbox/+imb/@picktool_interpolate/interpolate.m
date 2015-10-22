function [vals x_new] = interpolate(obj,x_old,y_old,x_new)
%
% Interpolate for picktool_enter
%

tool_idx = get(obj.panel.interp_mode_pdmenu,'Value');

if length(x_old) ~= length(unique(x_old))
  % This is technically okay, but it is probably caused by a bug elsewhere
  % (i.e. generally we do not want to points with the same GPS time).
  % Type "dbcont" to continue.
  warning('Two points were found with exactly the same GPS time.');  
  [x_old unique_idxs] = unique(x_old);
  y_old = y_old(unique_idxs);
end

switch tool_idx
  case 1 % linear interpolate
    vals = interp1(x_old,y_old,x_new,'linear','extrap');
  case 2 % spline interpolate
    vals = interp1(x_old,y_old,x_new,'spline','extrap');
  case 3 % max track interpolate
    
  case 4 % leading edge track interpolate
    
end