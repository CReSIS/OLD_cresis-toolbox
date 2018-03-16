function [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,cur_layer)
% [manual_idxs,auto_idxs,point_idxs] = find_matching_pnts(obj,param,cur_layer)
%
% Utility function to find which point indexes were selected (useful
% for child classes' left_click_and_drag function)
%
% param is standard param to left_click_and_drag plus these fields:
%  x_bounds: Integer which controls how param.x is interpretted
%    0: no x-bounds
%    1: must be a manual point
%    2: must lie in the box drawn by user
%    3: conditions 1 and 2
%  y_bounds = Same options as x_bounds for param.y

mask = logical(ones(size(param.layer.x)));
if mod(param.x_bounds,2)
  % Check for "must be a manual point"
  mask = mask & param.layer.type{cur_layer} == 1 & ~isnan(param.layer.y{cur_layer});
end
if mod(floor(param.x_bounds/2),2)
  % Check for "must lie in the box drawn by user"
  mask = mask & param.layer.x >= param.x(1) & param.layer.x <= param.x(2);
end
if mod(param.y_bounds,2)
  mask = mask & param.layer.type{cur_layer} == 1;
end
if mod(floor(param.y_bounds/2),2)
  mask = mask & param.layer.y{cur_layer} >= param.y(1) & param.layer.y{cur_layer} <= param.y(2);
end

mask(find(mask,1):find(mask,1,'last')) = true;
manual_idxs = find(mask & ~isnan(param.layer.y{cur_layer}) & param.layer.type{cur_layer} == 1);
auto_idxs = find(mask & (isnan(param.layer.y{cur_layer}) | param.layer.type{cur_layer} ~= 1));
point_idxs = find(mask);

end
