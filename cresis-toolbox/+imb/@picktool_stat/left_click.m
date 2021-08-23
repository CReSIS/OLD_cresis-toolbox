function cmds = left_click(obj,param)

%% Get search range from tool param window
rbin_range_str = get(obj.panel.manual_rangeTB,'String');
try
  % Assumes the value entered is a matlab expression that can be evaluated
  search_range = eval(sprintf('[%s]', rbin_range_str));
  if length(search_range) == 1
    search_range = -search_range:search_range;
  else
    search_range = round(min(search_range)):round(max(search_range));
  end
catch ME
  warning('Search range parameter is not valid, using default range');
  search_range = -5:5;
end

if isempty(search_range)
  search_range = -5:5;
end

param.search_range = search_range;
image_x = param.image_x;
image_y = param.image_y;
% cur_layers = param.cur_layers;
x = param.x;
y = param.y;
search_range = param.search_range;

% fprintf('Insert point %f, %f\n', x, y);

% Clamp entered points to the x-boundaries
xlims = [min(image_x) max(image_x)];
if x < xlims(1)
  x = xlims(1);
elseif x > xlims(2)
  x = xlims(2);
end
  
% Find the closest point in the image
[~,image_x_idx] = min(abs(image_x-x));
dy = image_y(2)-image_y(1);
image_y_idx = 1 + round((y-image_y(1))/dy);
if image_y_idx < 1 
  max_idx = image_y_idx;
  image_y_idx = 1;
elseif image_y_idx > size(param.image_c,1)
  image_y_idx = size(param.image_c,1);
  max_idx = image_y_idx;
else
  % Prevent search range from searching outside bounds of param.image_c
  search_range = search_range(image_y_idx+search_range >= 1 ...
    & image_y_idx+search_range < size(param.image_c,1));
  [~,max_idx] = max(param.image_c(image_y_idx+search_range,image_x_idx));
  max_idx = search_range(max_idx) + image_y_idx;
end

% [~,point_idxs] = min(abs(param.layer.x-x));
new_y = interp1(1:length(image_y),image_y,max_idx,'linear','extrap');
pnt = [image_x_idx,image_y_idx];
param.val = param.image_c(pnt(2),pnt(1));
fprintf('X = %f Y = %f Val = %f\n', x,new_y,param.val)
cmds = [];

return
