function cmds = left_click_and_drag(obj,param)
% cmds = left_click_and_drag(obj,param)
%
% Stat tool

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

image_x = param.image_x;
image_y = param.image_y;

x = param.x;
y = param.y;
cmds = [];

%% Get tool stat method
stat_idx = get(obj.panel.stat_modePM,'Value');
stat_idx_y = get(param.echowin.left_panel.yaxisPM,'Value');
param.x_bounds = 3;
param.y_bounds = 1;

% Clamp entered points to the x-boundaries
xlims = [min(image_x) max(image_x)];
if x(1) < xlims(1)
  x(1) = xlims(1);
elseif x(2) > xlims(2)
  x(2) = xlims(2);
end
% f = max(min_allowed_f,f)
% f(f< min_allowed_f) = min_allowed_f


%% Find the closest point in the image - initial point
[~,image_x_idx_init] = min(abs(image_x-x(1)));
dy = image_y(2)-image_y(1);
if stat_idx_y == 2
  image_y_idx_init = 1 + round((y(2)-image_y(1))/dy);
else
  image_y_idx_init = 1 + round((y(1)-image_y(1))/dy);
end

if image_y_idx_init < 1
  max_idx_init = image_y_idx_init;
  image_y_idx_init = 1;
elseif image_y_idx_init > size(param.image_c,1)
  max_idx_init = image_y_idx_init;
  image_y_idx_init = (size(param.image_c,1));
else
  % Prevent search range from searching outside bounds of param.image_c
  search_range = search_range(image_y_idx_init+search_range >= 1 ...
    & image_y_idx_init+search_range < size(param.image_c,1));
  [~,max_idx_init] = max(param.image_c(image_y_idx_init+search_range,image_x_idx_init));
  max_idx_init = search_range(max_idx_init) + image_y_idx_init;
end

% [~,point_idxs_init] = min(abs(param.layer.x-x(1)));
new_y_init = interp1(1:length(image_y),image_y,max_idx_init,'linear','extrap');
pnt_init = [image_x_idx_init,image_y_idx_init];
param.val_init = param.image_c(pnt_init(2),pnt_init(1));

%% Find the closest point in the image - Final Point
[~,image_x_idx_final] = min(abs(image_x-x(2)));
dy = image_y(2)-image_y(1);
if stat_idx_y == 2
  image_y_idx_final = 1 + round((y(1)-image_y(1))/dy);
else
  image_y_idx_final = 1 + round((y(2)-image_y(1))/dy);
end
if image_y_idx_final < 1 || image_y_idx_final > size(param.image_c,1)
  max_idx_final = image_y_idx_final;
  image_y_idx_final = size(param.image_c,1);
else   
  search_range = search_range(image_y_idx_final+search_range >= 1 ...
    & image_y_idx_final+search_range < size(param.image_c,1));
  [~,max_idx_final] = max(param.image_c(image_y_idx_final+search_range,image_x_idx_final));
  max_idx_final = search_range(max_idx_final) + image_y_idx_final;
end

% [~,point_idxs_final] = min(abs(param.layer.x-x(2)));
new_y_final = interp1(1:length(image_y),image_y,max_idx_final,'linear','extrap');
pnt_final = [image_x_idx_final,image_y_idx_final];
param.val_final = param.image_c(pnt_final(2),pnt_final(1));

if image_y_idx_init > image_y_idx_final
  SelArea = param.image_c(image_y_idx_final:image_y_idx_init,image_x_idx_init:image_x_idx_final);
else
  SelArea = param.image_c(image_y_idx_init:image_y_idx_final,image_x_idx_init:image_x_idx_final);
end
rline = image_x_idx_init;
cur_frm = num2str(find(param.echowin.eg.image_gps_time(rline) >= param.echowin.eg.start_gps_time,1,'last'));
if isempty(cur_frm)
  cur_frm = num2str(1);
end
frm = strcat(param.echowin.eg.cur_sel.day_seg,'_',cur_frm);
ip = ['[',num2str(x(1)),',',num2str(new_y_init),']'];
fp = ['[',num2str(x(2)),',',num2str(new_y_final),']'];
if stat_idx == 1
  SelArea = SelArea(:);
  max_val = sprintf('%.2f',nanmax(SelArea));
  min_val = sprintf('%.2f',nanmin(SelArea));
  avg_val = sprintf('%.5f',10*log10(nanmean(10.^(SelArea./10))));
  med_val = sprintf('%.2f',nanmedian(SelArea));
  fprintf('Frame ID: %s\nSelected Region (X,Y): %s to %s\nMax Value: %s\nMin Value %s\nAverage Value: %s\nMedian Value: %s\n',frm,ip,fp,max_val,min_val,avg_val,med_val)
  
elseif stat_idx == 2
  h_fig = figure;
  clf(h_fig(1));
  h_axes(1) = axes('parent',h_fig(1));
  histogram(h_axes(1),SelArea)
  title(h_axes(1),strcat(frm,{' '},ip,' to',{' '},fp),'Interpreter','none');
  xlabel('Relative power (dB)')
  ylabel('Frequency')
  
elseif stat_idx == 3
  max_val = nanmax(SelArea);
  min_val = nanmin(SelArea);
  avg_val = 10*log10(nanmean(10.^(SelArea./10)));
  med_val = nanmedian(SelArea);
  h_fig = figure;
  clf(h_fig(1));
  h_axes(1) = axes('parent',h_fig(1));
  rline = image_x_idx_init:1:image_x_idx_final;
  plot(h_axes(1),rline,max_val);
  hold(h_axes(1),'on')
  plot(h_axes(1),rline,min_val);
  hold(h_axes(1),'on')
  plot(h_axes(1),rline,avg_val);
  hold(h_axes(1),'on')
  plot(h_axes(1),rline,med_val);
  hold(h_axes(1),'on')
  xlabel('Range line')
  ylabel('Relative power (dB)')
  title(h_axes(1),strcat(frm,{' '},ip,' to',{' '},fp),'Interpreter','none')
  legend('Max Value', 'Min Value', 'Average Value', 'Median Value')
 end

return

