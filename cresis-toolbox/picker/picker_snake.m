function vals = picker_snake(indices,values)
% vals = picker_snake(indices,values)
%
% Snake layer tracking tool called from picker_pick_button.
%
% Author: John Paden

global gCtrl;
global hui;

% Interpolate data onto main xaxis
A = get(hui.pickfig.image.h,'CData');
%A = interp1(gCtrl.pick.xaxis,A.',1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time)).';

% Only snake on the selection
A = A(:,indices(1):indices(end));

tool_param1_str = get(hui.fig.ctrl_panel.tool_param2_TE,'String');
search_range = [];
try
  % Assumes the ylimit_str is a matlab expression that can be evaluated
  search_range = eval(sprintf('[%s]', tool_param1_str));
  if length(search_range) == 1
    search_range = -search_range:search_range;
  else
    search_range = round(min(search_range)):round(max(search_range));
  end
end
if isempty(search_range)
  search_range = -3:3;
end

for idx = 1:length(indices)
  dataPnts(idx).row = round(interp1(gCtrl.pick.time*1e6, ...
    1:length(gCtrl.pick.time),values(idx),'linear','extrap'));
  dataPnts(idx).col = indices(idx) - indices(1) + 1;
  dataPnts(idx).method = 's';
  dataPnts(idx).snake.search_range = search_range;
end

layer = tracker_snake(A,dataPnts);

vals = interp1(1:length(gCtrl.pick.time),gCtrl.pick.time*1e6,layer,'linear','extrap');

return;
