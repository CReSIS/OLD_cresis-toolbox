function vals = picker_threshold(indices,values,first_bin,last_bin)
% vals = picker_threshold(indices,values,first_bin,last_bin)
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

tool_param1_str = get(hui.fig.ctrl_panel.tool_param1_TE,'String');
threshold = [];
num_lower = [];
try
  % Assumes the ylimit_str is a matlab expression that can be evaluated
  param_field = eval(sprintf('[%s]', tool_param1_str));
  threshold = param_field(1);
  num_lower = param_field(2);
end
if isempty(threshold)
  threshold = 10;
end
if isempty(num_lower)
  num_lower = 10;
end

noise_pow = mean(mean(A(first_bin:first_bin+5,:)));


for idx = 1:length(indices)
  dataPnts(idx).row = round(interp1(gCtrl.pick.time*1e6, ...
    1:length(gCtrl.pick.time),values(idx),'linear','extrap'));
  dataPnts(idx).col = indices(idx) - indices(1) + 1;
  dataPnts(idx).method = 't';
  dataPnts(idx).thresh.value = noise_pow + threshold;
  dataPnts(idx).thresh.num_lower = num_lower;
  dataPnts(idx).snake.search_range = -20:20;
end

layer = track_layer(A,dataPnts);

vals = interp1(1:length(gCtrl.pick.time),gCtrl.pick.time*1e6,layer,'linear','extrap');

return;
