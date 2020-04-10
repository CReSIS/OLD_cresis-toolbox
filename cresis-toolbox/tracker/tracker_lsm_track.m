function labels = tracker_lsm_track(data,param)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Load one frame at a time
frm = param.layer_tracker.tracker.frm;
   lsm_tic = tic;

%   if isfield(param.layer_tracker.track.lsm, 'surf_layer_params') && ~isempty(param.layer_tracker.track.lsm.surf_layer_params)
%     surf_layer_params = param.layer_tracker.track.lsm.surf_layer_params;
%   else
%     surf_layer_params = [];
%   end
%
%   if isfield(param.layer_tracker.track.lsm, 'bot_layer_params') && ~isempty(param.layer_tracker.track.lsm.bot_layer_params)
%     bot_layer_params = param.layer_tracker.track.lsm.bot_layer_params;
%   else
%     bot_layer_params = [];
%   end
%
%   try
%     Surface = opsLoadLayers(param,surf_layer_params);
%     Bottom = opsLoadLayers(param,bot_layer_params);
%   catch ME
%     warning(ME.getReport);
%     keyboard
%   end

obj          = tomo.LSMObject_tuning({lp(data.Data)});
obj.setLSMOptions('y', param.layer_tracker.track.lsm.y, 'dy', param.layer_tracker.track.lsm.dy, 'outerIter', param.layer_tracker.track.lsm.numOuterIter);

[flag, Labels.top, Labels.bot, matrix_x, matrix_y] = obj.runLSM();
temp = size(matrix_y,3);
for i =1:2
  for idx = 1:temp
    try
      matrix_y(i,:,idx) = interp1(matrix_x(i,:,idx),matrix_y(i,:,idx),1:length(data.Bottom),'linear','extrap');
    catch ME
      flag(idx) = 0;
      continue
    end
  end
end


%   lsm_toc = toc(lsm_tic);
%   figure; imagesc(lp(data.Data)); colormap(1 - gray(256)); hold on;
%   plot(matrix_y(1,:,temp), 'g'); plot(matrix_y(2,:,temp), 'r');
%   legend('Ice-surface', 'Ice-bottom');
%   keyboard
  
  labels = matrix_y;

