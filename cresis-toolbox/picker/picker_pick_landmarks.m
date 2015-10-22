function picker_pick_landmarks
% % picker_pick_landmarks
% %
% % Support function for the creation of landmark handles for each individual
% % landmark. This function also supports handling the handles for switching
% % frames.
% %
% % Author: Aric Beaver
% 
% 
% % Create a new function which handles switching frames: this function will 
% % delete all the landmark handles from the last frame and create new handles 
% % for the currently loaded frame
% 

global gCtrl;
global hui;

if isempty(gCtrl.source.landmarks.fn) || ~exist(gCtrl.source.landmarks.fn,'file')
  landmarks = [];
else
  load(gCtrl.source.landmarks.fn)
end

% ===================================================================
% Create landmark handles
% ===================================================================
% Find applicable landmarks for current frame based on LB menuString

if ~isempty(landmarks)
  lm_frm_idx = [];
  pick_frm_id = gCtrl.source.frm_id(gCtrl.source.cur_pick,1:end); 
  cur_frm_lm_idxs = [];
  for lm_idx = 1:length(landmarks)
    lm_frm_id = char(landmarks(lm_idx).frm_id);
    if strcmp(lm_frm_id,pick_frm_id)
      if strcmp(landmarks(lm_idx).season_name,gCtrl.season_name)
        cur_frm_lm_idxs(end+1) = lm_idx; 
      end  
    end
  end

  % Get x/y data for plotting landmarks
  if ~isempty(cur_frm_lm_idxs)
    for lm_idx = 1:length(cur_frm_lm_idxs)

      tmp_x1 = interp1(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time,...
        1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time),landmarks(cur_frm_lm_idxs(lm_idx)).gpstime_start,'linear','extrap');
      tmp_x2 = interp1(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time, ...
        1:length(gCtrl.source.geo{gCtrl.source.cur_pick}.GPS_time), landmarks(cur_frm_lm_idxs(lm_idx)).gpstime_stop,'linear','extrap');
      
%       x1(lm_idx) = landmarks(cur_frm_lm_idxs(lm_idx)).rline_start;
%       x2(lm_idx) = landmarks(cur_frm_lm_idxs(lm_idx)).rline_stop;
      
      x1(lm_idx) = tmp_x1;
      x2(lm_idx) = tmp_x2;
      y1(lm_idx) = landmarks(cur_frm_lm_idxs(lm_idx)).rbin_start;
      y2(lm_idx) = landmarks(cur_frm_lm_idxs(lm_idx)).rbin_stop;
    end
    figure(hui.pickfig.handle)
    hold on;
    hui.pickfig.landmarks_h = [];
    for lm_idx = 1:length(cur_frm_lm_idxs)
      hui.pickfig.landmarks_h(cur_frm_lm_idxs(lm_idx)) = plot([x1(lm_idx) x1(lm_idx) x2(lm_idx) x2(lm_idx) x1(lm_idx)],...
        [y1(lm_idx) y2(lm_idx) y2(lm_idx) y1(lm_idx) y1(lm_idx)],'Color','y','LineStyle','--','Visible','on');
    end
    hold off;
  else
    % plot nothing
  end
else
  hui.pickfig.landmarks_h = [];
end

return
  
  
  
