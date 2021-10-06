% script runOpsGetFramesWithinPolygon
%
% Example script for running opsGetFramesWithinPolygon
%

%% User Settings
sys = 'rds';

param = [];
param.properties.location = 'arctic';
param.properties.bound = 'POLYGON((-46.348462741297276 67.1603370256662,-49.34129895497196 66.14631499808071,-45.32811799251264 65.89532478450859,-46.348462741297276 67.1603370256662))';

%% Query database
[status,message] = opsGetFramesWithinPolygon(sys,param);

%% Interpret the return
if status == 1
  % Create a list of unique segments
  for frm_idx = 1:length(message.frame)
    % Get the segment name portion of the frame name
    seg{frm_idx} = message.frame{frm_idx}(1:11);
  end
  [seg,idx1,idx2] = unique(seg);
  
  % Print list of segments
  fprintf('day_seg\tnum frm\tfrms\n');
  for seg_idx = 1:length(seg)
    % Print segment ID (day_seg)
    fprintf('%s\t', seg{seg_idx});
    % Get frame IDs associated with this segment
    frm_ids = message.frame(idx2==seg_idx);
    % Print number of frames
    fprintf('%d\t', length(frm_ids));
    % Print list of frames in nice matlab array format
    %   Instead of [1,2,3,6], this prints [1:3,6]
    fprintf('[');
    state = 0; % State machine for printing
    for frm_ids_idx = 1:length(frm_ids)
      frm = str2double(frm_ids{frm_ids_idx}(13:end));
      switch state
        case 0 % Initial state
          state = 1;
          fprintf('%d',frm);
        case 1 % New sequential range started
          if frm == old_frm + 1
            state = 2;
          else
            fprintf(',%d',frm);
          end
        case 2 % Within a sequential range
          if frm ~= old_frm + 1
            fprintf(':%d',old_frm);
            fprintf(',%d',frm);
            state = 1;
          end
      end
      old_frm = frm;
    end
    if state == 2 % Print off the last unfinished range
      fprintf(':%d',frm);
    end
    fprintf(']\n');
  end
end
