% script runOpsGetFramesWithinPolygon
%
% Example script for running opsGetFramesWithinPolygon
%

%% User Settings
sys = 'rds';

param = [];
param.properties.location = 'arctic';
param.properties.bound = 'POLYGON((-41.385225931851004 74.18070913812296,-40.953933829180954 74.21354044957413,-40.7835769291146 74.10614286784377,-41.261468378407606 74.07456002598235,-41.32999713155497 74.0983784767899,-41.385225931851004 74.18070913812296))';

%% Query database
[status,message] = opsGetFramesWithinPolygon(sys,param);

%% Interpret the return
if status == 1
  % Create a list of unique segments
  for frm_idx = 1:length(message.frame)
    % Get the segment name portion of the frame name
    seg{frm_idx} = message.frame{frm_idx}(1:11);
  end
  seg = unique(seg);
  
  % Print list of segments
  for seg_idx = 1:length(seg)
    fprintf('%s\n', seg{seg_idx});
  end
end
