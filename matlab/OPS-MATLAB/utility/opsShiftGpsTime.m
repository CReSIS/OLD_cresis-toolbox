function [ status,message ] = opsShiftGpsTime( sys,param )
%
% [status,data] = opsShiftGpsTime(sys,param)
%
% Shift the GPS time of layer points in the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.offset = number
%     properties.location = string ('arctic' or 'antarctic')
%     properties.season = string
%     properties.segment = string
%     properties.lyr_name = string ('surface','bottom', etc...)
%
% Output:
%   data: see opsGetLayerPoints
%   message: status message
%
% Author: Weibo Liu & Kyle W. Purdon 

% authenticate the user
[authParam,~,~] = opsAuthenticate(struct('properties',[]));

% get the layer points
[~,data] = opsGetLayerPoints(sys,param);

% shift the gps time
shifted_gps_time = data.properties.gps_time + param.properties.offset;

% interpolate twtt onto the original gps time
if sum(isfinite(data.properties.twtt)) >= 2
  out_twtt = interp1(shifted_gps_time,data.properties.twtt,data.properties.gps_time);
else
  status = 1;
  msg = 'All points NaN, skipping shift';
  fprintf('Gps Time Shift Successful:%d: %s\n', status, msg);
  message = sprintf('Gps Time Shift Successful:%d: %s', status, msg);
  return
end

% construct the param for point creation
createParam.properties.point_path_id = data.properties.point_path_id;
createParam.properties.username = authParam.properties.userName;
createParam.properties.twtt = out_twtt;
createParam.properties.type = data.properties.type;
createParam.properties.quality = data.properties.quality;
if ~isnan(param.properties.quality)
  % Overriding quality setting
  createParam.properties.quality(:) = param.properties.quality;
end
createParam.properties.lyr_name = param.properties.lyr_name;

% create the layer points
try
  [status,msg] = opsCreateLayerPoints(sys,createParam);
catch ME
   for idx=1:length(ME.stack); disp(ME.stack(idx)); end;
   error('opsCreateLayerPoints failed: %s %s', ME.message);
end
 
if status == 1
    fprintf('Gps Time Shift Successful:%d: %s\n', status, msg);
    message = sprintf('Gps Time Shift Successful:%d: %s', status, msg);
else
    warning('Gps Time Shift Failed:%d: %s', status, msg);
    message = sprintf('Gps Time Shift Failed:%d: %s', status, msg);
end

end