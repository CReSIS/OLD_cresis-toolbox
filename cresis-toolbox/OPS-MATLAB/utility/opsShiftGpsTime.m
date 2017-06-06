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
out_twtt = interp1(shifted_gps_time,data.properties.twtt,data.properties.gps_time);

% construct the param for point delete
 deleteParam.properties.start_point_path_id = min (data.properties.point_path_id);
 deleteParam.properties.stop_point_path_id = max (data.properties.point_path_id);
 deleteParam.properties.max_twtt = max (data.properties.twtt);
 deleteParam.properties.min_twtt = min (data.properties.twtt);
 deleteParam.properties.lyr_name = param.properties.lyr_name;
 
 % delete original layer points
 try
  [status,msg] = opsDeleteLayerPoints(sys,deleteParam);
 catch ME
   error('Something wrong with point delete: %d:%s',status, msg);
 end

% construct the param for point creation
idx = ~isnan(out_twtt);
createParam.properties.point_path_id = data.properties.point_path_id(idx);
createParam.properties.username = authParam.properties.userName;
createParam.properties.twtt = out_twtt(idx);
createParam.properties.type = data.properties.type(idx); % need to verify
createParam.properties.quality = data.properties.quality(idx); % need to verify
createParam.properties.lyr_name = param.properties.lyr_name;

% create the layer points
try
  [status,msg] = opsCreateLayerPoints(sys,createParam);
catch ME
   error('Something wrong with layer point creation: %d:%s', status, msg);
end
 
if status == 1
    fprintf('Gps Time Shift Successful\n');
    message = 'Gps Time Shift Successful';
else
    fprintf('Gps Time Shift Failed\n');
    message = 'Gps Time Shift Failed';
end

end