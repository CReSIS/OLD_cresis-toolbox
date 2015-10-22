function [status,message] = opsShiftTwtt(sys,param)
%
% [status,data] = opsShiftTwtt(sys,param)
%
% Shift the TWTT of layer points in the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.location = string ('arctic' or 'antarctic')
%     properties.season = string
%     properties.segment = string
%     properties.lyr_name = string ('surface','bottom', etc...)
%     properties.offset = double (twtt offset)
%
% Output:
%   status: 0 (error), 1 (success), or 2 (warning)
%   message: status message
%
% Author: Kyle W. Purdon, Weibo Liu, Trey Stafford

% authenticate the user
[authParam,~,~] = opsAuthenticate(struct('properties',[]));

% Create SQL query to shift twtt
query = sprintf('UPDATE %s_layer_points SET twtt = twtt + %s WHERE layer_id = (SELECT id FROM %s_layers WHERE name = ''%s'') AND point_path_id IN (SELECT pp.id FROM %s_point_paths pp JOIN %s_segments seg ON pp.segment_id=seg.id WHERE seg.name= ''%s'') RETURNING 1;', sys,param.properties.offset,sys,param.properties.lyr_name,sys,sys,param.properties.segment);
% Query the DB to update layer points with twtt shift
[status,data] = opsQuery(query);

% Check status and return message
if status == 1
    fprintf('TWTT Shift Successful\n');
    message = 'TWTT Shift Successful';
else
    fprintf('TWTT Shift Failed\n');
    message = 'TWTT Shift Failed';
end

end