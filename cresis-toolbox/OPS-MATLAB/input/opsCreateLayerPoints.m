function [status,message] = opsCreateLayerPoints(sys,param)
%
% [status,message] = opsCreateLayerPoints(sys,param)
%
% Loads a layer into the database.
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   param: structure with fields
%     properties.point_path_id = integer array 
%     properties.twtt = double array
%     properties.type = integer arry (0,1 or 2)
%     properties.quality = integer array (1,2 or 3)
%     properties.lyr_name OR properties.lyr_id
%       lyr_name: string ('surface','bottom', etc...)
%       lyr_id = scalar integer (database ID)
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%   message: status message
%
% Note a NaN in twtt will force the layer point to be deleted in the database.
%
% Author: Kyle W. Purdon

% SET UP DEFAULT VALUES FOR GROUP
if ~isfield(param.properties,'lyr_group_name')
  if isfield(param.properties,'lyr_name')
    if any(strcmp(param.properties.lyr_name,{'surface','bottom'}))
      param.properties.lyr_group_name = 'standard';
    else
      param.properties.lyr_group_name = 'uncatagorized';
    end
  else
    if any(param.properties.lyr_id == [1 2]) % surface/bottom
      param.properties.lyr_group_name = 'standard';
    else
      param.properties.lyr_group_name = 'uncatagorized';
    end
  end
end

% CONSTRUCT THE JSON STRUCTURE
[param,~,~] = opsAuthenticate(param);
param.properties.twtt(~isfinite(param.properties.twtt)) = NaN;
jsonStruct = struct('properties',param.properties);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
  jsonStr = tojson(jsonStruct);
catch ME
  jsonStr = savejson('',jsonStruct,'FloatFormat','%2.10f','NaN','null');
end

% SEND THE COMMAND TO THE SERVER
opsCmd;
if gOps.profileCmd
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr 'view' 'createLayerPoints'});
else
  [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'create/layer/points'),gOps.dbUser,gOps.dbPswd,...
    'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

end