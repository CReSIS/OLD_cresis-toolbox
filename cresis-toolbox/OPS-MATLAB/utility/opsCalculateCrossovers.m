function [status,data] = opsCalculateCrossovers()
%
% [status,data] = opsCalculateCrossovers()
%
% Retrieves segments which have not had crossovers calculated and
% calculates crossovers on them.
%
% Input:
%   none
%
% Output:
%   status: integer (0:Error,1:Success,2:Warning)
%
% Author: Reece Mathews

[authParam,~,~] = opsAuthenticate(struct('properties',[]));

% SEND THE COMMAND TO THE SERVER
opsCmd;
sys = 'rds';

query = sprintf('SELECT id from %s_segments where crossover_calc=false;', sys);
[~,data] = opsQuery(query);
segments = [data{:}];
jsonStruct = struct('segments',segments);

% CONVERT THE JSON STRUCTURE TO A JSON STRING
try
    jsonStr = tojson(jsonStruct);
catch ME
    jsonStr = savejson('',jsonStruct,'FloatFormat','%d','NaN','null');
end

if ~isinteger(segments(1))
  disp 'No crossovers to calculate';
  return
end

% SEND THE COMMAND TO THE SERVER
opsCmd;
if gOps.profileCmd
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr 'view' 'calculateCrossovers'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'calculateCrossovers'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'app' sys 'data' jsonStr});
end

% DECODE THE SERVER RESPONSE
[status,message] = jsonResponseDecode(jsonResponse);

end
