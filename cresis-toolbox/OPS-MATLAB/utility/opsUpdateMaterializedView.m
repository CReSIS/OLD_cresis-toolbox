function [status,message] = opsUpdateMaterializedView(update)

% REFRESH MATERIALIZED VIEW xxx; on the database.
% 	
% 	Input:
% 		update: (string) an SQL clause 'REFRESH MATERIALIZED VIEW xxx'
% 		
% 	Output:
% 		status: (integer) 0:error 1:success
		
% Author: Weibo Liu

% INITIATE THE DATABASE TRANSACTION
opsCmd;
if gOps.profileCmd
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'profile'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'query' update 'view' 'updateMaterializedView'});
else
    [jsonResponse,~] = opsUrlRead(strcat(gOps.serverUrl,'updateMaterializedView'),gOps.dbUser,gOps.dbPswd,...
        'Post',{'query' update});
end

% DECODE THE SERVER RESPONSE
[status,decodedJson] = jsonResponseDecode(jsonResponse);

% CREATE THE DATA OUPUT STRUCTURE OR MESSAGE
message = decodedJson;
return