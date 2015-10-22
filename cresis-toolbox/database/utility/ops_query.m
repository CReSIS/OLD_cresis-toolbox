function [status,data] = ops_query(query)
%
% [status,data] = ops_query(query)
%
% Retrieves query results from the database.
%
% Input:
%   query: string SQL (SELECT) query
%
% Output:
%   data: cell array with varying fields
%     Converted JSON query results.
%
% Example:
%   query = 'SELECT * FROM rds_seasons WHERE season_name="2011_Greenland_P3";';
%   [status,data] = ops_query(query);
%   season_id = data{1};
%   season_name = data{2};
%   descriptions = data{3};
%
% WARNING:
%     If your query returns nothing an error will be thrown even if
%   the query worked. Always use the RETURN statement in your SQL to confirm
%   a succesfull query.
%
% Author: Kyle W. Purdon

% INITIATE THE DATABASE TRANSACTION
ops_sys_cmd;
if profile_cmd
  [json_response,~] = cr_urlread(strcat(server_url,'profile'),db_user,db_pswd,...
    'Post',{'query' query 'view' 'query'});
else
  [json_response,~] = cr_urlread(strcat(server_url,'query'),db_user,db_pswd,...
    'Post',{'query' query});
end

try
  % DECODE THE QUERY RESPONSE (NEW JSON METHOD)
  response = fromjson(json_response);
  status = response.status;
  % THROW ERROR IF QUERY FAILS
  if status == 0
    error(response.data{1});
  end
  if status == 2
    warning(response.data{1});
  end
  for idx=1:length(response.data)
    for idx2=1:length(response.data{1})
      data{idx2,idx} = response.data{idx}{idx2};
    end
  end
catch ME
  % DECODE THE QUERY RESPONSE (OLD JSON METHOD)
  response = loadjson(json_response);
  status = response.status;
  % THROW ERROR IF QUERY FAILS
  if status == 0
    error(response.data);
  end
  if status == 2
    warning(response.data);
  end
  data = json_wrapper(response.data);
end
return