function [status,decoded_json] = json_response_decode(json_response)
%
% [status,decoded_json] = json_response_decode(json_response)
%
% Decodes a database json_response object. Formats all responses in
% fromjson() format.
%
% Input:
%   json_response: json string returned from django (the database)
%
% Output:
%   an error may be thrown, OR
%
%   decoded_json: a matlab structure decoded from json_response
%
% Author: Kyle W. Purdon


try
  % FIRST TRY TO LOAD WITH FROMJSON() 
  response = fromjson(json_response);
  status = response.status;
  message = response.data;
catch ME
  % IF FROMJSON() FAILS REVER TO LOADJSON()
  warning('fromjson() failed, loadjson() will be used.')
  response = loadjson(json_response);
  status = response.status;
  message = response.data;
end

% THROW ERROR IF GET FAILS
if status == 0
  error(message);
elseif status == 2
  decoded_json = message;
  warning(message);
  return;
else
  if ~strcmp(message(1),'{')
    % RETURN THE MESSAGE (THERE IS NO DATA)
    decoded_json = message;
  else
    try
      % CONVERT TO THE OUTPUT PARAM (NEW JSON METHOD)
      decoded_json = fromjson(message);
    catch ME
      % CONVERT TO THE OUTPUT PARAM (OLD JSON METHOD)
      lj = loadjson(message);
      lj_fns = fieldnames(lj);
      for fidx = 1:length(lj_fns)
        data = getfield(lj,lj_fns{fidx});
        if isnumeric(data)
          % CONVERT NUMERIC ARRAY TO CELL ARRAY
          lj = setfield(lj,lj_fns{fidx},num2cell(data)');
        elseif ischar(data)
          if length(data) > 1
            % CONVER CHARACTER ARRAY TO CELL ARRAY
            lj = setfield(lj,lj_fns{fidx},cellstr(data));
          else
            % CONVERT STRING TO CELL
            lj = setfield(lj,lj_fns{fidx},{data});
          end
        elseif iscell(data)
          % FLIP CELL TO CORRECT ORIENTATION
          lj = setfield(lj,lj_fns{fidx},data');
        else
          warning('unkown format');
          continue;
        end
      end
      decoded_json = lj;
    end
  end
end
end