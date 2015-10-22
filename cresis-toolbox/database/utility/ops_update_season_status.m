% OPS UPDATE SEASONS STATUS
%
% Modifies the status for a season in the database (public,private)
%
% Input:
%   sys: (string) sys name ('rds','accum','snow',...)
%   season/s: (cell) season names ({'2011_Greenland_P3','2013_Greenland_P3'})
%   status: (integer) 0/1 (0:public,1:private)
%
% Output:
%   none
%
% Author: Kyle W. Purdon

sys = 'rds';
season = {'2011_Greenland_P3','2013_Greenland_P3'};
status = 0; %(0:public,1:private)

% AUTOMATED SECTION BELOW THIS LINE

if status == 0
  statusStr = 'public';
elseif status == 1
  statusStr = 'private';
else
  error('status must be 0 or 1');
end

seasonStr = '';
for sidx = 1:length(season)
  if sidx == 1
    seasonStr = cat(2,seasonStr,strcat('(','''',season{sidx},''''));
  elseif sidx == length(season) && length(season) > 2
    seasonStr = cat(2,seasonStr,strcat('''',season{sidx},''''));
  else
    seasonStr = cat(2,seasonStr,strcat(',','''',season{sidx},'''',')'));
  end
end
  
query = sprintf('UPDATE %s_seasons SET status=''%s'' WHERE season_id IN (SELECT season_id FROM %s_seasons WHERE season_name IN %s) RETURNING season_name,status;',sys,statusStr,sys,seasonStr);
[s,d] = ops_query(query);

fprintf('STATUS: %d\n',s);
d