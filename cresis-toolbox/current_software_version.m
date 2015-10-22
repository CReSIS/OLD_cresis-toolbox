function sw_version = current_software_version
% sw_version = current_software_version
%
% Returns current version of the toolbox

sw_version.ver = '2.0';
sw_version.date_time = '20130102_120000';
sw_version.cur_date_time = datestr(now);
try
  % Get the SVN version (requires a command line svn program to be in the
  % system path)
  global gRadar;
  [status,result] = system(sprintf('svn info %s',gRadar.path));
  sw_version.rev = str2double(strtok(result(regexp(result,'Revision: ')+10:end)));
  sw_version.URL = strtok(result(regexp(result,'URL: ')+5:end));
catch
  sw_version.rev = NaN;
  sw_version.URL = NaN;
end

return;
