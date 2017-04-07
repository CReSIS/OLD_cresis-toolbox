function sw_version = current_software_version
% sw_version = current_software_version
%
% Returns current version of the toolbox

sw_version.ver = '2.0';
sw_version.date_time = '20130102_120000';
sw_version.cur_date_time = datestr(now);
try
  if 0
    % Get the SVN version (requires a command line svn program to be in the
    % system path)
    global gRadar;
    [status,result] = system(sprintf('svn info %s',gRadar.path));
    sw_version.rev = str2double(strtok(result(regexp(result,'Revision: ')+10:end)));
    sw_version.URL = strtok(result(regexp(result,'URL: ')+5:end));
    
  else
    % Get the Git version (requires a command line git program to be in the
    % system path)
    global gRadar;
    % Ensure gRadar.path has no file separators at the end so that fileparts
    % returns the parent directory which should be the repository
    % directory. The extra "cat" at the end was done to fix an odd bug in
    % the matlab/git interface that manifested itself once.
    [status,result]=system(sprintf('export TERM=linux; git --git-dir="%s" log -1 | cat', ...
      fullfile(fileparts(gRadar.path(1:find(gRadar.path~='/' & gRadar.path~='\',1,'last'))),'.git')));
    sw_version.rev = strtok(result(regexp(result,'commit ')+7:end));
    sw_version.URL = '';
  end
catch
  sw_version.rev = NaN;
  sw_version.URL = NaN;
end

return;
