function customverctrl(fileNames, arguments)

action = arguments{strcmp('action',arguments(:,1)),2};

if ~strcmp(action,'checkin')
  return;
end

sys_cmd = sprintf('svn status %s', fileNames{1})
fprintf('  running %s\n', sys_cmd);
[status,result] = system(sys_cmd);

comments = arguments{strcmp('comments',arguments(:,1)),2};
if result(1) == '?'
  sys_cmd = sprintf('svn add %s', fileNames{1})
  fprintf('  running %s\n', sys_cmd);
  [status,result] = system(sys_cmd);
  sys_cmd = sprintf('svn ci --message "%s" %s', comments, fileNames{1})
  fprintf('  running %s\n', sys_cmd);
  system(sys_cmd);
elseif result(1) == 'A' || result(1) == 'M'
  sys_cmd = sprintf('svn ci --message "%s" %s', comments, fileNames{1})
  fprintf('  running %s\n', sys_cmd);
  system(sys_cmd);
else
  fprintf('Does not need to be checked in:\n  %s', result);
end

return

