global gRadar
sys_cmd = sprintf('svn update %s', gRadar.path);
fprintf('  running %s\n', sys_cmd);
system(sys_cmd);

return;