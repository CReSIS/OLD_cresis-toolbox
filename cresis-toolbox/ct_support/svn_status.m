global gRadar
sys_cmd = sprintf('svn status -u %s', gRadar.path);
fprintf('  running %s\n', sys_cmd);
system(sys_cmd);

return;