function svn(varargin)

input_string = '';
for idx = 1:length(varargin)
  input_string = cat(2,input_string,' ',varargin{idx});
end
sys_cmd = sprintf('svn%s', input_string);
fprintf('  running %s\n', sys_cmd);
system(sys_cmd);

return
