function [pathstr, name, ext] = ct_fileparts(file)

file_ispc = ~isempty(regexpi(file,'[\:\\]'));

if ispc && file_ispc || ~ispc && ~file_ispc
  [pathstr, name, ext] = fileparts(file);
  return
elseif ispc && ~file_ispc
  fn_parts = regexp(file,'\/');
else
  fn_parts = regexp(file,'\\');
end

% Cut out directory
if ~isempty(fn_parts)
  pathstr = file(1:fn_parts(end)-1);
  file = file(fn_parts(end)+1:end);
else
  pathstr = '';
end

% Split name and extension
period_idx = find(file=='.',1,'last');
if ~isempty(period_idx)
  ext = file(period_idx:end);
  name = file(1:period_idx-1);
else
  ext = '';
  name = file;
end
