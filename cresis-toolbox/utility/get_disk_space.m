function free = get_disk_space(path)
% free = get_disk_space(path)
%
% Returns the amount of free space left on the drive that path points to

if ~exist(path,'dir')
  [path] = fileparts(path);
end
if nargin < 1 || isempty(path)
  path= '.';
end

free = java.io.File(path).getFreeSpace();

