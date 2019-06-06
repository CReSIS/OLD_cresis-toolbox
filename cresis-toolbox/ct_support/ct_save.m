function ct_save(fn,varargin)

global gRadar;
if ~isfield(gRadar,'min_disk_space') || isempty(gRadar.min_disk_space)
  min_disk_space = 10e9;
else
  min_disk_space = gRadar.min_disk_space;
end

free = get_disk_space(fn);

if free < min_disk_space
  error('Insufficient disk space (%g MB free, minimum allowed % MB).', free/1e6, min_disk_space/1e6);
end

save(fn,varargin{:});
