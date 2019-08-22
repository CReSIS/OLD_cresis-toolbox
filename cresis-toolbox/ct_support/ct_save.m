function ct_save(fn,varargin)

global gRadar;
if ~isfield(gRadar,'min_disk_space') || isempty(gRadar.min_disk_space)
  min_disk_space = 10e9;
else
  min_disk_space = gRadar.min_disk_space;
end

free = get_disk_space(fn);

% Assume free == 0 is an error in get_disk_space and so ignore it because
% the save will fail anyway if this is the case.
if free > 0 && free < min_disk_space
  error('Insufficient disk space (%g MB free, minimum allowed % MB).', free/1e6, min_disk_space/1e6);
end

fn_dir = fileparts(fn);
if ~exist(fn_dir,'dir')
  mkdir(fn_dir);
end
cmd = sprintf('save(''%s''%s)',fn,sprintf(',''%s''',varargin{:}));
evalin('caller',cmd);
