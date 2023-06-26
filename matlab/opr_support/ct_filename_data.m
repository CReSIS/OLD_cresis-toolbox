function fns = ct_filename_data(param,base_fn,fns)
% fns = ct_filename_data(param,base_fn,fns)
%
% Returns a standardized filename for data files or if fns is not declared
% then the standardized path to the filenames is given.
% 1. Handles absolute and relative path conversions, default paths, and
%    file separator differences (e.g. windows uses '\' and unix uses '/')
% 2. There are three modes of operation:
%  - base_fn is an absolute path: the param.data_path and standardized
%    path are not used (base_fn replaces both)
%  - base_fn is a relative path: the param.data_path is used, but the
%    standardized path/directory structure is not used (base_fn is used)
%  - base_fn is empty: the param.data_path and standardized path are used
%
% param = control structure to data processor
%  .radar_name (e.g. mcords)
%  .season_name (e.g. 2009_antarctica_DC8)
%  .day_seg (e.g. 20091020_03)
% base_fn = parameter filename provided
% fns = cell array of data files (optional)
%   if fns is not specified, then just the path to the data files
%   is returned
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis

global gRadar;
param = merge_structs(gRadar,param);

if ~exist('fns','var')
  not_a_cell = true;
  fns{1} = '';
else
  not_a_cell = false;
end

[output_dir,radar_type] = ct_output_dir(param.radar_name);

for file_idx = 1:length(fns)
  if isempty(base_fn)
    % Generate the default path
    fns{file_idx} = fullfile(param.data_path, output_dir, param.season_name, ...
      param.day_seg(1:8), fns{file_idx});
  elseif base_fn(1) == filesep || (ispc && ~isempty(strfind(base_fn,':')))
    % This is already an absolute path
    fns{file_idx} = fullfile(base_fn, fns{file_idx});
    continue;
  else
    % Append the current path to the support path
    fns{file_idx} = fullfile(param.data_path, base_fn, fns{file_idx});
    fns{file_idx}(fns{file_idx} == '/' | fns{file_idx} == '\') = filesep;
  end
end

if not_a_cell
  fns = fns{1};
end

return;
