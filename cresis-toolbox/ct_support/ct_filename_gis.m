function fn = ct_filename_gis(param,fn)
% fn = ct_filename_gis(param,fn)
%
% Returns a standardized filename for temporary files.
% 1. Handles absolute and relative path conversions, default paths, and
%    file separator differences (e.g. windows uses '\' and unix uses '/')
% 2. There are three modes of operation:
%  - base_fn is an absolute path: the param.tmp_path and standardized
%    path are not used (base_fn replaces both)
%  - base_fn is a relative path: the param.tmp_path is used, but the
%    standardized path/directory structure is not used (base_fn is used)
%  - base_fn is empty: the param.gis_path and standardized path are used
%
% param = control structure to data processor
% fn = parameter filename provided
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis

global gRadar;
param = merge_structs(gRadar,param);

if ~exist('fn','var')
  fn = [];
end

if ~isempty(fn) && (fn(1) == filesep || (ispc && (~isempty(strfind(fn,':\')) || ~isempty(strfind(fn,':/')))))
  % This is already an absolute path
  return
else
  % Append the current path to the support path
  fn = fullfile(param.gis_path, fn);
end

fn(fn == '/' | fn == '\') = filesep;

return;
