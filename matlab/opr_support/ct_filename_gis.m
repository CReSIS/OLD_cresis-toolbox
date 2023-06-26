function fn = ct_filename_gis(param,fn)
% fn = ct_filename_gis(param,fn)
%
% Returns a standardized filename for GIS files.
% 1. Handles absolute and relative path conversions, default paths, and
%    file separator differences (e.g. windows uses '\' and unix uses '/')
% 2. There are two modes of operation:
%  - fn is an absolute path: the param.gis_path is not used.
%  - fn is a relative path: fn is appended to param.gis_path
%
% param: control structure from parameter spreadsheet, if empty this
% function will try to use the global variable gRadar
%  .gis_path: Path to GIS files
% fn: GIS filename (either absolute path or relative to param.gis_path)
%
% Legacy format:
% param: The "fn" from above.
% fn: NOT USED
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis, ct_filename_param

if ischar(param)
  % Legacy format
  fn = param;
  param = [];
end

global gRadar;
param = merge_structs(gRadar,param);

if ~exist('fn','var')
  fn = [];
end

if ~isempty(fn) && (fn(1) == filesep || (ispc && (~isempty(strfind(fn,':\')) || ~isempty(strfind(fn,':/')))))
  % This is already an absolute path
  return
else
  % Append the current path to the GIS path
  if ~isfield(param,'gis_path')
    error('gis_path is missing from global variable gRadar');
  end
  fn = fullfile(param.gis_path, fn);
end

fn(fn == '/' | fn == '\') = filesep;

return;
