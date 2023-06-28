function fn = ct_filename_param(param,fn)
% fn = ct_filename_param(param,fn)
%
% Returns a standardized filename for param spreadsheet files.
% 1. Handles absolute and relative path conversions, default paths, and
%    file separator differences (e.g. windows uses '\' and unix uses '/')
% 2. There are two modes of operation:
%  - fn is an absolute path: the param.param_path is not used.
%  - fn is a relative path: fn is appended to param.param_path
%
% param: control structure from parameter spreadsheet, if empty this
% function will try to use the global variable gRadar
%  .param_path: Path to parameter spreadsheet files
% fn: parameter spreadsheet filename (either absolute path or relative
%   to param.param_path)
%
% Legacy format:
% param: The "fn" from above.
% fn: NOT USED
%
% Author: Kyle Purdon, John Paden
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
  % Append the current path to the param path
  if ~isfield(param,'param_path')
    error('param_path is missing from global variable gRadar');
  end
  fn = fullfile(param.param_path, fn);
end

fn(fn == '/' | fn == '\') = filesep;

return;
