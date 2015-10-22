function fn = ct_filename_param(fn)
% fn = ct_filename_param(fn)
%
% Returns the standard filename for a param file based on gRadar.param_path
%
% fn = parameter filename provided (ex.
%
% Author: Kyle Purdon

global gRadar;

if ~exist('fn','var')
  fn = [];
end

if ~isempty(fn) && (fn(1) == filesep || (ispc && (~isempty(strfind(fn,':\')) || ~isempty(strfind(fn,':/')))))
  % This is already an absolute path
  return
else
  % Append the current path to the support path
  fn = fullfile(gRadar.param_path, fn);
end

fn(fn == '/' | fn == '\') = filesep;

return;
