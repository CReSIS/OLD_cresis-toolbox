function fn = ct_filename_out(param,fn,type,generic_data_flag)
% fn = ct_filename_out(param,fn,type,generic_data_flag)
%
% Returns a standardized filename for output files.
% 1. Handles absolute and relative path conversions, default paths, and
%    file separator differences (e.g. windows uses '\' and unix uses '/')
% 2. There are three modes of operation:
%  - base_fn is an absolute path: the param.out_path and standardized
%    path are not used (base_fn replaces both)
%  - base_fn is a relative path: the param.out_path is used, but the
%    standardized path/directory structure is not used (base_fn is used)
%  - base_fn is empty: the param.out_path and standardized path are used
%
% param = control structure to data processor
%  .radar_name (e.g. mcords)
%  .season_name (e.g. 2009_antarctica_DC8)
%  .day_seg (e.g. 20091020_03)
% fn = parameter filename provided
% type = string describing type of output (really just appended to the
%   output path), e.g. 'CSARP_standard, 'CSARP_qlook', etc.
% generic_data_flag = if enabled, the segment is excluded
%
% Author: John Paden
%
% See also: ct_filename_data, ct_filename_out, ct_filename_support,
%  ct_filename_tmp, ct_filename_gis

global gRadar;
param = merge_structs(gRadar,param);

if ~exist('generic_data_flag','var') || isempty(generic_data_flag)
  generic_data_flag = 0;
end

[output_dir,radar_type] = ct_output_dir(param.radar_name);

if isempty(fn)
  % Generate the default path
  if generic_data_flag
    fn = fullfile(param.out_path, output_dir, param.season_name, ...
      type);
  else
    fn = fullfile(param.out_path, output_dir, param.season_name, ...
      type, param.day_seg);
  end
elseif fn(1) == filesep || (ispc && (~isempty(strfind(fn,':\')) || ~isempty(strfind(fn,':/'))))
  % This is already an absolute path
  if ~generic_data_flag
    fn = fullfile(fn, param.day_seg);
  end
else
  [fn_dir fn_name] = fileparts(fn);
  fn_name = sprintf('CSARP_%s', fn_name);
  fn = fullfile(fn_dir,fn_name);
  % Generate the default path with the modified name
  if generic_data_flag
    fn = fullfile(param.out_path, output_dir, param.season_name, fn);
  else
    fn = fullfile(param.out_path, output_dir, param.season_name, fn, ...
      param.day_seg);
  end
end

fn(fn == '/' | fn == '\') = filesep;

return;
