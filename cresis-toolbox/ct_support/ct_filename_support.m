function fn = ct_filename_support(param,fn,type,generic_data_flag)
% fn = ct_filename_support(param,fn,type,generic_data_flag)
%
% Returns a standardized filename for support files.
% 1. Handles absolute and relative path conversions, default paths, and
%    file separator differences (e.g. windows uses '\' and unix uses '/')
% 2. There are three modes of operation:
%  - base_fn is an absolute path: the param.support_path and standardized
%    path are not used (base_fn replaces both)
%  - base_fn is a relative path: the param.support_path is used, but the
%    standardized path/directory structure is not used (base_fn is used)
%  - base_fn is empty: the param.support_path and standardized path are used
%
% param = control structure to data processor
%  .radar_name (e.g. mcords)
%  .season_name (e.g. 2009_antarctica_DC8)
%  .day_seg (e.g. 20091020_03)
% fn = parameter filename provided (leave blank to get the default file
%   path)
% type = string containing the type of data (e.g. vectors, gps, records,
%   frames, radar_config)
% generic_data_flag = if enabled, the radar name and segment are excluded.
%   This field is optional. Default is false.
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

if strcmpi(type,'records') && isfield(param,'records') && isfield(param.records,'records_fn') && isempty(fn)
  fn = param.records.records_fn;
end

if strcmpi(type,'frames') && isfield(param,'records') && isfield(param.records,'frames_fn') && isempty(fn)
  fn = param.records.frames_fn;
end

if isempty(fn)
  if ~isfield(param,'day_seg') || isempty(param.day_seg)
    % Generate the default path
    if generic_data_flag
      fn = fullfile(param.support_path ,type, param.season_name);
    else
      fn = fullfile(param.support_path, type, output_dir, ...
        param.season_name);
    end
  else
    % Generate the default filename
    if generic_data_flag
      fn = fullfile(param.support_path ,type, param.season_name, ...
        sprintf('%s_%s.mat', type, param.day_seg(1:8)));
    else
      fn = fullfile(param.support_path, type, output_dir, ...
        param.season_name, sprintf('%s_%s.mat', type, param.day_seg));
    end
  end
elseif fn(1) == filesep || (ispc && (~isempty(strfind(fn,':\')) || ~isempty(strfind(fn,':/'))))
  % This is already an absolute path
  return
else
  % Append the current path to the support path
  fn = fullfile(param.support_path, fn);
end

fn(fn == '/' | fn == '\') = filesep;

return;
