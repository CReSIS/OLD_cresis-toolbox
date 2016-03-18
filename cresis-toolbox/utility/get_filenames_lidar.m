function lidar_fns = get_filenames_lidar(param,filename_middle)
% lidar_fns = get_filenames_lidar(param,filename_middle)
%
% Gets a cell array of absolute filename strings for LIDAR files. Works
% with AWI LIDAR data.
%
% Input:
%   param: structure controlling which LIDAR files are loaded
%     .season_name: string containing the season that the data is from
%   filename_middle: string containing expression for selecting the
%     filenames (usually a date in the YYYYMMDD format)
%
% Output:
%   lidar_fns: (cell of strings) lidar filenames (same as get_filenames output)
%
% Examples:
%   lidar_fns = get_filenames_lidar(struct('season_name','2015_Greenland_Polar6'),'20150913')
%
% Author: John Paden
%
% See also: get_filenames.m, get_filenames_atm.m

global gRadar;

if ~isfield(param,'data_support_path') || isempty(param.data_support_path)
  data_support_path = gRadar.data_support_path;
else
  data_support_path = param.data_support_path;
end

in_base_path = fullfile(data_support_path,'AWI_lidar',param.season_name);

lidar_fns = get_filenames(in_base_path, 'ALS_L1B', filename_middle, '.nc');

end
