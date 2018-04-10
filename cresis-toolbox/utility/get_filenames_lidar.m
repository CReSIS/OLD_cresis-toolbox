function lidar_fns = get_filenames_lidar(param,lidar_source,gps_time)
% lidar_fns = get_filenames_lidar(param,lidar_source,gps_time)
%
% Gets a cell array of absolute filename strings for LIDAR files. Works
% with AWI LIDAR and DTU LIDAR data.
%
% Input:
%   param: structure controlling which LIDAR files are loaded
%     .season_name: string containing the season that the data is from
%   gps_time: start and stop times in ANSI-C time (seconds since Jan 1,
%     1970). Used to determine which files to get.
%   lidar_source: string containing 'awi' or 'dtu'
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

if ~isfield(param,'data_support_path') || isempty(param.data_support_path)
  global gRadar;
  data_support_path = gRadar.data_support_path;
else
  data_support_path = param.data_support_path;
end
  
in_base_path = fullfile(data_support_path,sprintf('LIDAR_%s',upper(lidar_source)),param.season_name);

if strcmpi(lidar_source,'awi')
  lidar_fns = get_filenames(in_base_path, 'ALS_L1B', filename_middle, '.nc');
  start_time = epoch_to_datenum(gps_time(1));
  end_time = epoch_to_datenum(gps_time(end));
  lidar_fns = {};
  cur_time = start_time;
  while cur_time <= end_time
    new_fns = get_filenames(in_base_path, 'ALS_L1B_%s', sprintf('%.0f_',day_of_year), '', '_1x1.ver');
    lidar_fns(end+(1:length(new_fns)),1) = new_fns;
    cur_time = cur_time + 1;
  end

else
  % dec.hour(UTC) latitude longitude  elevation  amplitude  #points/swath  GPS.h  range
  [year month day] = datevec(epoch_to_datenum(gps_time(1)));
  day_of_year_start = datenum(year,month,day) - datenum(year,0,0);
  [year month day] = datevec(epoch_to_datenum(gps_time(end)));
  day_of_year_end = datenum(year,month,day) - datenum(year,0,0);
  lidar_fns = {};
  for day_of_year = day_of_year_start:day_of_year_end
    new_fns = get_filenames(in_base_path, sprintf('%.0f_',day_of_year), '', '_1x1.ver');
    lidar_fns(end+(1:length(new_fns)),1) = new_fns;
  end
end

end
