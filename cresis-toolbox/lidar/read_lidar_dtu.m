function [lidar] = read_lidar_dtu(param)
% [lidar] = read_lidar_dtu(param)
%
% Reads the LIDAR data take by the Twin Otter for the hfrds2 radar
% param = struct that controls reading of file(s)
%   .time_reference = 'gps' or 'utc' (should always be 'utc' if from AWI)
%   .file_type = '.ver' Can be either .scn or .ver
%   .year = 2016 Desired campaign year (only 2016 is currently valid)
%   .month = 11 Desired campaign month (only 11/November is currently valid)
%   .day = [1,2,7,8,10,11,12] Desired campaign days of the month 
%   .dates Optional field can be filled with datetime values to find the get multiple days worth of data
% 	.season = '2016_Greenland_TOdtu'Campaign season 
%	.data_dir = 'X:/metadata' Optional field that is hardcoded as X:/metadata if it is not included

% lidar = struct of position and LIDAR data, each N x 1 vectors
%   where N is the number of records in the file(s). The fields are:
%  .gps_time = GPS time in dec.hour(UTC)
%  .lat = latitude (deg)
%  .lon = longitude (deg)
%  .surface = WGS-84 surface elevation (m)

%{
From read-me file


LIDAR *.scn files: dec.hour(UTC) latitude longitude  elevation  amplitude  #points/swath
where position refers to the reflecting point on the ground.

Be aware that there are some returns from low clouds on day 313.

The  *.ver files only contain the central scan point, format is: 
dec.hour(UTC) latitude longitude  elevation  amplitude  #points/swath  GPS.h  range
%}


%% Check the param structure
%Check for the data_dir field
basedir = 'X:/metadata';
if ~isfield(param,'data_dir')
    data_dir = basedir;
else
    if isempty(param.data_dir)
        data_dir = basedir;
    else
        data_dir = param.data_dir;
    end
end
fulldir = fullfile(data_dir,param.season,'LIDAR');

%Check for the dates field
if ~isfield(param,'dates')
    param.dates = datetime(param.year,param.month,param.day);
end

%% Read Lidar Data from the TOdtu data
%Iterate through the dates
    %Initialize the loading index
    load_id = 1;
for t = param.dates
    %Calculate the day of the year in order to grab the correct lidar files
    file_doy = num2str(day(t,'dayofyear'));

    fns = get_filenames(fulldir,file_doy,'',param.file_type);

    for fn_id = 1:length(fns)
        %Read the file
        A = dlmread(fns{fn_id});
        lidar(load_id).date = datestr(t);
        lidar(load_id).gps_time = A(:,1); %In UTC based on read_me
        lidar(load_id).lat = A(:,2);
        lidar(load_id).lon = A(:,3);
        lidar(load_id).surface = A(:,4); %This could be column 4 or 6 but the data is larger so this column is assumed to be height above MSL or GPS elevation
        load_id = load_id +1;
    end

    %Check that lidar was populated
    if isempty(fns)
        warning('There is no data for the specified date %s or in the prescribed directory %s.',t,fulldir)
    end
end
end