%  This script is used to convert CReSIS Radar Depth Sounder Data from the
%  Old format to the new format. This is the base script and has to be
%  modified to point to the appropriate folders
%
% Author: Shashanka Jagarlapudi

% base_path: Points to the folder with the datafiles in the old format

% day_segs: Segments which need to be converted. In the old format,
% multiple segments for the same day are stored with extensions of A,B,C or
% a,b,c. The script checks for these and renames them appropriately.

% gps_path: Points to the folder containing the GPS file for the
% appropriate mission. The GPS file is used to obtain Elevation which is
% missing in a few datafiles

% For some old datasets, the Surface and Bottom are in depth axis instead
% of in time axis as per the new format. If this is the case, uncomment the
% lines in the script which does this conversion.

error('Load and check onle file of the old data to make changes to this script')

base_path = '/cresis/web/cresis_data/datafiles/';

day_segs = '20091028a';
% day_segs = '20091028b';
% day_segs = '20091028c';

if ~isempty(str2num(day_segs(end)))
    day_seg = sprintf('%s_01',day_segs);
elseif double(day_segs(end)) == 97 || double(day_segs(end)) == 65
    day_seg = sprintf('%s_01',day_segs(1:end-1));
elseif double(day_segs(end)) == 98 || double(day_segs(end)) == 66
    day_seg = sprintf('%s_02',day_segs(1:end-1));
elseif double(day_segs(end)) == 99 || double(day_segs(end)) == 67
    day_seg = sprintf('%s_03',day_segs(1:end-1));
end

seg_path = fullfile(base_path,day_segs);
fns = get_filenames(seg_path,'','','.mat');
for idx = 1:length(fns);
    if ~isempty(str2num(day_segs(end)))
        file_save_name = sprintf('Data_%s_01_%03d',day_segs,idx);
    elseif double(day_segs(end)) == 97 || double(day_segs(end)) == 65
        file_save_name = sprintf('Data_%s_01_%03d',day_segs(1:end-1),idx);
    elseif double(day_segs(end)) == 98 || double(day_segs(end)) == 66
        file_save_name = sprintf('Data_%s_02_%03d',day_segs(1:end-1),idx);
    elseif double(day_segs(end)) == 99 || double(day_segs(end)) == 67
        file_save_name = sprintf('Data_%s_03_%03d',day_segs(1:end-1),idx);
    end
    
    file_save_path = fullfile(base_path,day_seg,file_save_name);
    file_save_folder = fullfile(base_path,day_seg);
    if ~exist(file_save_folder)
        mkdir(file_save_folder)
    end
    load (fns{idx});
    
    % For very old datasets like 1990's Greenland data, the naming of the
    % variables is also different. Use the following lines of code to
    % convert them. The Surface and Bottom are stored as pixels and they do
    % not have a time axis. These lines of code converts them into the new
    % format. Uncomment if running these datasets
    
    Nt = size(A,1);
    dt = 1./18.75e6;
    Time = 0:dt:dt.*(Nt-1);
    cf = thick./(bot-top);
    cf = cf.*sqrt(3.15).*2./3e8;
    top = top.*cf;
    bot = bot.*cf;
    Data = A;
    Latitude = latitude;
    Longitude = longitude;
    Surface = top;
    Thickness = thick;
    Bottom = bot;
    
    % If Elevation is missing, grab it from the GPS data
    tmp = load(spintf('%s/gps_%s.mat',gps_path,day_segs));
    Elevation = interp1(tmp.UTC_time,tmp.elev,UTC_time);
    
    % If Surface and Bottom are in Depth axis
    Surface = interp1(Depth,Time,Surface);
    Bottom = interp1(Depth,Time,Bottom);
    
    % Adding GPS_time
    leap_sec = utc_leap_seconds(UTC_Time(1));
    for t_idx = 1:length(UTC_Time)
        GPS_time(t_idx) = UTC_Time(t_idx) + leap_sec;
    end
    
    save (file_save_path,'Data','Latitude','Longitude','Thickness','Surface','Bottom','Time','UTC_time','GPS_time','Elevation')
end
