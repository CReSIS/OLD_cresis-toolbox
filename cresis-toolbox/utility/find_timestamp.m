% find_timestamp
%
% Determine which DDS.xml file was used for a particular segment of data
%
% INPUT:
%   data_dirs: List of segments to find the appropriate DDX.xml file
%   xml_dir: Location of candidate DDS.xml files
%
% OUTPUT:
%   Timestamp of the appropriate DDS.xml file
%
% Author: Logan Smith

param.type = 'd';
data_dirs = get_filenames('/cresis/data3/MCoRDS/2011_Greenland_P3/20110409B/board0/seg_08/','','','',param);
xml_fns = get_filenames('/cresis/data3/MCoRDS/2011_Greenland_P3/20110409B/','','20110409','.xml');

for xml_idx=1:length(xml_fns)
  xml_timestamps(xml_idx) = str2double(xml_fns{xml_idx}(end-9:end-4));
end

for dir_idx = 1:length(data_dirs)
  data_fns = get_filenames(data_dirs{dir_idx},'','','');
  [path name ext] = fileparts(data_fns{1});
  data_timestamp = str2double(name(20:25));
  
  right_stamps(dir_idx) = find(xml_timestamps <= data_timestamp,1,'last');
end
fprintf('\n Use the DDS file with timestamp: %d\n',xml_timestamps(right_stamps))