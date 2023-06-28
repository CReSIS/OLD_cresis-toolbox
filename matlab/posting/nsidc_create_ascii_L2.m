function nsidc_create_ascii_L2(in_fn,data_fn)
% nsidc_create_ascii_L2(in_fn,data_fn)
%
% Convert a .csv file to .txt file in ASCII File Format Convention.
% Called from nsidc_delivery_script.m
%
% file_path = the full path of given .csv file
%
% This function aims to convert a L2 data in .csv file (with 9 columns) to
% a standard .txt file.
%
% Author: Yi Zhu
%
% See also: type "nsidc_help.m"


%% Load the data from .csv file to a 1x9 cell array
fid_csv = fopen(in_fn);
csv_data = textscan(fid_csv,'%f%f%f%f%f%f%f%f%f','delimiter',',','headerlines',1);
fclose(fid_csv);

%% Manully input Header documentation about the data
data_fn_dir = fileparts(data_fn);
if ~exist(data_fn_dir,'dir')
  mkdir(data_fn_dir);
end
[fid,msg] = fopen(data_fn,'w+');
if fid < 1
  fprintf('Could not open file %s\n', data_fn);
  error(msg);
end

fprintf(fid, '# The MCoRDS L2 Ice Thickness data set contains measurements for Elevation, Surface, Bottom and Thickness.\n');
fprintf(fid, '# Format is 9 columns separated by commas. Columns are:\n');
fprintf(fid, '# \n');
fprintf(fid, '# LAT           Latitude 	Degrees North, Antenna position on the aircraft w.r.t WGS-84 and ITRF2008 from postprocessed GPS data.\n');
fprintf(fid, '# LON           Longitude 	Degrees East, Antenna position on the aircraft w.r.t WGS-84 and ITRF2008 from postprocessed GPS data.\n');
fprintf(fid, '# TIME          UTC Time 	Seconds of day. Note: When aligning with GPS time tagged data, account for leap seconds.\n');
fprintf(fid, '#               Note: time tages can be obtained by using the YYYY MM DD fields from the FRAME and adding TIME.\n');
fprintf(fid, '# THICK         Ice Thickness: Bottom minus Surface. Constant dielectric of 3.15 (no firn) is assumed for converting propagation delay into range. -9999 indicates no thickness available. 	Meters\n');
fprintf(fid, '# ELEVATION 	Elevation of GPS antenna referenced to WGS-84 Ellipsoid and ITRF2008 . 	Meters\n');
fprintf(fid, '# FRAME 	(YYYYMMDDSSFFF) Fixed length numeric field. YYYY = year, MM = month, DD = day, SS = segment FFF = frame.\n');
fprintf(fid, '# SURFACE       Range to Ice Surface. Actual surface height is Elevation minus this number. 	Meters\n');
fprintf(fid, '# BOTTOM        Range to Ice Bottom. Actual ice bottom height is Elevation minus this number. Constant dielectric of 3.15 (no firn) is assumed for converting propagation delay into range. -9999 indicates no thickness available. 	Meters\n');
fprintf(fid, '# QUALITY 	1: High confidence pick\n');
fprintf(fid, '#               2: Medium confidence pick\n');
fprintf(fid, '#               3: Low confidence pick\n');
fprintf(fid, '# \n');
fprintf(fid, '# ');
fprintf(fid,'%s,%s,%s,%s,%s,%s,%s,%s,%s\n', 'LAT','LON','TIME','THICK',...
        'ELEVATION','FRAME','SURFACE','BOTTOM','QUALITY');
      

%% Print out the data line by line into the .txt file
bad_mask = ~isfinite(csv_data{9});
if sum(bad_mask) > 0
  warning('NaN quality found in %d points', sum(bad_mask));
  if sum(bad_mask) > 1
    keyboard;
  end
end
for data_idx = find(~reshape(bad_mask,[1 length(bad_mask)]))
  fprintf(fid, '%8.6f,%8.6f,%8.4f,%8.2f,%8.4f,%8.0f,%8.2f,%8.2f,%.0f\n',...
                    csv_data{1}(data_idx),...
                    csv_data{2}(data_idx),...
                    csv_data{3}(data_idx),...
                    csv_data{4}(data_idx),...
                    csv_data{5}(data_idx),...
                    csv_data{6}(data_idx),...
                    csv_data{7}(data_idx),...
                    csv_data{8}(data_idx),...
                    csv_data{9}(data_idx));                                            
end

fclose(fid);

end
