function nsidc_create_spatial_L2(in_fn,out_fn)
% nsidc_create_spatial_L2(in_fn,out_fn)
%
% Generate a .spatial file containing Longitude and Latitude information
% from an existing .csv data file. Same filename, different extension.
% Called from nsidc_delivery_script.m
%
% file_path = the full path of given .csv file
%
% This function aims to extract two columns(Lon and Lat) from a L2 data in 
% .csv file to a .spatial file.
%
% Author: Yi Zhu
%
% See also: type "nsidc_help.m"

%% Load the data from .csv file to a 1x9 cell array
fid_csv = fopen(in_fn);
csv_data = textscan(fid_csv,'%f%f%f%f%f%s%f%f%f','delimiter',',','headerlines',1);
fclose(fid_csv);

%% Write data into .spatial file
out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end

[fid,msg] = fopen(out_fn,'w+');

if fid < 1
  fprintf('Could not open file %s\n', out_fn);
  error(msg);
end

for data_idx = 1:length(csv_data{1})
  fprintf(fid, '%f\t%f\n',...
               csv_data{2}(data_idx),...   % Longitude
               csv_data{1}(data_idx));     % Latitude                                         
end

fclose(fid);

end
