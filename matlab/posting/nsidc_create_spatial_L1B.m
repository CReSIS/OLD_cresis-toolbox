function nsidc_create_spatial_L1B(in_fn,out_fn)
% nsidc_create_spatial_L1B(in_fn,out_fn)
%
% Generate a .spatial file containing Longitude and Latitude information
% from an existing .mat data file. Same filename, different extension.
% Called from nsidc_delivery_script.m
%
% file_path = the full path of given .mat file
%
% This function aims to extract two columns(Lon and Lat) from a L1B data in 
% .mat file to a .spatial file.
%
% Author: Yi Zhu
%
% See also: type "nsidc_help.m"

%% Generate the .spatial file. Obtain filename from original .csv data file
metaData = load(in_fn);

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


for data_idx = 1:length(metaData.Longitude)
  fprintf(fid, '%f\t%f\n',...
               metaData.Longitude(data_idx),...   % Longitude
               metaData.Latitude(data_idx));      % Latitude                                         
end

fclose(fid);

end
