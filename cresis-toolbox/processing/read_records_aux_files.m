function [records] = read_records_aux_files(records_fn,recs)
% [records] = read_records_aux_files(records_fn,recs)
% 
% Reads specific fields out of the auxilliary netcdf file for processing.
% The netcdf file is used because subsets of variables can be loaded.
%
% recs = 2x1 vector [start record, stop record]
%   Set stop record to inf to read to end of file
%   "inf" can be used for start or stop record to mean last record
%
% Examples:
%  records_fn = '/cresis/projects/dev/csarp_support/records/kuband2/2012_Greenland_P3/records_20120421_01.mat';
%  records = read_records_aux_files(records_fn,[1 4000]);
%  check = load(records_fn);
%
% Author: John Paden
%
% See also: create_records_aux_files, create_records_aux_files_fmcw_accum, 
%           create_records_mcords

if ~exist('recs')
  recs = [1 inf];
end

[path name] = fileparts(records_fn);
cdf_fn = fullfile(path, sprintf('%s.nc', name));

if recs(1) == inf
  finfo = ncinfo(cdf_fn);
  % Determine number of records and set recs(1) to this
  recs(1) = finfo.Variables(find(strcmp('lat',{finfo.Variables.Name}))).Size(2);
end

records.lat = ncread(cdf_fn,'lat',[1 recs(1)],[1 recs(2)-recs(1)+1]);
records.lon = ncread(cdf_fn,'lon',[1 recs(1)],[1 recs(2)-recs(1)+1]);
records.elev = ncread(cdf_fn,'elev',[1 recs(1)],[1 recs(2)-recs(1)+1]);
records.roll = ncread(cdf_fn,'roll',[1 recs(1)],[1 recs(2)-recs(1)+1]);
records.pitch = ncread(cdf_fn,'pitch',[1 recs(1)],[1 recs(2)-recs(1)+1]);
records.heading = ncread(cdf_fn,'heading',[1 recs(1)],[1 recs(2)-recs(1)+1]);
records.gps_time = ncread(cdf_fn,'gps_time',[1 recs(1)],[1 recs(2)-recs(1)+1]);
records.surface = ncread(cdf_fn,'surface',[1 recs(1)],[1 recs(2)-recs(1)+1]);
try
  records.adc_phase_corr_deg = ncread(cdf_fn,'adc_phase_corr_deg',[recs(1), 1],[recs(2)-recs(1)+1, 6]);
catch ME
  %warning(ME.getReport());
end
records.offset = ncread(cdf_fn,'offset',[1 recs(1)],[inf recs(2)-recs(1)+1]);
tmp = netcdf_to_mat(cdf_fn,[],'^gps_source$');
records.gps_source = tmp.gps_source;
records.settings(1).wfs_records = ncread(cdf_fn,'settings(1).wfs_records');

tmp = netcdf_to_mat(cdf_fn,[],'^relative_rec_num');
records.relative_rec_num = tmp.relative_rec_num;
tmp = netcdf_to_mat(cdf_fn,[],'^relative_filename');
records.relative_filename = tmp.relative_filename;

tmp = netcdf_to_mat(cdf_fn,[],'^param_records');
records.param_records = tmp.param_records;

tmp = netcdf_to_mat(cdf_fn,[],'^settings(1\).wfs(');
records.settings.wfs = tmp.settings.wfs;

return;



