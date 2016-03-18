function lidar = read_lidar_netcdf(awi_fns, param)
% lidar = read_lidar_netcdf(awi_fns, param)
%
% Read's AWI Netcdf LIDAR data for single beam laser and radar altimetry.
% Need to look at the specific data product to know if time reference is GPS or UTC.
%
% awi_fns = filename(s) of AWI lidar netcdf files
% param = struct that controls reading of file(s)
%   .time_reference = 'gps' or 'utc' (should always be 'utc' if from AWI)
%
% lidar = struct of position and LIDAR data, each N x 1 vectors
%   where N is the number of records in the file(s). The fields are:
%  .gps_time = GPS time in seconds since Jan 1, 1970 epoch (sec)
%  .lat = latitude (deg)
%  .lon = longitude (deg)
%  .surface = WGS-84 surface elevation (m)
%  .rms = roughness, generally root mean squared (m^2)
%
% Example: See bottom of file
%  fn = '/cresis/snfs1/dataproducts/metadata/AWI_lidar/2015_Greenland_Polar6/ALS_L1B_20150911T134604_141010_cb.nc';
%  param = struct('time_reference','utc');
%  param.nc_field = {'TIME','LATITUDE','LONGITUDE','ELEVATION','MJD','FOOT_ROUGH'};
%  param.nc_type = {'v','v','v','v','v','v'};
%  param.types = {'sec','lat_deg','lon_deg','elev_m','mjd_18581117','rms'};
%  param.scale = [1 1 1 1 1 1];
%  param.custom_flag = [0 0 0 0 0 1];
%  param.reshape_en = [1 1 1 1 1 1];
%  lidar = read_lidar_netcdf(fn,param);
%
% Author: John Paden
%
% See also get_filenames_atm.m, get_filenames_awi.m, plot_lidar.m, read_lidar_atm.m

if ~exist('param','var')
  param = struct();
end

if ischar(awi_fns)
  awi_fns = {awi_fns};
end

for file_idx = 1:length(awi_fns)
  awi_fn = awi_fns{file_idx};

  if file_idx == 1
    lidar = read_gps_netcdf(awi_fn,param);
  else
    tmp_lidar = read_gps_netcdf(awi_fn,param);
    for fieldname = fieldnames(lidar).'
      fieldname = fieldname{1};
      lidar.(fieldname) = [lidar.(fieldname) tmp_lidar.(fieldname)];
    end
  end

end

lidar = rmfield(lidar,'roll');
lidar = rmfield(lidar,'pitch');
lidar = rmfield(lidar,'heading');
lidar.surface = lidar.elev;
lidar = rmfield(lidar,'elev');

[~,new_idxs] = unique(lidar.gps_time);
for fieldname = fieldnames(lidar).'
  fieldname = fieldname{1};
  lidar.(fieldname) = lidar.(fieldname)(new_idxs);
end

return;
