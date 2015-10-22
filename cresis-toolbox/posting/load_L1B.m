function mdata = load_L1B(fn)
% function load_L1B(fn)
%
% fn = filename string and must contain correct extension ('.nc' or '.mat')
%
% Loads L1B cresis echogram files.
% L1B files can be in netcdf or mat format.
% They can be "compressed" (compression is used for Snow and Kuband radar data
% to make them smaller... main difference is that the surface is tracked
% so that the range gate stored in the file changes on each range line).
%
% Example:
%  fn = 'IRMCR1B_V01_20130408_01_020.nc';
%  mdata = load_L1B(fn);
%
%  fn = 'Data_20111107_02_191.mat';
%  mdata = load_L1B(fn);
%
% Author: John Paden
%
% See also: plot_L1B.m, uncompress_echogram.m

[fn_dir,fn_name,fn_ext] = fileparts(fn);

if strcmpi(fn_ext,'.nc')
  mdata = netcdf_to_mat(fn);
  
  mdata.Latitude = mdata.lat;
  mdata = rmfield(mdata,'lat');
  mdata.Longitude = mdata.lon;
  mdata = rmfield(mdata,'lon');
  mdata.Elevation = mdata.alt;
  mdata = rmfield(mdata,'alt');
  mdata.Data = mdata.amplitude;
  mdata = rmfield(mdata,'amplitude');
  mdata.GPS_time = mdata.time;
  mdata = rmfield(mdata,'time');
  mdata.Time = mdata.fasttime;
  mdata = rmfield(mdata,'fasttime');
  
elseif strcmpi(fn_ext,'.mat')
  mdata = load(fn);
  
else
  error('Unsupported extention %s', fn_ext);
end

mdata = uncompress_echogram(mdata);

return;
