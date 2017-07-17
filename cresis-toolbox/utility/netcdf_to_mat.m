function mat = netcdf_to_mat(cdf_fn,mat_fn,field_mask)
% mat = netcdf_to_mat(cdf_fn,mat_fn,field_mask)
%
% Reads NetCDF file created using netcdf_from_mat and writes it to
% a structure and/or Matlab file.
%
% cdf_fn: string with NetCDF4 filename created using netcdf_from_mat.m
% mat_fn: optional input argument, if supplied, the contents of the
%   NetCDF file will be written into this Matlab file
%   If empty, it is ignored.
% field_mask: optional "regular expression" for which fields to read
%   regexpi.m is used
% mat: optional output argument, cdf contents are written to this struct
%
% Examples:
%  cdf_fn = '/cresis/projects/dev/csarp_support/records/kuband2/2012_Greenland_P3/records_20120421_01.nc';
%  mat = netcdf_to_mat(cdf_fn)
%  mat = netcdf_to_mat(cdf_fn,[],'^lat$')
%  mat = netcdf_to_mat(cdf_fn,[],'^settings(1\).wfs([0-9]*\).')
%
% Author: John Paden
%
% See also: run_netcdf_from_mat.m, netcdf_type.m, netcdf_to_mat.m

if ~exist('mat_fn','var')
  mat_fn = [];
end

if ~exist('field_mask','var')
  field_mask = '.*';
end

try
  ncid = netcdf.open(cdf_fn,'NOWRITE');
catch ME
  warning('Exception during file opening');
  ME
  keyboard
end

%% Read all variables out of the netcdf file
[ndims,nvars] = netcdf.inq(ncid);
mat = [];
for var_idx = 0:nvars-1
  [varname] = netcdf.inqVar(ncid,var_idx);
  if ~isempty(regexpi(varname,field_mask))
    matlab_class = netcdf.getAtt(ncid,var_idx,'matlab_class');
    if strcmp(matlab_class,'char')
      %% Read out string (it is zero-padded in netcdf file, so we have
      % to get its actual length)
      matlab_size = netcdf.getAtt(ncid,var_idx,'matlab_size');
      eval(sprintf('mat.%s = char(netcdf.getVar(ncid,var_idx,[0 0],matlab_size));',varname));
    elseif strcmp(matlab_class,'function_handle')
      eval(sprintf('mat.%s = str2func(''%s'');',varname,netcdf.getVar(ncid,var_idx)));
    elseif strcmp(matlab_class,'cell_string')
      %% Read out cell array of strings:
      % 1. These are stored as a 2-D char array
      % 2. These are zero-padded in netcdf file, so we have to get their
      %    actual lengths)
      M = netcdf.getVar(ncid,var_idx);
      for cell_idx = 1:size(M,2)
        string_length = M(end-1)*256 + M(end);
        eval(sprintf('mat.%s{cell_idx} = M(1:string_length,cell_idx).'';',varname));
      end
    else
      matlab_size = netcdf.getAtt(ncid,var_idx,'matlab_size');
      eval(sprintf('mat.%s = %s(netcdf.getVar(ncid,var_idx));',varname,matlab_class));
      if strcmp(varname,'fasttime')
          mat.fasttime = mat.fasttime*1e-6;
      elseif strcmp(varname,'amplitude')
          mat.amplitude = 10.^(mat.amplitude/10);
      end
    end
  end
end

netcdf.close(ncid);

if ~isempty(mat_fn)
  save(mat_fn,'-struct','mat');
  if nargout == 0
    clear mat;
  end
end

return;
