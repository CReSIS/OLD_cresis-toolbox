function type = netcdf_type(type)
% type = netcdf_type(type)
%
% Converts Matlab type to NetCDF type
%
% type: if type is a string, then it must be a valid Matlab class or the
%   string 'cell_string'
%
% Outputs:
% type: numeric (double) representing netcdf equivalent class (using
%   netcdf.getConstant)
%
% Author: John Paden
%
% See also: run_netcdf_from_mat.m, netcdf_type.m, netcdf_to_mat.m

if strcmpi('int32',type)
  type = netcdf.getConstant('NC_INT');
elseif strcmpi('uint32',type)
  type = netcdf.getConstant('NC_UINT');
elseif any(strcmpi({'uint16','logical'},type))
  type = netcdf.getConstant('NC_USHORT');
elseif strcmpi('int16',type)
  type = netcdf.getConstant('NC_SHORT');
elseif any(strcmpi({'uint8','logical'},type))
  type = netcdf.getConstant('NC_UBYTE');
elseif strcmpi('int8',type)
  type = netcdf.getConstant('NC_BYTE');
elseif strcmpi('double',type)
  type = netcdf.getConstant('NC_DOUBLE');
elseif strcmpi('single',type)
  type = netcdf.getConstant('NC_FLOAT');
elseif strcmpi('single',type)
  type = netcdf.getConstant('NC_FLOAT');
elseif any(strcmpi({'char','inline','function_handle','cell_string'},type))
  type = netcdf.getConstant('NC_CHAR');
else
  type = NaN;
end

end
