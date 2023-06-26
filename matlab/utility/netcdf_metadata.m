function mat_metadata = netcdf_metadata(type)

if strcmpi(type,'L1B')
  mat_metadata.dim_names = {'fast_time','slow_time','MAT_UNIT'};
  mat_metadata.dim_lens = [NaN NaN 1];
  
  var_idx = 0;
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'Latitude';
  mat_metadata.var(var_idx).cdf_name = 'lat';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'degrees';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Latitude of measurement phase center';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'Longitude';
  mat_metadata.var(var_idx).cdf_name = 'lon';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'degrees';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Longitude of measurement phase center';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'Elevation';
  mat_metadata.var(var_idx).cdf_name = 'elev';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'meters';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Elevation of measurement phase center';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'GPS_time';
  mat_metadata.var(var_idx).cdf_name = 'gps_time';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'seconds';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'GPS time of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'Data';
  mat_metadata.var(var_idx).cdf_name = 'data';
  mat_metadata.var(var_idx).type = 'float';
  mat_metadata.var(var_idx).dims = {'fast_time','slow_time'};
  mat_metadata.var(var_idx).fill = single(NaN);
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'Watt/Watt';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Radar echogram (fast-time x slow-time), linear power scale';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'Time';
  mat_metadata.var(var_idx).cdf_name = 'time';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'fast_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'seconds';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Fast time axis';
  
elseif strcmpi(type,'records')
  mat_metadata.dim_names = {'slow_time','num_boards','MAT_UNIT'};
  mat_metadata.dim_lens = [NaN NaN 1];
  mat_metadata.var = struct([]);
  
  var_idx = 0;
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'lat';
  mat_metadata.var(var_idx).cdf_name = 'lat';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'degrees';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Latitude of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'lon';
  mat_metadata.var(var_idx).cdf_name = 'lon';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'degrees';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Longitude of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'elev';
  mat_metadata.var(var_idx).cdf_name = 'elev';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'meters';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Elevation of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'roll';
  mat_metadata.var(var_idx).cdf_name = 'roll';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'radians';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Roll of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'pitch';
  mat_metadata.var(var_idx).cdf_name = 'pitch';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'radians';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Pitch of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'heading';
  mat_metadata.var(var_idx).cdf_name = 'heading';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'radians';
  mat_metadata.var(var_idx).attr(1).name = 'datum';
  mat_metadata.var(var_idx).attr(1).val = 'WGS84';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'Heading of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'gps_time';
  mat_metadata.var(var_idx).cdf_name = 'gps_time';
  mat_metadata.var(var_idx).type = 'double';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = NaN;
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'seconds';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'GPS time of measurement';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'offset';
  mat_metadata.var(var_idx).cdf_name = 'offset';
  mat_metadata.var(var_idx).type = 'uint32';
  mat_metadata.var(var_idx).dims = {'num_boards','slow_time'};
  mat_metadata.var(var_idx).fill = uint32(2^32-1);
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'seconds';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'record offset';
  
  var_idx = var_idx + 1;
  mat_metadata.var(var_idx).name = 'settings.nyquist_zone';
  mat_metadata.var(var_idx).cdf_name = 'settings.nyquist_zone';
  mat_metadata.var(var_idx).type = 'uint8';
  mat_metadata.var(var_idx).dims = {'MAT_UNIT','MAT_UNIT','slow_time'};
  mat_metadata.var(var_idx).fill = uint8(2^8-1);
  mat_metadata.var(var_idx).attr(1).name = 'units';
  mat_metadata.var(var_idx).attr(1).val = 'seconds';
  mat_metadata.var(var_idx).attr(1).name = 'description';
  mat_metadata.var(var_idx).attr(1).val = 'record offset';
  
else
  mat_metadata.dim_names = {};
  mat_metadata.dim_lens = [];
  mat_metadata.var = struct([]);
end

end
