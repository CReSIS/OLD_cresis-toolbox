% script run_netcdf_from_mat
%
% Convert Matlab file outputs into NetCDF format
%
% Author: John Paden
%
% See also: netcdf_from_mat.m, netcdf_metadata.m, netcdf_to_mat.m, run_netcdf_from_mat

params = read_param_xls('/process/scripts/params/mcords_param_2012_Antarctica_DC8.xls');

cdf_dir = '/process/netcdf/';

outputs = {'CSARP_csarp-combined','CSARP_layerData'};
cdf_names = {'standard','layerData'};
outputs_post_dir = 'CSARP_post';

mat_metadata = netcdf_metadata('L1B');

for param_idx = 1:length(params)
  param = params(param_idx);
  cmd = param.cmd;
  if isfield(cmd,'generic') && cmd.generic
    for output_idx = 1:length(outputs)
      out_dir = fullfile(ct_filename_out(param,'','',1),outputs_post_dir, ...
        outputs{output_idx},param.day_seg);
      fprintf('  Output %s\n', out_dir);
      
      mat_fns = get_filenames(out_dir,'Data_','','.mat');

      for mat_idx = 1:length(mat_fns)
        mat_fn = mat_fns{mat_idx};
        [tmp mat_fn_name mat_fn_ext] = fileparts(mat_fn);
        cdf_fn_name = sprintf('%s_%s_%s_%s.nc', param.season_name, ...
          param.radar_name, cdf_names{output_idx}, mat_fn_name);
        cdf_fn = fullfile(cdf_dir,cdf_fn_name);
        fprintf('%s to\n  %s\n', mat_fn, cdf_fn);
        netcdf_from_mat(cdf_fn,mat_fn,mat_metadata)
      end
    end
  end
end

return;





