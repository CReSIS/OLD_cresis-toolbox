function create_records_aux_files(records_fn,print_flag)
% create_records_aux_files(records_fn,print_flag)
%
% Creates the auxilliary files to the main records file.  These files
% are used by tasks running on the scheduler because they allow a much
% smaller memory footprint to load (currently Matlab .mat files require
% that entire variables be loaded rather than subsets).
%
% The auxilliary file is a NetCDF file.
%
% Examples: See bottom of file for example to create files in a batch.
%
% Author: John Paden
%
% See also: read_records_aux_files, create_records_aux_files_fmcw_accum,
%           create_records_mcords

if ~exist('print_flag','var') || isempty(print_flag)
  print_flag = true;
end

% =====================================================================
% Create records netcdf file
%   -- This file allows quick loading of specific records rather than
%      having to load the entire file.
% =====================================================================

[path name] = fileparts(records_fn);
cdf_fn = fullfile(path, sprintf('%s.nc', name));
if print_flag
  fprintf('  Creating file %s\n  from %s\n', cdf_fn, records_fn);
end
netcdf_from_mat(cdf_fn,records_fn);

return;


% =====================================================================
% =====================================================================
% Examples
% =====================================================================
% =====================================================================

records_path = '/cresis/scratch1/mdce/csarp_support/records/mcords/2009_Antarctica_DC8/';
records_fns = get_filenames(records_path,'records','','.mat');

for file_idx = 1:length(records_fns)
  create_records_aux_files(records_fns{file_idx});
end
