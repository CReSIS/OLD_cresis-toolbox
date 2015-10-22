% Script: remove_variable.m
%
% Removes a specific set of variables from a list of files

%% User Settings
fns = get_filenames('/cresis/snfs1/dataproducts/public/data/accum/2013_Greenland_P3/CSARP_qlook','Data_','','.mat',struct('recursive',1));

variables_to_remove = {'Depth'}

%% Automated Section
for fn_idx = 1:length(fns)
  fprintf('%s\n', fns{fn_idx});
  for var_idx = 1:length(variables_to_remove)
    status = remove_variable_support(fns{fn_idx},variables_to_remove{var_idx});
    if status ~= 0
      fprintf('  remove_variable_support failed to delete variable (%d)\n', status);
    end
  end
end

