% script create_records_save_workspace
%
% Saves the workspace to a temporary location

%% Save workspace in case there is a failure
if isfield(param,'tmp_path') && ~isempty(param.tmp_path)
  fn = ct_filename_ct_tmp(param,'','records','workspace');
  fn = [fn '.mat'];
  fprintf('Saving workspace %s (%s)\n', fn, datestr(now));
  try
    [fn_path,name] = fileparts(fn);
    if ~exist(fn_path,'dir')
      mkdir(fn_path);
    end
    save(fn);
  catch
    fprintf('  Saving workspace failed\n');
  end
end