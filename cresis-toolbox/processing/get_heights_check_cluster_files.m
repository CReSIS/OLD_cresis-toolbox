% script get_heights_check_cluster_files
%
% Temporary debug function to track down compute system failure. Bad files
% are getting written to the file system during qlook generation.
%   Dec 2015
%
% Function called by get_heights.m before it cleans up torque jobs.
%
% Author: John Paden

if 1
  if isempty(param.get_heights.qlook.out_path)
    test_qlook_out_path = 'qlook';
  else
    test_qlook_out_path = param.get_heights.qlook.out_path;
  end
  test_dir = ct_filename_out(param,test_qlook_out_path,'')
  test_yymmdd = param.day_seg(1:8)
  
else
  
  % test_yymmdd = param.day_seg(1:8)
  % test_dir =
  % /cresis/snfs1/dataproducts/ct_data/snow/2009_Antarctica_DC8/CSARP_qlook/20091102_05
  % test_yymmdd =
  % 20091102
  
  % test_dir =
  % /cresis/snfs1/dataproducts/ct_data/snow/2011_Antarctica_DC8/CSARP_qlook/20111112_01
  % test_yymmdd =
  % 20111112
  
  test_dir ='/cresis/snfs1/dataproducts/ct_data/snow/2010_Antarctica_DC8/CSARP_deconv/20101030_01'
  test_yymmdd ='20101030'
end

dbstack_info = dbstack;
fprintf('%s: %s (%s)\n', dbstack_info(1).name, test_yymmdd, datestr(now,'HH:MM:SS'));

test_fns = get_filenames(test_dir,test_yymmdd,'','.mat','recursive');
% test_fns = get_filenames(test_dir,test_yymmdd,'','.mat','recursive');

bad_files = false;
move_files = true;
for fn_idx = 1:length(test_fns)
  test_fn=test_fns{fn_idx};
  try
    test = load(test_fn,'GPS_time');
  catch
    bad_files = true;
    fprintf('%s\n', test_fn);
    test_fn_info = dir(test_fn);
    fprintf('  %s\n', test_fn_info.date);
    
    test_fn_dir = fileparts(test_fn);
    [~,test_fn_dir_name] = fileparts(test_fn_dir);
    
    frm = str2double(test_fn_dir_name(9:11));
    block_idx = find(strcmp(get_filenames(test_fn_dir,'','','.mat'),test_fn));
    fprintf('  %d %d\n', frm, block_idx);
    
    if move_files
      movefile(test_fn,fullfile(gRadar.tmp_path,'bad_get_heights_files'));
    end
  end
end

if bad_files
  warning('Bad files were found and printed to the screen. Email helpdesk to let them know. To fix the bad files, determine which jobs need to be rerun with torque_exec_job(ctrl,[LIST OF JOB IDS]). Once they have been recreated, then type dbcont.')
  keyboard
end

return;
