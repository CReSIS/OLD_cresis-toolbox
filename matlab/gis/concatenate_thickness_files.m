function concatenate_thickness_files(in_fns_dir,in_fns_name,out_fn,delim)
% concatenate_thickness_files(in_fns_dir,in_fns_name,out_fn,delim)
%
% Concatenates txt/csv files specified by in_fns and saves them in the
% out_fn.
% 
% Example:
%   in_fns_dir = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_post/csv/20091118_01'
%   in_fns_name = '*.csv'
%   out_fn = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_post/csv/20091118_01.csv';
%   concatenate_thickness_files(in_fns_dir,in_fns_name,out_fn,',');
%
% Author: Logan Smith, John Paden, Shashanka Jagarlapudi

fns = get_filenames(in_fns_dir,in_fns_name,'','');

fprintf('  Creating CSV %s\n', out_fn);
fid_out = fopen(out_fn,'w');
for fn_idx = 1:length(fns)
  fn = fns{fn_idx};
  fid = fopen(fn);
  if fn_idx > 1
    fgets(fid);
  end
  buffer = fread(fid,inf,'char');
  fclose(fid);
  fwrite(fid_out,buffer,'char');
end
fclose(fid_out);

return;



