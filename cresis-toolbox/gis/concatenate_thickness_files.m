function concatenate_thickness_files(in_search_str,out_fn,delim)
% concatenate_thickness_files(in_search_str,out_fn,delim)
%
% Concatenates txt/csv files specified by in_search_str and saves them in the
% out_fn.
% 
% Example:
%   in_search_str = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_post/csv/20091118_01/*.csv'
%   out_fn = '/cresis/scratch2/mdce/mcords/2009_Antarctica_DC8/CSARP_post/csv/20091118_01.csv';
%   concatenate_thickness_files(in_search_str,out_fn,',');
%
% Author: Logan Smith, John Paden, Shashanka Jagarlapudi

sys_cmd = sprintf('cat %s  | grep -v LAT | sed 1i"LAT%sLON%sTIME%sTHICK%sELEVATION%sFRAME%sSURFACE%sBOTTOM%sQUALITY" >%s', ...
  in_search_str,delim,delim,delim,delim,delim,delim,delim,delim,out_fn);

system(sys_cmd);

return;



