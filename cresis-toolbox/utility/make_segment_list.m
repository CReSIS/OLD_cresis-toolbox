function segs = make_segment_list(xls_fn,param)
% segs = make_segment_list(xls_fn,param)
%
% xls_fn = parameter spreadsheet (.xls) filename or directory (for multiple
%          .csv files)
% param = structure controlling operation of the function
%  .all = boolean; when true it will process all segments, ignoring the
%         param.vectors.verification field
% segs = cell vector of segments
%
% Examples:
%  segs = make_segment_list('/users/paden/scripts/matlab/kuband_param_2011_Greenland_P3_dravid2.xls');
%  segs = make_segment_list('cresis/scratch2/mdce/mcords/2011_Antarctica_DC8/CSARP_post/csv');
%
% Author: John Paden
% Contributions by: Steve Foga

if ~exist('param','var') || ~isfield('param','all')
  param.all = false;
end

if ~isdir(xls_fn) % If an individual param file is called into the function
params = read_param_xls(xls_fn);

segs = {};


  for seg_idx = 1:length(params)
    if param.all || isempty(strfind(lower(params(seg_idx).vectors.verification),'do not process'))
      segs{end+1} = params(seg_idx).day_seg;
    end
  end
  
else % Else a directory (where CSV files exist) is called into the function
  files = get_filenames(xls_fn,'Data_','','.csv');  
  segs = cell(length(files),1)';
  for file_idx = 1:length(files)
    segs{file_idx} = files{file_idx}(end-14:end-4);
  end
  
  
end
return;
