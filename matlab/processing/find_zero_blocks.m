% script find_zero_blocks
%
% DEALS WITH ERRORS IN 2012 ANTARCTICA DC8 RAW DATA
%
% This function searches through all the raw files and finds zero blocks.
% The zero blocks were caused by an error in the data archival scripts.
% This function's output is used by records_remove_zero_blocks.
%
% See also: records_remove_zero_blocks.m

% fn = '/N/dc/projects/cresis/2012_Chile_DC8/mcords/20121025/board0/mcords2_0_20121025_135855_00_0001.bin'

% params = read_param_xls('/N/u/jpaden/Quarry/scripts/branch/params-cr1/snow_param_2012_Antarctica_DC8.xls');
% params = read_param_xls('/users/paden/scripts/branch/params-cr1/snow_param_2012_Antarctica_DC8.xls');
params = read_param_xls('/users/paden/scripts/branch/params-cr1/kuband_param_2012_Antarctica_DC8.xls');
% params = read_param_xls('/users/paden/scripts/branch/params-cr1/mcords_param_2012_Antarctica_DC8.xls');
adc = 1;

for param_idx = 1:length(params)
  param = params(param_idx);
  
  if param.cmd.generic
    fprintf('=======================================================\n');
    fprintf('Finding zeros %s\n', param.day_seg);
    fprintf('=======================================================\n');
    [base_dir,adc_folder_name,fns,file_idxs] = get_segment_file_list(param,adc);
    
    zero_blocks_fn = ct_filename_support(param,'','zero_blocks');
    zero_blocks_fn_dir = fileparts(zero_blocks_fn);
    if ~exist(zero_blocks_fn_dir,'dir')
      mkdir(zero_blocks_fn_dir);
    end
    
    zero_blocks = [];
    vectors_idx = 0;
    for file_num = 1:1:length(file_idxs)
      fn = fns{file_idxs(file_num)};
      
      if 1
        [tmp fn_name] = fileparts(fn);
        fid = fopen(fn);
        A = fread(fid,inf,'uint32');
        A = A';
        fclose(fid);
        
        % Example
        % A = [1 1 0 0 0 1 1 0 1 0 0 0 1 0]
        
        % Get a list of every index that has zero
        % idxs = [3     4     5     8    10    11    12    14    15]
        idxs = find(A == 0);
        
        % Get a list where consecutive zeros occur
        % con_idxs = [1     2     5     6     8]
        diff_idxs = diff(idxs);
        con_idxs = find(diff_idxs == 1);
        
        % Get a list of breaks in consecutive zeros
        % zero_block_base = [1     2     4]
        diff_idxs = diff(con_idxs);
        zero_block_base = 1+[0 find(diff_idxs~= 1)];
        
        % zero_idxs = [3 10 14]
        zero_idxs = idxs(con_idxs(zero_block_base));
        % num_zeros = [3 3 2]
        num_zeros = 1 + [diff(zero_block_base) length(con_idxs)-zero_block_base(end)+1];
        
        good_idxs = find(num_zeros > 32);
        zero_idxs = zero_idxs(good_idxs);
        num_zeros = num_zeros(good_idxs);
        total_zeros = sum(num_zeros)*4;
        total_data = numel(A)*4;
        percent_bad = total_zeros/total_data;
        fprintf('  %40s %.1f%% (%s)\n', fn_name, ceil(100*percent_bad), datestr(now));
        zero_blocks.fns{file_num} = fn;
        zero_blocks.zero_idxs{file_num} = zero_idxs;
        zero_blocks.num_zeros{file_num} = num_zeros;
        zero_blocks.total_data(file_num) = total_data;
      end
    end
    zero_blocks.param_zero_blocks = param;
    save(zero_blocks_fn,'-struct','zero_blocks');
  end
end

