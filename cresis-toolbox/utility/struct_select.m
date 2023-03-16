function records = struct_select(records,rec_len,idxs,depth)
% function records = struct_select(records,rec_len,start_idx,stop_idx,depth)
% Select idxs of records' fields of length rec_len
% Recursively looks up to depth levels in the  parent struct
% Written for use with radiometric calibration's crossover qloox.m
%
% Usage:
%   records = struct_select(records,rec_len,idxs,0);
%
% Author: Hara Madhav Talasila
%
% See also: 

rec_fields = fieldnames(records);

if depth >= 0
  for idx = 1:length(rec_fields)
    if eval(sprintf('length(records.%s)',rec_fields{idx})) == rec_len
      eval(sprintf('records.%s = records.%s(idxs);',rec_fields{idx},rec_fields{idx}));
    elseif eval(sprintf('isstruct(records.%s)',rec_fields{idx})) && eval(sprintf('length(records.%s)',rec_fields{idx})) == 1
      eval(sprintf('records.%s = struct_select(records.%s,rec_len,idxs,depth-1);',rec_fields{idx},rec_fields{idx}));
    end
  end
end

end