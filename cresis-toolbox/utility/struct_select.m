function rrr = struct_select(rrr,matching_len,idxs,depth)
% function rrr = struct_select(rrr,rec_len,start_idx,stop_idx,depth)
% Select idxs of rrr' fields of length rec_len
% Recursively looks up to depth levels in the  parent struct
% Written for use with radiometric calibration's crossover qloox.m
%
% Usage:
%   rrr = struct_select(rrr,rec_len,idxs,0);
%
% Author: Hara Madhav Talasila
%
% See also:

rec_fields = fieldnames(rrr);

if depth >= 0
  for idx = 1:length(rec_fields)
    if eval(sprintf('isstruct(rrr.%s)',rec_fields{idx})) && eval(sprintf('length(rrr.%s)',rec_fields{idx})) == 1
      eval(sprintf('rrr.%s = struct_select(rrr.%s,matching_len,idxs,depth-1);',rec_fields{idx},rec_fields{idx}));
    elseif eval(sprintf('isvector(rrr.%s)',rec_fields{idx})) && ...
        eval(sprintf('length(rrr.%s)',rec_fields{idx})) == matching_len
      eval(sprintf('rrr.%s = rrr.%s(idxs);',rec_fields{idx},rec_fields{idx}));
    elseif eval(sprintf('ismatrix(rrr.%s)',rec_fields{idx})) % && ...
      matching_dim = eval(sprintf('find(size(rrr.%s) == matching_len)',rec_fields{idx})); %== 2
      if ~isempty(matching_dim) && matching_dim == 2
        eval(sprintf('rrr.%s = rrr.%s(:,idxs);',rec_fields{idx},rec_fields{idx}));
      end
    else
    end
  end
end

end