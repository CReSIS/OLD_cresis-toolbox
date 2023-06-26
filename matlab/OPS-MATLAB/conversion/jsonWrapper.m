function outData = jsonWrapper(data)
%
% outData = jsonWrapper(data)
%
% Wraps the loadjson output to match fromjson
%
% Input:
%   data: data from old json (loadjson)
%
% Output:
%   new_data: same as data, but in fromjson format.
%
% Author: Haiji Wang, Kyle W. Purdon

if ischar(data)
  outData = cell(1,size(data,1));
  for idx = 1:size(data,1)
    outData{idx} = data(idx,:);
  end
elseif isnumeric(data)
  data = data.';
  outData = cell(size(data));
  for idx = 1:numel(data)
    outData{idx} = data(idx);
  end
elseif iscell(data)
  if iscell(data{1})
    outData = cell(size(data{1},2),size(data,2));
    for result_idx = 1:size(data,2)
      for field_idx = 1:size(data{1},2)
        outData{field_idx,result_idx} = data{result_idx}{field_idx};
      end
    end
  end
end
return;
