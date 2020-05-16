function list_out = merge_filelists(list_in,list_override)
% list_out = merge_filelists(list_in,list_override)
%
% Merges two lists of files (cell vectors).  Any filename present
% in both lists will only be passed in from the list_override list
% The idea is that the paths to the file will be different.
%
% Used to create file dependencies list for cluster/parallel
% processing programs.
%
% list_in: cell vector of file paths
% list_override: cell vector of file paths
% list_out: cell vector of file paths
% 
% Author: John Paden
%
% See also: 

list_out = list_override;

% Create list of filenames only
list_out_fns = {};
for out_idx = 1:length(list_out)
  [path list_out_fns{out_idx}] = fileparts(list_out{out_idx});
end

for in_idx = 1:length(list_in)
  [path name] = fileparts(list_in{in_idx});
  if isempty(strmatch(name,list_out_fns,'exact'))
    % File is not present in output list, so we need to add it:
    list_out{end+1} = list_in{in_idx};
  end
end

return;
