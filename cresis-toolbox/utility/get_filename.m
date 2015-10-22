function [filenames] = get_filename(filepath,filename_start,filename_middle,filename_end,varargin)
% [filenames] = get_filename(filepath,filename_start,filename_middle,filename_end,varargin)
%
% Calls get_filenames, but only returns one matchign file and 
% array.
%
% Author: John Paden

[filenames] = get_filenames(filepath,filename_start,filename_middle,filename_end,varargin{:});

if length(filenames) > 1
  for idx = 1:length(filenames)
    fprintf('  (%d) %s\n', idx, filenames{idx});
  end
  idx = input('Choose a file: ');
  filenames = filenames{idx};
elseif length(filenames) == 1
  filenames = filenames{1};
else
  error('No files found\n');
end

return

