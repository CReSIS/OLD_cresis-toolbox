function yn = check_acords_header_diff(h1,idx1,h2,idx2,check_fields,verbose)
% yn = check_acords_header_diff(h1,idx1,h2,idx2,check_fields)
%
% Detects differences between headers in a single ACORDS file.
%
%  h1: first data header for comparison
%  idx1: header index to check
%  h2: second data header for compaison (could be same as first)
%  idx2: header index to check
%  check_fields: vector of numbers indicating the index of the header
%     parameters to be compared
%
%  yn: boolean to indicate whether differences are present (0 = no
%  differences, 1 = there are differences)
%
% Author: Logan Smith
%

error(nargchk(5,6,nargin,'struct'));
if nargin == 5
  verbose = 0;
end

yn = 0;
names = fieldnames(h1);

for n_idx = check_fields
  if verbose
    fprintf('Checking %s (field %d...\n',names{n_idx},n_idx)
    fprintf('h2(%d).%s-h1(%d).%s\n',idx2,names{n_idx},idx1,names{n_idx})
    eval(sprintf('h2(%d).%s-h1(%d).%s',idx2,names{n_idx},idx1,names{n_idx}))
  end
  if eval(sprintf('h2(%d).%s-h1(%d).%s',idx2,names{n_idx},idx1,names{n_idx}))
    yn = 1;
  end
end


return