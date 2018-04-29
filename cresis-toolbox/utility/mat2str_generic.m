function str = mat2str_generic(M)
% str = mat2str_generic(M)
%
% Expanded capability over Matlab's mat2str. Supports cell arrays.
% Also supports simplifcation of regular structures in matrices such as
% [1 1 1 1 1] getting turned into ones(1,5).
%
% Author: John Paden

if numel(size(M)) > 2
  error('Only 1D and 2D matrices are supported');
end

if iscell(M)
  % '{[1 2; 1 3; 1 4; 1 5],[]; [], [3 4 3]}'
  str = '{';
  for row = 1:size(M,1)
    if row > 1
      str = cat(2,str,';');
    end
    for col = 1:size(M,2)
      if col > 1
        str = cat(2,str,',');
      end
      str = cat(2,str,mat2str_generic(M{row,col}));
    end
  end
  str = cat(2,str,'}');
  
elseif isnumeric(M)
  diff_M = diff(M);
  if numel(M) <= 4
    str = mat2str(M);
  elseif ~isempty(M) && all(M(:) == M(1))
    str = sprintf('%g*ones(%s)', M(1), mat2str(size(M)));
  elseif (size(M,1) == 1 || size(M,2) == 1) && numel(M)>2 && all(diff_M == diff_M(1))
    if diff_M(1) == 1
      str = sprintf('[%g:%g]',M(1),M(end));
    else
      str = sprintf('[%g:%g:%g]',M(1),diff_M(1),M(end));
    end
    if size(M,2) == 1
      str = cat(2,str,'.''');
    end
  else
    regular = false;
    if all(size(M) > [1 1])
      % Check each column for regularity
      regular = true;
      for col = 1:size(M,2)
        Mcol = M(:,col);
        diff_M = diff(Mcol);
        if ~all(Mcol(:) == Mcol(1)) && ~all(diff_M == diff_M(1))
          regular = false;
          break;
        end
      end
    end
    if regular
      str = '[';
      for col = 1:size(M,2)
        Mcol = M(:,col);
        diff_M = diff(Mcol);
        if col > 1
          str = cat(2,str,',');
        end
        if all(Mcol(:) == Mcol(1))
          if Mcol(1) == 1
            str = cat(2,str,sprintf('ones(%s)', mat2str(size(Mcol))));
          elseif Mcol(1) == 0
            str = cat(2,str,sprintf('zeros(%s)', mat2str(size(Mcol))));
          else
            str = cat(2,str,sprintf('%g*ones(%s)', Mcol(1), mat2str(size(Mcol))));
          end
        elseif all(diff_M == diff_M(1))
          if diff_M(1) == 1
            str = cat(2,str,sprintf('(%g:%g).''',Mcol(1),Mcol(end)));
          else
            str = cat(2,str,sprintf('(%g:%g:%g).''',Mcol(1),diff_M(1),Mcol(end)));
          end
        end
      end
      str = cat(2,str,']');
    else
      str = mat2str(M);
    end
  end
end

end