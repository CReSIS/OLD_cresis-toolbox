function merged = merge_vectors(master, slave)
% merged = merge_vectors(master, slave)
%
% Merges two vectors and removes discontinuities at the merge points.
% The master vector takes precedence. The merged output will be equal to
% the master everywhere the master is finite. Whereever the master is not
% defined, the slave vector will be used, but a "correction" will be added
% to the slave vector to remove discontinuities with the master.
%
% master: N length numeric vector
% slave: N length numeric vector
%
% merged: N length merged vector
% 
% Examples:
% 
% master = [1 2 3 NaN NaN 6 7];
% slave  = [0 1 2  3   4  5 6];
% merged = merge_vectors(master, slave)
% 
% master = [1 NaN NaN 4 NaN 6 7];
% slave  = [0 1 2  3   4  5 6];
% merged = merge_vectors(master, slave)
% 
% master = [NaN NaN 3 4 5 6 7];
% slave  = [ 0   1  2 3 4 5 6];
% merged = merge_vectors(master, slave)
% 
% master = [NaN NaN 3 4 5 NaN NaN];
% slave  = [ 0   1  2 3 4  5   6];
% merged = merge_vectors(master, slave)
% 
% master = [1 2 3 NaN NaN 6 7];
% slave  = [0 1 2  3.5   5  6.5 8];
% merged = merge_vectors(master, slave)
%
% master = [1 2 3 NaN NaN 6 7];
% slave  = [0 1 NaN  3   4  NaN 6];
% merged = merge_vectors(master, slave)
%
% master = [1 2 3 NaN NaN 6 7];
% slave  = nan(size(master));
% merged = merge_vectors(master, slave)
%
% Author: John Paden

merged = master;
slave = interp_finite(slave, NaN);

start_idx = 1;
% Find the next bad ~isfinite point
merge_idx = find(~isfinite(merged(start_idx:end)),1);
while ~isempty(merge_idx)
  merge_idx = merge_idx + start_idx - 1;
  % Find the next good isfinite point
  end_merge_idx = find(isfinite(merged(merge_idx:end)),1);
  
  %% Merge this bad section
  % merge_idx: index containing a bad point, merge_idx-1 contains a good
  %   point since merge_idx was the next first bad point
  % end_merge_idx: index containing a good point (or empty if there are no
  %   good points left in the vector
  if merge_idx == 1 && isempty(end_merge_idx)
    % master is all bad so just set merged to slave
    end_merge_idx = length(merged)+1;
    merged = slave;
  elseif merge_idx == 1
    % master has a bad section at the beginning, the correction is just
    % a constant term based on the discontinuity on the right side of the
    % bad section
    end_merge_idx = end_merge_idx + merge_idx - 1;
    merged(merge_idx:end_merge_idx-1) = slave(merge_idx:end_merge_idx-1) ...
      - (slave(end_merge_idx) - merged(end_merge_idx));
  elseif isempty(end_merge_idx)
    % master has a bad section at the end, the correction is just
    % a constant term based on the discontinuity on the left side of the
    % bad section
    end_merge_idx = length(merged)+1;
    merged(merge_idx:end_merge_idx-1) = slave(merge_idx:end_merge_idx-1) ...
      - slave(merge_idx-1) + merged(merge_idx-1);
  else
    % master has a bad section at the middle, the correction is a linear
    % interpolation from the discontinuity on the left side to the
    % discontinuity on the right side
    end_merge_idx = end_merge_idx + merge_idx - 1;
    merged(merge_idx:end_merge_idx-1) = slave(merge_idx:end_merge_idx-1) ...
      - interp1([merge_idx-1 end_merge_idx], ...
      [slave(merge_idx-1) - merged(merge_idx-1), slave(end_merge_idx) - merged(end_merge_idx)], ...
      merge_idx:end_merge_idx-1);
  end
  start_idx = end_merge_idx;
  % Find the next bad ~isfinite point
  merge_idx = find(~isfinite(merged(start_idx:end)),1);
end
