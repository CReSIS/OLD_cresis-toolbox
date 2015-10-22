function gap_idxs = data_gaps_check(master_dist, slave_dist, thresh_dist, trim_dist)
% gap_idxs = data_gaps_check(master_dist, slave_dist, thresh_dist, trim_dist)
%
% Returns indices corresponding to gaps in a dataset that is being
% interpolated. This is useful when you interpolate one dataset (e.g.
% ATM data) onto another master time reference (e.g. radar data)
% and there might be large gaps in the ATM dataset that should be
% treated differently. Process:
%   1. Interpolate data onto master axis using interp1 with extrapolation
%   2. Set all the "gap_idxs" from this function to NaN
%
% master_dist = this is the master distance axis (e.g. radar along track)
% slave_dist = this is the slave distance axis (e.g. ATM along track)
% thresh_dist = this is the threshold distance that will be considered
%   a gap in the data (e.g. 50 m)
% trim_dist = if a gap is found, setting this trim input argument will
%   cause all indices where extrapolation occurs to be considered part
%   of the gap past this value in distance (e.g. trim == 0 means no
%   extrapolation into a gap will be allowed, 20 m would allow extrapolation
%   up to 20 m from the closest slave point). trim is assumed to be less
%   than threshold.
%
% See example at bottom of file for how this works.
%
% Author: John Paden
%
% See also: data_gaps_check_mex

% Check for threshold condition to find gaps
gap_idxs = logical(zeros(size(master_dist)));
for idx = 1:length(master_dist)
  gap_idxs(idx) = min(abs(master_dist(idx) - slave_dist)) > thresh_dist;
end

% Trim the gaps
if exist('trim_dist','var')
  idx = 1;
  while idx <= length(master_dist)
    if gap_idxs(idx)
      % We have found a gap
      
      %% Trim beginning of gap
      % Find last slave point before this master point
      slave_idx = find(slave_dist < master_dist(idx), 1, 'last');
      if ~isempty(slave_idx)
        % Trim points
        gap_idxs(find(master_dist(1:idx-1) > slave_dist(slave_idx) + trim_dist)) = 1;
      end
      
      %% Trim end of gap
      % Find the end of the gap
      end_idx = find(gap_idxs(idx+1:end) == 0,1);
      if isempty(end_idx)
        % We are done since this gap goes to the end of the master axis
        break;
      else
        end_idx = idx + end_idx - 1;
      end
      
      % Find slave point at the end of this gap
      slave_idx = find(slave_dist >= master_dist(end_idx), 1, 'first');
      % Trim points
      new_gap_idxs = find(master_dist(end_idx:end) < slave_dist(slave_idx) - trim_dist);
      new_gap_idxs = end_idx - 1 + new_gap_idxs;
      gap_idxs(new_gap_idxs) = 1;

      idx = new_gap_idxs(end) + 1;
    else
      % No gap, keep going
      idx = idx + 1;
    end
  end
end

return;

% Create test data
master = [0 1 2 3 4 5 6 7 8 9 10];
master_dist = master*10 + cumsum(randn(size(master)) * 4);
master_dist = [-4.3195    6.4773   10.3932   17.4987   25.1257   36.7310   50.4995   61.7015   70.2092   83.4712   96.6667];
slave = [-0.5 4 4.5 5.5 9 9.5 9.9 10 10.5 11 11.5]
slave_dist = interp1(master,master_dist,slave);
slave = slave(~isnan(slave_dist));
slave_dist = slave_dist(~isnan(slave_dist));

% Call data_gaps_check
thresh_dist = 15;
trim_dist = 10;
gap_idxs = data_gaps_check(master_dist, slave_dist, thresh_dist, trim_dist)
gap_idxs = data_gaps_check_mex(master_dist, slave_dist, thresh_dist, trim_dist)

figure(2); clf;
plot(master,master_dist,'x')
hold on;
plot(slave,slave_dist,'ro')
plot(master(gap_idxs),master_dist(gap_idxs),'kx');
hold off; 
xlabel('Master Axis');
ylabel('Distance');
grid on;
xlim([master(1) master(end)]);
master_dist
slave_dist
