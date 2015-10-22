function out = initialization(DCM,wb_param)
% out = initialization(DCM,wb_param)
%
% Function used to initialize wideband estimator called in array_proc.m.
% initialization uses the alternating projection approach to initialize.
% The cost function is evaluating over a coarse grid to get a rough
% estimate of the minimimum and then a 2nd order polynomial fit is used to
% refine that estimate.
%
% Inputs:
%   DCM = complex space-time data covariance matrix with dimensions (P*W) x
%     (P*W) where P is the number of sensors and W is the widening factor
%     passed into the combine spreadsheet and accessed through
%     array_param.W.
%
%   wb_param = control struction containing the following fields:
%     .Nsig   = number of sources specified by array_param.Nsig field,
%     .src_limits = 1 x wb_param.Nsig cell containing bounds for the 
%               search, specified as DOAs in radians, of the form:
%               wb_param.src_limits{source index} = [lowerbound upperbound]
%                 or the default:
%               wb_param.src_limits{source index} = ''
%     .y_pc   = column vector of y coordinates of each phase center in SAR
%               flight coordinate system,
%     .z_pc   = column vector of z coordinates of each phase center in SAR
%               flight coordinate system,
%     .fc     = center frequency,
%     .fs     = sampling frequency of decimated SAR outputs,
%     .h      = impulse response time series computed in
%               combine_wf_chan_task.m
%     .t0     = start time of the impulse response time series,
%     .dt     = sampling interval of the impulse response time series,
%     .theta  = coarse grid used to evaluate cost function.  The number of
%               elements in S.theta is set by the Nsv field in the combine
%               spreadsheet.  S.theta is the FFTSHIFTED theta vector output
%               by the array_proc_sv call in array_proc.
%     .src_limits
%
% Outputs:
%   out = Nsig x 1 vector containing initial doa estimates in radians.
%
% Author:  Theresa Stumpf
%
% See Also:  array_proc.m, initialization.m, compute_cost.cpp
% =========================================================================

if ~isfield(wb_param,'src_limits');
  for src_idx = 1:wb_param.Nsig
    wb_param.src_limits{src_idx} = '';
  end
end


physical_constants;

t0    = wb_param.t0;
dt    = wb_param.dt;
h     = wb_param.h(:);


% =========================================================================
%% Initialization of first source
% =========================================================================

% Setup search space
% -----------------------------------------------------------------------
%
% If wb_param.src_limits was specified, only look at thetas bounded by
% wb_param.src_limits.  Otherwise use the entire coarse grid stored in wb_param.theta.
if ~isempty(wb_param.src_limits{1})
  indexes = (wb_param.theta >= wb_param.src_limits{1}(1)) & (wb_param.theta <= wb_param.src_limits{1}(end));
  search_theta = wb_param.theta(indexes);
else
  search_theta = wb_param.theta;
end

% Evalue cost cunction over coarse grid whose step size is set by
% array_param.Nsv
J   = NaN(size(search_theta));
for idx = 1:length(search_theta);
  theta   = search_theta(idx);
  uy      = sin(theta).';
  uz      = cos(theta).';
  tau     = (2/c)*(-wb_param.y_pc*uy - wb_param.z_pc*uz);
  val     = compute_cost(tau, DCM, h, t0, dt);
  J(idx)  = -10*log10(abs(val));
end

% Minimize Cost Function
% -----------------------------------------------------------------------
% Find minimium of J.  If min falls on an edge, then just use the edge.
% Otherwise, find minimum index and adjacent neighbors.  Then use a
% second order polynomial fit to J to refine estimate.
[~,min_idx] = nanmin(J);

if min_idx == 1 || min_idx == length(search_theta)
  theta0(1,1) = search_theta(min_idx);
else
  theta_rng = search_theta(min_idx - 1: min_idx + 1);
  Jvals     = J(min_idx - 1: min_idx + 1);
  p         = polyfit(theta_rng(:),Jvals(:),2);
  a         = p(1);
  b         = p(2);
  theta0(1,1) = -b/(2*a);
end


% =========================================================================
%% Initialization of additional sources
% =========================================================================
if wb_param.Nsig > 1
  
  % Loop over sources and estimate each direction
  for src_idx = 2:wb_param.Nsig;
    clear J uy uz min_idx search_theta
    
    % If wb_param.src_limits was specified, only look at thetas bounded by
    % wb_param.src_limits.  Otherwise use the entire coarse grid stored in wb_param.theta.
    if ~isempty(wb_param.src_limits{src_idx})
      indexes = (wb_param.theta >= wb_param.src_limits{src_idx}(1)) & (wb_param.theta <= wb_param.src_limits{src_idx}(end));
      search_theta = wb_param.theta(indexes);
    else
      search_theta = wb_param.theta;
    end
    
    % Avoid evaluating cost function around previous sources
    prev_doa = theta0(theta0 >= search_theta(1) & theta0 <= search_theta(end));
    if ~isempty(prev_doa)
      for prev_idx = 1:length(prev_doa)
        mask_doa = prev_doa(prev_idx);
        
        bad_index = find(search_theta > prev_doa(prev_idx),1) - 1;
        
        if ~isempty(bad_index)
          
          bad_rng = bad_index - 1:bad_index + 1;
          bad_mask = bad_rng < 1 | bad_rng > length(search_theta);
          bad_rng = bad_rng(~bad_mask);
          search_theta(bad_rng) = NaN;    
        end
        
%         if mask_doa == search_theta(1) 
%           bad_rng = 1;
%         elseif mask_doa == search_theta(end)
%           bad_rng = length(search_theta);      
%         else
%           bad_index = find(search_theta > prev_doa(prev_idx),1) - 1;
%           if ~isempty(bad_index)
%             bad_rng = bad_index - 1: bad_index+1;
%           else
%             bad_rng = '';
%           end
%         end
%         
%         if ~isempty(bad_rng)
%           bad_mask = bad_rng<1 | bad_rng > length(search_theta);
%           bad_rng = bad_rng(bad_mask);
%           search_theta(bad_rng) = NaN;
%         end
      end
    end
    
    search_theta  = search_theta(~isnan(search_theta));
    
    % Evalue cost cunction over coarse grid whose step size is set by
    % array_param.Nsv
    J = NaN(size(search_theta));
    for idx = 1:length(search_theta);
      theta_i = search_theta(idx);
      theta   = [theta_i;theta0];
      uy      = sin(theta).';
      uz      = cos(theta).';
      tau     = (2/c)*(-wb_param.y_pc*uy - wb_param.z_pc*uz);
      val     = compute_cost(tau,DCM,h,t0,dt);
      J(idx)  = -10*log10(abs(val));
    end
    
    [~,min_idx] = nanmin(J);
    
    if min_idx == 1 || min_idx == length(search_theta)
      theta0(src_idx,1) = search_theta(min_idx);
    elseif any(isnan(J(min_idx - 1:min_idx + 1)))
      theta0(src_idx,1) = search_theta(min_idx);
    else
      theta_rng = search_theta(min_idx - 1: min_idx + 1);
      Jvals     = J(min_idx - 1: min_idx + 1);
      p         = polyfit(theta_rng(:),Jvals(:),2);
      a         = p(1);
      b         = p(2);
      theta0(src_idx,1) = -b/(2*a);
    end
    theta0 = sort(theta0);
    
      
%     
%     mask(min_idx) = false;
%     loop_vec      = 1:length(wb_param.theta);
%     
%     if 1 && wb_param.Nsig == 2
%       
%       if theta0(1,1) > 0
%         theta_mask = logical(wb_param.theta < 0);
%         loop_vec    = loop_vec(theta_mask);
%       elseif theta0(1,1) < 0
%         theta_mask = logical(wb_param.theta > 0);
%         loop_vec  = loop_vec(theta_mask);
%       else
%         loop_vec = loop_vec(mask);
%       end
%     else
%       loop_vec = loop_vec(mask);
%     end
%     
%     
%     clear J theta uy uz min_idx
%     J       = NaN(size(wb_param.theta));
%     for idx = loop_vec;
%       theta_i = wb_param.theta(idx);
%       theta   = [theta_i;theta0];
%       uy      = sin(theta).';
%       uz      = cos(theta).';
%       tau     = (2/c)*(-wb_param.y_pc*uy - wb_param.z_pc*uz);
%       val     = compute_cost(tau,DCM,h,t0,dt);
%       J(idx)  = -10*log10(abs(val));
%     end
%     
%     [~,min_idx] = nanmin(J);
    
%     if min_idx == 1 || min_idx == length(wb_param.theta)
%       theta0(src_idx,1) = wb_param.theta(min_idx);
%     elseif any(isnan(J(min_idx - 1:min_idx + 1)))
%       theta0(src_idx,1) = wb_param.theta(min_idx);
%     else
%       theta_rng = wb_param.theta(min_idx - 1: min_idx + 1);
%       Jvals     = J(min_idx - 1: min_idx + 1);
%       p         = polyfit(theta_rng(:),Jvals(:),2);
%       a         = p(1);
%       b         = p(2);
%       theta0(src_idx,1) = -b/(2*a);
%     end
%     theta0 = sort(theta0);
  end
  
end
out.doa = theta0;

end






