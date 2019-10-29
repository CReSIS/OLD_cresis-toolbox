function [J, grad] = moe_logistic_regression_cost_fn(theta_q,param)
x_q = param.training_input;
y_q = param.training_output;
if isfield(param,'model') && ~isempty(param.model)
  h = param.model;
else
  h = @(theta,x) (1./(1+exp(-(theta.'*x))));
end

Nexamples = size(x_q,2);

I = ones(length(y_q),1);

if isfield(param,'reg_param') && ~isempty(param.reg_param)
  lambda = param.reg_param;
else
  lambda = 0e2;
end

reg_term = lambda * (1/(2*Nexamples) * sum(theta_q(2:end).^2)); 

% inf * 0 = NaN. So, to avoid having inf we can do something like setting
% any inf values to a big number, such as 10e10; or set any value that is
% =0 inside the log to a small number.

% L0 = I.' - h(theta_q,x_q);
% L0 = real(log(L0));
% 
% L1 = h(theta_q,x_q);
% L1 = real(log(L1));

% L0_inf_percentage = length(find(isinf(L0)))/length(L0);
% L1_inf_percentage = length(find(isinf(L1)))/length(L1);

% if lambda >= 0 && lambda <= 10
%   a = 50/100;
% elseif lambda > 10 && lambda <= 100
%   a = 40/100;
% elseif lambda > 100 && lambda <= 500
%   a = 30/100;
% else
%   a = 20/100;
% end
% 
% if L0_inf_percentage >= a || L1_inf_percentage >= a || (L0_inf_percentage >= a && L1_inf_percentage >= a)
%   warning('\nMore than %d percent of the examples are Inf ... Consider using larger regularization parameter\n',a*100)
%   keyboard;
% end
  
if 1
  % Set values of L=0 to 1e-20
  hh = h(theta_q,x_q);
  hh(hh==1) = 1-1e-10;
  L0 = I.' - hh;%h(theta_q,x_q);
%   L0(L0>=0 & L0<=1e-30) = 1e-30;
  L0 = real(log(L0));
  
  L1 = hh;%h(theta_q,x_q);
%   L1(L1>=0 & L1<=1e-30) = 1e-30;
  L1 = real(log(L1));
elseif 0
  % Set values of L=-inf to -1e10
  L0 = I.' - h(theta_q,x_q);
  L0 = real(log(L0));
  
  L1 = h(theta_q,x_q);
  L1 = real(log(L1));
  
%   L0(isinf(L0)) = -1e20;
%   L1(isinf(L1)) = -1e20;
  
  L0_inf_idxs = find(isinf(L0));
  L0(L0_inf_idxs) = sign(L0(L0_inf_idxs))*(1e20);
  
  L1_inf_idxs = find(isinf(L1));
  L1(L1_inf_idxs) = sign(L1(L1_inf_idxs))*(1e20);
elseif 0
  % Ignore examples that have L=-inf (i.e. set L=-inf to 0)
  L0 = I.' - h(theta_q,x_q);
  L0 = real(log(L0));
  
  L1 = h(theta_q,x_q);
  L1 = real(log(L1));
  
  L0(isinf(L0)) = 0;
  L1(isinf(L1)) = 0;
end

J = -1/Nexamples * sum(y_q.' .* L1 + (I-y_q).' .* L0) + reg_term;

x_q = x_q.';
grad(1) = (1/Nexamples) * sum((hh.'-y_q).*x_q(:, 1));

for idx = 2:length(theta_q)
  grad(idx) = (1/Nexamples) * sum((hh.'-y_q).*x_q(:, idx)) + (lambda * theta_q(idx) / Nexamples);
end

grad = grad.';
% if isnan(J)
%   keyboard
% end

% if isinf(J)
%   J = 1e10;
% end
return

