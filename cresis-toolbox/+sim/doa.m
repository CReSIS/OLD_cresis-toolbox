function [results, DCM_runs] = doa(param,param_override)
% results = doa(param,param_override)
%
% Function for simulating direction of arrival algorithms
%
% Author: John Paden, Theresa Stumpf, Gordon Ariho

%% General Setup

% =====================================================================
param = merge_structs(param, param_override);

fprintf('=====================================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.cluster.type, datestr(now));
fprintf('=====================================================================\n');

%% Input checks
% =========================================================================

if ~isfield(param.src,'ft_wind') || isempty(param.src.ft_wind)
  param.src.ft_wind = @boxcar;
end

%% Setup
% =========================================================================

% Get the list of array processing methods
%    array_proc_methods;

% param.src.fc: center frequency
param.src.fc = (param.src.f0 + param.src.f1)/2;
% param.src.fs: sampling frequency
param.src.fs = abs(param.src.f1 - param.src.f0);

% Get phase center positions of antenna array
if isfield(param.src,'phase_center')
  param.src.y_pc   = param.src.phase_center(2,:).';
  param.src.z_pc   = param.src.phase_center(3,:).';
  
elseif isfield(param.src,'lever_arm') && isfield(param.src.lever_arm,'fh')
  [phase_center] = param.src.lever_arm.fh(param.src.lever_arm.args{:});
  param.src.y_pc   = phase_center(2,:).';
  param.src.z_pc   = phase_center(3,:).';
elseif isfield(param.src,'y_pc') && isfield(param.src,'z_pc')
  % Do nothing. Already defined
else
  warning('Phase center information are not defined')
  keyboard
end

% Allocate all the output matrices
theta_est = [];
if isfield(param,'M') && ~isempty(param.M)
  % Model order estimation simulation. M=Nc-1 usually.
  for method = param.method.list
    theta_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),param.Nc-1);
    hessian_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),param.Nc-1);
    amp_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),param.Nc-1);
    cost_func{method} = zeros(param.monte.runs,size(param.monte.SNR,1));
  end
else
  for method = param.method.list
    theta_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),size(param.monte.SNR,2));
    hessian_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),size(param.monte.SNR,2));
    amp_est{method} = zeros(param.monte.runs,size(param.monte.SNR,1),size(param.monte.SNR,2));
    cost_func{method} = zeros(param.monte.runs,size(param.monte.SNR,1));
  end
end


% Setup the non-linear constraint function handle (this keeps theta values
% from getting too close and causing matrix singularity problems with the
% steering vector matrix)
% doa_nonlcon_fh = eval(sprintf('@(x) doa_nonlcon(x,%f);', param.method.theta_guard));

%% Cluster set up
% =========================================================================

global gRadar;

param = merge_structs(gRadar, param);
ctrl = cluster_new_batch(param);

argsin = param;
num_args_out = 4;
ctrl.cluster.max_mem_mode = 'auto';
cluster_compile({'sim_doa_task.m'},ctrl.cluster.hidden_depend_funs,ctrl.cluster.force_compile,ctrl);
dparam = [];
dparam.task_function = 'sim_doa_task';
dparam.argsin{1} = argsin;
block = true;
sparam.task_function = 'sim_doa_task';
sparam.argsin{1} = argsin;
sparam.num_args_out = num_args_out;
sparam.cpu_time = param.cpu_time;

for test_idx = 1:size(param.monte.SNR,1)
  param.test_idx = test_idx;
  dparam.argsin{1} = param;
  ctrl = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
 end
ctrl = cluster_save_dparam(ctrl);
if block
  % Blocking, submit and wait
  fprintf('Submitting %s\n', ctrl.batch_dir);
  ctrl_chain = {{ctrl}};
  ctrl_chain = cluster_set_chain(ctrl_chain,'cluster.cpu_time_mult',2);
  cluster_save_chain(ctrl_chain); 
  ctrl_chain = cluster_run(ctrl_chain,block);
  [in,out] = cluster_print(ctrl_chain{1}{1}.batch_id,[1:size(param.monte.SNR,1)],0);
  %cluster_cleanup(ctrl);
else
  % Non-blocking, do not submit
  out = ctrl;
end

% combining job outputs
for idx = 1:size(param.monte.SNR,1) %
  theta_est_1(:,idx,:) = out{idx}.argsout{1}{65539}(:,end,:);
  theta_est_2(:,idx,:) = out{idx}.argsout{1}{end}(:,end,:);
  
  hessian_est_1(:,idx,:) = out{idx}.argsout{3}{65539}(:,end,:);
  hessian_est_2(:,idx,:) = out{idx}.argsout{3}{end}(:,end,:);
  
  cost_func_1(:,idx) = out{idx}.argsout{4}{65539}(:,end);
  cost_func_2(:,idx) = out{idx}.argsout{4}{end}(:,end);
  
end
% populating the output with combined job outputs
A = {[]};
theta_est = repmat(A,1,65543);
theta_est{65539} = theta_est_1;
theta_est{end} = theta_est_2;

hessian_est = repmat(A,1,65543);
hessian_est{65539} = hessian_est_1;
hessian_est{end} = hessian_est_2;

cost_func = repmat(A,1,65543);
cost_func{65539} = cost_func_1;
cost_func{end} = cost_func_2;

rng_args = out{1}.argsout{2};
%
%% Test Loop: loop over each test
% =========================================================================
% Each test is a particular SNR, DOA, or snapshot configuration
%
% Copy outputs into output argument structure
if isfield(param,'Nsig_tmp') && ~isempty(param.Nsig_tmp)
  % For model order estimation simulation.
  results.theta_est = theta_est;
  results.rng_args = rng_args;
  
  if param.doa_example == 1
    results.hessian_est = hessian_est;
    results.cost_func = cost_func;
  end
  results.DCM = DCM_runs;
else
  results.theta_est = theta_est;
  results.rng_args = rng_args;
  results.hessian_est = hessian_est;
  results.cost_func = cost_func;
end
end

