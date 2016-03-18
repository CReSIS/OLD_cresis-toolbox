% Estimate the number of sources
%
%
N = 11*3;
Q = 1; % Number of sources
SNR = 15;
p = 8;


param.src.Nsnap = N;
param.src.SNR = SNR*ones(1,Q);
param.src.DOAs = 20*(0:Q-1);
param.src.fs = 2e6;
param.src.f0 = 194e6;
param.src.f1 = 196e6;
param.src.fc = 195e6;
physical_constants;
lambda_fc = c/param.src.fc;
param.src.y_pc = lambda_fc/4 * (0:p-1).';
param.src.z_pc = 0*(0:p-1).';
param.src.widening_factor = 1;
param.src.ft_wind = @hanning;
param.src.noise_en = true;
param.method.wb_td.widening_factor = 1;
param.method.wb_fd.filter_banks = 1;

num_runs = 1000;
q_est_aic = zeros(1,num_runs);
q_est_dml = zeros(1,num_runs);
for run = 1:num_runs

[Data, DCM, imp_resp, DCM_fd] = sim.doa_wideband_data(param);

eigen_values = sort(eig(DCM).','descend');

% N: number of snapshots
% p: number of sensors
% k: estimated number of sources
% D: eigenvalues

  
AIC = [];
MDL = [];
k_vals = 0:p-1;
for k = k_vals
%   AIC(k+1) = -log( (prod(eigen_values(k+1:p).^(1/(p-k))) ...
%     ./ mean(eigen_values(k+1:p))).^((p-k)*N) ) + k*(2*p-k);
  AIC(k+1) = -log( (prod(eigen_values(k+1:p).^(1/(p-k))) ...
    ./ mean(eigen_values(k+1:p))) ) * ((p-k)*N) + k*(2*p-k);
  
%   MDL(k+1) = -log( (prod(eigen_values(k+1:p).^(1/(p-k))) ...
%     ./ mean(eigen_values(k+1:p))).^((p-k)*N) ) + 0.5*k*(2*p-k)*log(N);
  MDL(k+1) = -log( (prod(eigen_values(k+1:p).^(1/(p-k))) ...
    ./ mean(eigen_values(k+1:p))) ) * ((p-k)*N) + k*(2*p-k)*(0.5*log(N));
end

abs(AIC);
[~,q_est] = min(AIC);
q_est = k_vals(q_est);
q_est_aic(run) = q_est;

abs(MDL);
[~,q_est] = min(MDL);
q_est = k_vals(q_est);

q_est_dml(run) = q_est;

end

figure(1); clf;
subplot(2,1,1);
res = hist(q_est_dml,[0:p])
bar(0:p,res/num_runs*100);
grid on;
xlim([-0.5 p+0.5]);
title('DML');
ylabel('Percentage');
subplot(2,1,2);
res = hist(q_est_aic,[0:p])
bar(0:p,res/num_runs*100);
grid on;
xlim([-0.5 p+0.5]);
title('AIC');
ylabel('Percentage');
xlabel('Estimated number of sources');
