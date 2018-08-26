function [sources_number,doa, log_penalty_cost_opt] = model_order_optimal(model_order_optimal_param)
% Function for estimating the model order or number of targets using optimal methods.
Nc                 = model_order_optimal_param.Nc;
Nsnap              = model_order_optimal_param.Nsnap;
eigval             = model_order_optimal_param.eigval;
eigvec             = model_order_optimal_param.eigvec;
model_order_method = model_order_optimal_param.method;
penalty_NT_opt         = model_order_optimal_param.penalty_NT_opt;
param_MOE          = model_order_optimal_param.param_MOE;
doa_mle            = model_order_optimal_param.doa_mle;

sv_param.src.y_pc  = model_order_optimal_param.y_pc;
sv_param.src.z_pc  = model_order_optimal_param.z_pc;
sv_param.src.fc    = model_order_optimal_param.fc;

M = length(doa_mle);
p = Nc;

clear MDL_cost log_func penalty

sum_term = [];
term1 = [];
for k = 0:M
    if k == 0
        proj = 0; %null matrix [];
    else
        A    = steering_mtx(doa_mle{k},sv_param);
        proj = A*((A'*A)\A');
    end
    I = eye(Nc);
    for p_idx = 1:Nc
        term1(k+1,p_idx) = eigval(p_idx).*norm((I-proj)*eigvec(:,p_idx))^2;
    end
    sum_term(k+1) = sum(term1(k+1,:));
    DOF(k+1) = k*((2*Nc)-k) + 1;  % degree of freedom of the space spanned by signal vectors
end
log_func = -2*(-Nsnap*p*(log((sum_term)./p)))  ;

%NORMALIZING OPTIMAL log-likelihood function
%(the normalization is to subtract away from all optimal results
%the difference between the suboptimal and optimal for q=0 case
%for SNR = 30 dB and N snapshots. This will make it so that the mean
%of the optimal case for q=0 and 30 dB SNR matches the suboptimal.
%We can then use the regular algorithms like AIC, MDL, etc. on the optimal results.)

if param_MOE.norm_allign_zero ==1
    log_func = log_func + param_MOE.norm_term_optimal;
else
    log_func = log_func + param_MOE.opt_norm_term;
end

switch model_order_method
    case 0
        % Numrically Tuned
        penalty = penalty_NT_opt;
    case 1
        % Akaike information criterion
        penalty = 2*DOF;
    case 2
        % Minimum descriptive length (MDL)
        
        % Bayesian information criterion
        penalty = (1/2)*2*DOF*log(Nsnap);
    case 3
        % Hannan and Quinn criterion
        penalty =  2*DOF*log(log(Nsnap));
    case 4
        % Corrected AIC
        %eq 9 --- results in paper matches with this eq
        penalty =  (2*Nsnap.*DOF)./(Nsnap-DOF-1);
        
        %eq 19
        %penalty =  (2*Nsnap.*DOF)./(Nsnap-Nc-[0:M]-1);% gives very good results
    case 5
        q = 0:M;
        % Vector corrected Kullback information criterion
        penalty = (Nsnap*Nc*(2.*q+Nc+1))./(Nsnap-Nc-q-1) + (Nsnap*Nc)./(Nsnap-q-((Nc-1)/2))...
            + (2*Nc.*q + Nc*Nc - Nc)./2 ;
    case 6
        % Weighted-average information criterion
        penalty = ((2*Nsnap*DOF).^2 + ((Nsnap-DOF-3).*DOF.*log(Nsnap)).^2 )...
            ./ (2*Nsnap*DOF.*(Nsnap-DOF-3) + (Nsnap-DOF-3).^2 .* DOF*log(Nsnap));
    otherwise
        error('Not supported')
end

cost = log_func + penalty;

[~, index] = min(cost);
sources_number = index-1;

doa = NaN *ones(1,M);

if sources_number > 0
    doa(1:sources_number) = doa_mle{sources_number};
end

log_penalty_cost_opt = [log_func penalty cost];
% log_penalty_cost_opt(model_order_method+1,:) = [log_func penalty cost];
return