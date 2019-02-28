function [sources_number, log_penalty_cost_subopt] = model_order_suboptimal(model_order_suboptimal_param)
% Function for estimating the model order or number of targets using suboptimal methods.
Nc                 = model_order_suboptimal_param.Nc;
Nsnap              = model_order_suboptimal_param.Nsnap;
eigval             = model_order_suboptimal_param.eigval;
model_order_method = model_order_suboptimal_param.method;
penalty_NT         = model_order_suboptimal_param.penalty_NT;
param_MOE          = model_order_suboptimal_param.param_MOE;

M = length(penalty_NT)-1;
for k = 0:M
    geom_mean   = sum(log(eigval(k+1:end)));
    arithm_mean = (Nc-k)*(log(sum((eigval(k+1:end))/(Nc-k))));
                       
    log_lq(k+1) = Nsnap*(geom_mean- arithm_mean);
    DOF(k+1) = k*((2*Nc)-k) + 1;  % degree of freedom of the space spanned by signal vectors
end

if param_MOE.norm_allign_zero ==1
    log_func_subopt = -2*log_lq + param_MOE.norm_term_suboptimal;
else
    log_func_subopt = -2*log_lq;
end

switch model_order_method
    case 0
        % SRAVYA penalty
        penalty = penalty_NT;
    case 1
        % Akaike information criterion
        penalty = 2*DOF;
    case 2
        % Hannan and Quinn criterion
        penalty =  2*DOF*log(log(Nsnap));
    case 3
        % MDL
        % Minimum descriptive length (MDL)
        % Bayesian information criterion
        penalty =  DOF*log(Nsnap);
    case 4
        % Corrected AIC
        %eq 9 --- results in paper matches with this eq
        penalty =  (2*Nsnap.*DOF)./(Nsnap-DOF-1);
        
        %eq 19
        %penalty =  (2*Nsnap.*DOF)./(Nsnap-Nc-[0:M]-1);% gives very good results
    case 5
        q = 0:M;
        % Vector corrected Kullback information criterion
        penalty =  (Nsnap*Nc*(2.*q+Nc+1))./(Nsnap-Nc-q-1) + (Nsnap*Nc)./(Nsnap-q-((Nc-1)/2))...
            + (2*Nc.*q + Nc*Nc - Nc)./2 ;
    case 6
        % Weighted-average information criterion
        %         A_tmp = (2*Nsnap.*DOF) ./ (Nsnap-DOF-1);
        %         B_tmp = DOF.*log(Nsnap);
        %         penalty = (A_tmp.^2 + B_tmp.^2) / (A_tmp + B_tmp);
        
        penalty = ((2*Nsnap*DOF).^2 + ((Nsnap-DOF-3).*DOF.*log(Nsnap)).^2 )...
            ./ (2*Nsnap*DOF.*(Nsnap-DOF-3) + (Nsnap-DOF-3).^2 .* DOF*log(Nsnap));
    otherwise
        error('Not supported')
end

cost = log_func_subopt + penalty;

log_penalty_cost_subopt = [log_func_subopt penalty cost];
% log_penalty_cost_subopt(model_order_method+1,:) = [log_func_subopt penalty cost];

%         figure(2);
%         subplot(2,4,model_order_method)
%         plot(0:M ,cost,'b*-')
%         hold on
%         plot(0:M ,-2*log_lq,'g*-')
%         hold on
%         plot(0:M ,penalty,'-r*')
%         legend('MDL cost','log-likelihood','penalty' ,'location', 'northeast')
%         set(gca,'xtick',0:M)
%         title( method_title(model_order_method,:));
%         ylabel('Cost');
%         xlabel('Estimated number of sources (k)');
%         hold off
%         grid on


[~, index] = min(cost);
sources_number = index-1;

%
%         figure(2);
%         plot(0:M ,cost,'b*-')
%         hold on
%         plot(0:M ,-2*log_lq,'g*-')
%         hold on
%         plot(0:M ,penalty,'-r*')
%         legend('MDL cost','log-likelihood','penalty' ,'location', 'east')
%         set(gca,'xtick',0:M)
%         ylabel('Cost');
%         xlabel('Estimated number of sources (k)');
%         grid on
% hold off
return