function [cost_overall,data] = optimizer_fit_NT(x0,optim_fit_param)
%P_md_coeff is a scalar to control Probability of missed detection, P_md, (underestimation)
%If P_md is very low, then false alarms will be higherï¿? i.e. weï¿½ll overestimate a lot.

%alternately the probability of false alarm, P_fa, can be set (overestimation),
%but not both P_md and P_fa together since they depend on one another.
param   = optim_fit_param.param;
results = optim_fit_param.results;
P_md_coeff = 1;
P_fa_coeff = 1;

N_snap = param.monte.Nsnap ;
SNR_db = param.SNR_db  ;  % CODE IT
runs = param.monte.runs;
M = param.M;

x0_intial = param.x0_intial;
x0_diff   = x0;

clear x0
for idx=1:M  %summing the differences to form penalty term
    x0(idx)= sum(x0_diff(1:idx));
end

x0 = x0_intial(1)+[0 x0];

for source_idx = 0:M    
    result = results{source_idx+1};
 
    log_func_all =     result.log_func_all ;    
    q_actual = source_idx;
    
    for N_idx=1:length(N_snap)
        for SNR_idx=1:length(SNR_db)
            SNR = SNR_db (SNR_idx);

            log_cost_temp = log_func_all{SNR_idx,N_idx};
            
            for runs_idx = 1:runs
                
%                 for k = 1:p-1
%                     
%                     penality_new(k) = x0(1)*(k^x0(2))*((200* (SNR/10)^x0(3)) + x0(4));
%                     
%                     penality_new(k) =x0(1)*(k^x0(2))*((x0(3)* (SNR/10)) + x0(4));
                penality_new = x0;
%                 
%             end
%             
%             penality_runs_new(runs_idx,:) = [0 penality_new];
%                 penality_runs_new(runs_idx,:) = penality_new;
                cost_runs_new(runs_idx,:) = log_cost_temp(runs_idx,:) + penality_new;%penality_runs_new(runs_idx,:);
                
                [~, index] = min(cost_runs_new(runs_idx,:));
                sources_number_new(SNR_idx,N_idx,runs_idx) = index-1;
            end
            
%             penality_all_new{SNR_idx,N_idx} = penality_runs_new;
%             cost_all_new{SNR_idx,N_idx} = cost_runs_new;
            
            % mean of runs result
%             cost_mean_new{SNR_idx,N_idx} = mean(cost_all_new{SNR_idx,N_idx});
            
%             [~, index] = min(mean(cost_all_new{SNR_idx,N_idx}));
%             q_estimated_mean(SNR_idx,N_idx) = index-1;
            %
            % subplot(length(SNR_db),2*length(N_snap),fig)
            % % errorbar(0:4 ,mean(log_func_all{SNR_idx,N_idx}),std(log_func_all{SNR_idx,N_idx}) )
            % plot(0:4 ,mean(log_func_all{SNR_idx,N_idx}));
            % hold on
            % plot(0:4 ,mean(penality_all_new{SNR_idx,N_idx}),'-r*')
            % hold on
            % plot(0:4 ,mean(cost_all_new{SNR_idx,N_idx}),'-b*')
            %
            %
            % title([ 'SNR db = ' num2str(SNR_db(SNR_idx ))''  , '   N = ' num2str(N_snap (N_idx))'' ])
            % ylabel('Cost');
            % xlabel('Estimated number of sources');
            % xlim([1 4])
            % % ylim([0 18e5])
            % grid on
            %
            % fig = fig + 1;
            % subplot(length(SNR_db),2*length(N_snap),fig)
            % histogram(sources_number_new(SNR_idx,N_idx,:))
            % xlim([-0.5 4.5])
            % ylim([0 100])
            % grid on
        end
    end
    
%     sources_number_all{source_idx+1} = sources_number_new;
    
    source_q_actual_new = find (sources_number_new == q_actual);
    q_times(source_idx+1) = length(source_q_actual_new); % For all SNRS and all runs at once.
    
    idx = 1:numel(sources_number_new);
    % COST_UPDATE_func(idx) = (sources_number_new(idx)- q_actual ).^2;
    % COST_UPDATE_sum(source_idx+1) = sum(COST_UPDATE_func);
    
    COST_UPDATE_func(idx) = (sources_number_new(idx) - q_actual);
    COST_UPDATE_func_left = (COST_UPDATE_func(find(COST_UPDATE_func < 0))).^2;
    COST_UPDATE_func_right = (COST_UPDATE_func(find(COST_UPDATE_func > 0))).^2;
    COST_UPDATE_sum(source_idx+1) = P_md_coeff*sum(COST_UPDATE_func_left) + P_fa_coeff*sum(COST_UPDATE_func_right);
end

tests= size(log_func_all,1)* size(log_func_all{1,1},1);

if 1
    
    cost_overall = sum(COST_UPDATE_sum)+sum(tests*ones(size(q_times))-q_times)*10;
    
end
return