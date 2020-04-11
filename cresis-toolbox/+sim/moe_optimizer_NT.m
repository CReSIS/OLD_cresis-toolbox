function [penalty_coeff, x0_all,param_debug_NT] = optimizer_NT(param)
% Generates Nt penalty coefficients
%
physical_constants;
N_snap = param.monte.Nsnap ;
M      = param.M;
%  p: number of elements on the array(number of sensors)
p  = param.Nc;
fc = (param.src.f0 +param.src.f1)/2;

% lambda: wavelength
lambda = c/fc;

d= lambda/2;
L =p*d;

k_spacing = [-(p-1)/2:1:(p-1)/2]*(lambda/L); % as implemented in paper for orthogonal sv
%k_spacing = [-p/2:1:(p/2)-1]*(lambda/L); % John

DOA_orthogonal = asin(k_spacing)*180/pi;
index_ref = ceil(p/2);

for source_idx = 0:M
    if   source_idx == 0
        q = [];
    else
        if mod(source_idx,2)==1 %odd
            index =((index_ref)-(source_idx-1)/2):1:(index_ref)+(source_idx-1)/2;
        else
            index =((index_ref)-(source_idx)/2):1:(index_ref)+(source_idx)/2;
            index(find(index==index_ref))= [];
        end
        q = DOA_orthogonal(index);
        clear index
    end
    
    if param.opt == 0
        result = sim.doa_example_NT(q,param);
        
    elseif   param.opt == 1
        result = sim.doa_example_NT_opt(q,param);
    end
    log_func_all =     result.log_func_all ;
    results{source_idx+1} = result;
end
n = N_snap;

q= [0:M];

m = q.*((2*p)-q) + 1;  % degree of freedom of the space spanned by signal vectors

AIC =  2*m;

MDL= (1/2)*2.*m*log(n);

HQ=2*m*log(log(n));

AICc = (2*n.*m)./(n-m-1);       %eq 9 --- results in paper matches with this eq

%AICc =(2*n.*m)./(n-p-q-1); % eq 19--- gives very good results

KICvc= (n*p*(2.*q+p+1))./(n-p-q-1) + (n*p)./(n-q-((p-1)/2))...
    + (2*p.*q + p*p - p)./2 ;

% A = 2*n*m ./ (n-m-1);
% B = m.*log(n);
% WIC = (A.^2 + B.^2) ./ (A + B);

WIC = ((2*n*m).^2 + ((n-m-3).*m.*log(n)).^2 )...
    ./ (2*n*m.*(n-m-3) + (n-m-3).^2 .* m*log(n));

x0_all = [AIC;HQ;MDL;AICc;KICvc;WIC ];
% save penalty_all x0_all

if param.opt == 0
    % if N_snap <= 15
    % x0_intial = MDL ;
    % else
    %   x0_intial = WIC ;
    % end
    
    x0_intial = MDL ; % I THINK PREVIOUS RESULTS ARE INITIALIZING WITH AICc
end

if param.opt == 1
    % if N_snap <= 15
    % x0_intial = MDL ;
    % else
    %   x0_intial = WIC ;
    % end
    x0_intial = MDL ;
    
end
x0_intial = x0_intial(1:M+1);

param.x0_intial = x0_intial;

penalty_coeff = diff(x0_intial) ;
if 0  %used when optimal methods are not normalized
    
    Aieq = -1*(horzcat(-1*eye(Nc-2),zeros(Nc-2,1))+    horzcat(zeros(Nc-2,1),eye(Nc-2)));
    bieq = [zeros(Nc-2,1)];
    
    % options =  psoptimset('TolX',1,'InitialMeshSize',150 ,'TolMesh',1,'MeshContraction',0.75)
    %  options =  psoptimset('TolX',0.5,'InitialMeshSize',200 ,'TolMesh',0.5,'MeshContraction',0.75)
    %
    %   coeff = patternsearch(@optimizer_fit_NT,1000*ones(length(x0),1),Aieq,bieq,[],[],x0,1500*ones(length(x0),1),[],options);
    options =  psoptimset('TolX',0.5,'InitialMeshSize',50 ,'TolMesh',0.5,'MeshContraction',0.75)
    
    coeff = patternsearch(@optimizer_fit_NT,[200 200 100 200 700 700],[],[],[],[],[200 200 100 200 700 700],[300 300 300 300 1500 1500],[],options);
end
if param.opt == 0  %SUBOPOTIMAL
    param_debug_NT = results;
    
    Aieq = [];%1*(horzcat(-1*eye(M-1),zeros(M-1,1))+    horzcat(zeros(M-1,1),eye(M-1)))
    bieq = [];%[zeros(M-1,1)]
    UB = []; %150*ones(length(x0),1)
    LB = penalty_coeff;
    options =  psoptimset('TolX',1,'InitialMeshSize',20 ,'TolMesh',1,'MeshContraction',0.75);
    
    optim_fit_param.results = results;
    optim_fit_param.param   = param;
    
    coeff = patternsearch(@(x)optimizer_fit_NT(x,optim_fit_param),penalty_coeff,Aieq,bieq,[],[],LB,UB,[],options);
    
elseif param.opt == 1  %OPTIMAL
    %   load params_debug_NT param_debug_NT
    param_debug_NT = results;
    %    save params_debug_NT param_debug_NT
    
    Aieq =  [];%1*(horzcat(-1*eye(M-1),zeros(M-1,1))+    horzcat(zeros(M-1,1),eye(M-1)))
    bieq =  [];%[zeros(M-1,1)]
    UB =[]; %300*ones(length(x0),1); %
    LB = penalty_coeff;
    options =  psoptimset('TolX',1,'InitialMeshSize',40 ,'TolMesh',1,'MeshContraction',0.75);
    
    optim_fit_param.results = results;
    optim_fit_param.param   = param;
    
    coeff = patternsearch(@(x)optimizer_fit_NT(x,optim_fit_param),penalty_coeff,Aieq,bieq,[],[],LB,UB,[],options);
end
%   x0 = coeff;

% x0_diff = coeff;

clear x0

for idx=1:M  %summing the differences to form penalty term
    penalty_coeff(idx)= sum(coeff(1:idx));
end
penalty_coeff = x0_intial(1)+[0 penalty_coeff];
return