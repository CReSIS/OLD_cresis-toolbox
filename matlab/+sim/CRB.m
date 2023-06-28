function [CRB_angular , CRB_spatial] = CRB(param,desired_src)
%This CRB function determines Cramer-Rao bound of MUSIC, MLE, WBMLE, and WBDCM.
physical_constants;
[phase_center]            = param.src.lever_arm.fh(param.src.lever_arm.args{:});
param.src.y_pc            = phase_center(2,:).';
param.src.z_pc            = phase_center(3,:).';
param.src.fc              = (param.src.f0 + param.src.f1)/2;
param.src.fs              = param.src.f1 - param.src.f0;
param.src.widening_factor = param.method.wb_td.widening_factor;
param.src.ft_wind         = @boxcar;

% w       = param.src.widening_factor;
fc      = param.src.fc;
fs      = param.src.fs;
yAnt    = param.src.y_pc;
zAnt    = param.src.z_pc;
sigma_n = 10^(0/10);

if ~isfield(param.method,{'list'})
    param.method.list = [2 7 8 9];
end

if sum(((param.method.list == 8) | (param.method.list == 9)))
    % WB method and WBMLE method
    Nk = param.method.wb_fd.filter_banks;
else
    % MUSIC method and MLE method
    Nk = 1;
end

k = 4*pi*(fc + fs*[0:floor((Nk-1)/2), -floor(Nk/2):-1]/Nk)/c;
d = abs(yAnt(end)-yAnt(end-1));  %inter-element spacing

for test_idx = 1:size(param.monte.SNR,1)
    
    param.src.SNR    = param.monte.SNR(test_idx,:);
    param.src.Nsnap  = param.monte.Nsnap(test_idx);
    param.src.DOAs   = param.monte.DOA(test_idx,:);
    %     [Data,DCM,imp_response,DCM_fd] = sim.doa_wideband_data(param);
    Nsrc             = numel(param.src.DOAs);
    Z_exact_Doron    = zeros(Nsrc);
    
    SNR              = 10.^(param.src.SNR./10);
    Rss = repmat(SNR',[1 Nsrc]) .* eye(Nsrc);
    
    for k_idx = 1:Nk
        A     = sqrt(1/length(yAnt)) * exp(1i*k(k_idx)*(-zAnt*cosd(param.src.DOAs) + yAnt*sind(param.src.DOAs)));
        D     = repmat(yAnt*k(k_idx),[1 Nsrc]);
        A_dot = (1i*D).*A;
        
        proj_mtx = A * inv(A'*A) * A';
        I        = eye(size(proj_mtx));
        
        g = 1./(k(k_idx)*d*cosd(param.src.DOAs)');
        G = diag(g);
        
        %         Rxx = DCM_fd((k_idx-1)*size(DCM_fd,2)+(1:size(DCM_fd,2)),:);
        %         Rss = inv(A'*A)*A'*(Rxx-sigma_n*eye(size(Rxx)))*A*inv(A'*A);
        
        H_exact_Doron = (Rss - inv(inv(Rss) + (1/sigma_n) * (A' * A))).';
        Y             = A_dot'*(I-proj_mtx)*A_dot;
        Z_exact_Doron = Z_exact_Doron + (real(Y.*H_exact_Doron))/sigma_n;
    end
    
    CRB_spatial_mtx_exact_Doron         = (inv(Z_exact_Doron))/(2*param.src.Nsnap);
    CRB_angular_mtx_exact_Doron         = G*CRB_spatial_mtx_exact_Doron*G;
    CRB_spatial_exact_Doron(:,test_idx) =  diag(CRB_spatial_mtx_exact_Doron);
    CRB_angular_exact_Doron(:,test_idx) = diag(CRB_angular_mtx_exact_Doron);
    
end
if ~exist('desired_src','var') || isempty(desired_src)
    CRB_spatial =  CRB_spatial_exact_Doron;
    CRB_angular =  CRB_angular_exact_Doron;
else
    if (desired_src > Nsrc) || (desired_src <= 0)
        error('Error using CRB function: desired_src must be a real integer between 1 and number of sources');
    end
    CRB_spatial =  CRB_spatial_exact_Doron(desired_src,:);
    CRB_angular =  CRB_angular_exact_Doron(desired_src,:);
end

return;