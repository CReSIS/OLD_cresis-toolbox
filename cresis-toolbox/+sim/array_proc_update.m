function [tomo_update] = array_proc_update(din,array_param,bin_idx,line_idx,doa)
% [tomo_update] = array_proc_update(din,array_param,bin_idx,line_idx,doa)
%
% Reruns array_proc for a single pixel with a new doa
%
% Todo: Need to support constrained searches on the DOA? For example, you
% change one source's DOA, but then want to search for the best location
% for the other DOA's.
%
% din: same input data as array_proc
% array_param: array_param structure returned from array_proc
% bin_idx,line_idx: specifies which pixel to update
% doa: new set of source DOA's to use for this pixel
% 
% Author: John Paden
%
% See also: run_array_proc_update.m

physical_constants;

rline_rng = array_param.rline_rng;

doa_param.fs              = array_param.wfs.fs;
doa_param.fc              = array_param.wfs.fc;
doa_param.Nsig            = array_param.Nsig;
doa_param.options         = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);
doa_param.doa_constraints = array_param.doa_constraints;
doa_param.theta_guard     = array_param.theta_guard;
doa_param.search_type     = array_param.init;

bin = array_param.bins(bin_idx);
line = array_param.lines(line_idx);

%% Steering Vector setup
for ml_idx = 1:length(array_param.imgs)
  % Make column vectors of y and z-positions
  for wf_adc_idx = 1:length(array_param.fcs{ml_idx})
    y_pos{ml_idx}(wf_adc_idx,1) = array_param.fcs{ml_idx}{wf_adc_idx}.pos(2,line);
    z_pos{ml_idx}(wf_adc_idx,1) = array_param.fcs{ml_idx}{wf_adc_idx}.pos(3,line);
  end
end
% Determine Steering Vector for beam-forming approaches
for ml_idx = 1:length(array_param.imgs)
  % Make column vectors of y and z-positions
  [~,array_param.sv{ml_idx}] = array_param.sv_fh(array_param.Nsv,array_param.wfs.fc,y_pos{ml_idx},z_pos{ml_idx});
end

doa_param.y_pc  = y_pos{1};
doa_param.z_pc  = z_pos{1};
doa_param.theta = array_param.theta;
doa_param.SV    = fftshift(array_param.sv{1},2);

% Number of fast-time samples in the din
Nt = size(din{1},1);

% Number of cross-track channels in the din
Nc = size(din{1},5);

% Number of subapertures
Na = size(din{1},3);

% Number of subbands
Nb = size(din{1},4);

%% MLE algorithm
switch array_param.method
  case 7
    
    dataSample  = din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:);
    dataSample  = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
    array_data  = dataSample.';
    Rxx         = (1/size(array_data,2)) * (array_data * array_data');
    doa_param.Rxx = Rxx;
    
    Jval = mle_cost_function(doa,doa_param);
    
    % Apply pseudoinverse and estimate power for each source
    k               = 4*pi*doa_param.fc/c;
    A               = exp(1i*k*(doa_param.y_pc*sin(doa(:)).' - doa_param.z_pc*cos(doa(:)).'));
    Weights         = inv(A'*A)*A';
    S_hat           = Weights*array_data;
    P_hat           = mean(abs(S_hat).^2,2);
    warning on;
    
    tomo_update.doa = doa;
    tomo_update.power = P_hat.';
    tomo_update.cost = Jval;
end

return;
