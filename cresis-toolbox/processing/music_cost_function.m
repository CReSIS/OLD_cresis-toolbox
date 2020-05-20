function S = music_cost_function(theta, param)
% S = music_cost_function(theta,param)
%
% Function used to compute MUSIC pseudospectrum. Minimized in array_proc
% under method 9 (music_doa) in order to refine estimate of local peaks.
%
% Inputs:
%   theta = 1 x 1 scalar where pseudospectrum is to be evaluated
% 
%   param = control structure containing the following fields:
%     .Nsrc   = 1 x 1 scalar specifying dimensionality of signal subspace,
%     .Rxx    = Nc x Nc complex valued sample covariance matrix,
%     .fc     = 1 x 1 scalar valued carrier frequency (Hz)
%     .y_pc   = Nc x 1 vector containing y coordinates of phase centers
%               in SAR FCS (meters),
%     .z_pc   = Nc x 1 vector containing z coordinates of phase centers in
%               SAR FCS (meters).
% Outputs:
%   S = 1x1 scalar specifying value of MUSIC pseudospectrum at theta.
%
% Author:   Theresa Stumpf
% 
% See also: array_proc.m, music_initialization.m
% =========================================================================


%% Evaluate steering vectors at test value theta
% =========================================================================
%
% NOTE: Currently only ideal steering vectors are supported. This could be 
% improved by passing in LUT through the param structure.
c = 2.997924580003452e+08; % physical_constants too slow
k = 4*pi*param.fc*param.sv_dielectric/c;
ky = k*sin(theta);
kz = k*cos(theta);
ky = ky(:).';
kz = kz(:).';
SV = sqrt(1/length(param.y_pc))*exp(1i*(-param.z_pc*kz + param.y_pc*ky));

%% Compute Cost
% =========================================================================
[V,D]     = eig(param.Rxx);
eigenVals = diag(D);
[eigenVals, noiseIdxs] = sort(eigenVals);
noiseIdxs = noiseIdxs(1:end - param.Nsrc);
Qn        = V(:,noiseIdxs);
S         = mean(abs(SV(:,:)' * Qn).^2,2);
S         = sum(S);
end