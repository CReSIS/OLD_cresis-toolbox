function [ra,bw,L] = sigma_beamwidth_sar_length(data_in,input_methods,param)
% function [ra,bw,L] = sigma_beamwidth_sar_length(data_in,input_methods)
%
% Input either desired along-track resolution, full-beamwidth in air, or
% SAR aperture length and this will calculate the other two values.
%
% INPUTS:
%   data_in: Desired value of along-track resolution, beamwidth, or SAR
%            aperture length
%   input_methods: Type of data_in (resolution, beamwidth, or length)
%   param.f = Signal frequency in Hz
%   param.h = Height of platform above media interface
%   param.d = Depth for desired along-track resolution
%   param.n0 = Permittivity of first medium
%   param.n1 = Permittivity of second medium
%
% Example:
%   [ra,bw,L] = sigma_beamwidth_sar_length([100],{'sar_length'});
%   ra =
%       6.0143
%   bw =
%       7.3232
%   L =
%       100
%
% Author: Logan Smith

error(nargchk(2,3,nargin,'struct'));
physical_constants

if ~isfield(param,'f')
  param.f = 195e6;
end
f = param.f;
if ~isfield(param,'h')
  param.h = 500;
end
h = param.h;
if ~isfield(param,'d')
  param.d = 500;
end
d = param.d;
if ~isfield(param,'n0')
  param.n0=1;
end
n0 = param.n0;
if ~isfield(param,'n1')
  param.n1=sqrt(er_ice);
end
n1 = param.n1;
lambda_air = c/f;
lambda_ice = c/n1/f;

% input_methods = {'sigma','beamwidth','sar_length'}; % 'sigma', 'beamwidth', 'sar_length'
for method_idx=1:length(input_methods)
  input_method = input_methods{method_idx};
  if strcmp(input_method,'sigma')
    %     ra = [1.046	2.055	4.093	8.178	16.351];
    if isvector(data_in)
      ra = data_in;
    else
      ra = data_in(:,method_idx);
    end
    
    %% Homogeneous medium
    % Calculate the SAR aperture length if the radar was at the surface.
    L_ra = lambda_ice./(2.*asin(ra./d));
    
    %% Two-media model
    x1 = L_ra/2;
    theta1 = atan(x1/d);
    theta0 = asin(n1/n0.*sin(theta1));
    x0 = h*tan(theta0);
    L = 2*(x0+x1);
    bw = 2*theta0*180/pi;
    
    fprintf('Sigma_x (m: Along-track Resolution) - INPUT\n')
    fprintf('%5.3f\t',ra(1:end-1))
    fprintf('%5.3f\n',ra(end))
    fprintf('Full Beamwdith (deg)\n')
    fprintf('%5.3f\t',bw(1:end-1))
    fprintf('%5.3f\n',bw(end))
    fprintf('SAR Length (m)\n')
    fprintf('%5.2f\t',L(1:end-1))
    fprintf('%5.2f\n',L(end))
  elseif strcmp(input_method,'beamwidth')
    %     full_beamwidth = [42.192	21.436	10.761	5.386	2.694]; % deg
    if isvector(data_in)
      full_beamwidth = data_in;
      bw = data_in;
    else
      full_beamwidth = data_in(:,method_idx);
      bw = data_in(:,method_idx);
    end
    theta0 = full_beamwidth/2*pi/180; % rad
    
    %% Two-media model
    x0 = h.*tan(theta0);
    theta1 = asin(n0./n1.*sin(theta0));
    x1 = d.*tan(theta1);
    L = 2*(x0+x1);
    
    L_ra = 2*x1;
    
    ra = sin(lambda_ice./(2.*L_ra)).*d;
    fprintf('Sigma_x (m: Along-track Resolution)\n')
    fprintf('%5.3f\t',ra(1:end-1))
    fprintf('%5.3f\n',ra(end))
    fprintf('Full Beamwdith (deg) - INPUT\n')
    fprintf('%5.3f\t',bw(1:end-1))
    fprintf('%5.3f\n',bw(end))
    fprintf('SAR Length (m)\n')
    fprintf('%5.2f\t',L(1:end-1))
    fprintf('%5.2f\n',L(end))
    
  elseif strcmp(input_method,'sar_length')
    %     L = [800 400 200 100 50];
    if isvector(data_in)
      L = data_in;
    else
      L = data_in(:,method_idx);
    end
    for idx=1:length(L)
      x = L(idx)/2;
      alpha = 0;
      
      [x0,x1] = two_media_model(x,d,h,er_ice,alpha);
      L_ra(idx) = 2*x1;
      
      theta0(idx) = atan(x0/h);
    end
    bw = 2.*theta0*180/pi;
    
    ra = sin(lambda_ice./(2.*L_ra)).*d;
    fprintf('Sigma_x (m: Along-track Resolution)\n')
    fprintf('%5.3f\t',ra(1:end-1))
    fprintf('%5.3f\n',ra(end))
    fprintf('Full Beamwdith (deg)\n')
    fprintf('%5.3f\t',bw(1:end-1))
    fprintf('%5.3f\n',bw(end))
    fprintf('SAR Length (m) - INPUT\n')
    fprintf('%5.2f\t',L(1:end-1))
    fprintf('%5.2f\n',L(end))
  end
  fprintf('\n')
end

return

