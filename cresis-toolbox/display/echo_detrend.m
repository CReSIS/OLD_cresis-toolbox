function data = echo_detrend(mdata, param)
% data = echo_detrend(mdata, param)
%
% The trend of the data is estimated using various methods and this trend
% is removed from the data.
%
% INPUTS:
%
% data = 2D input data matrix (linear power)
%
% param: struct controlling how the detrending is done
%  .units: units of the data, string containing 'l' (linear power) or 'd' (dB
%  log power)
%
%  .method: string indicating the method. The methods each have their own
%  parameters
%
%   'local': takes a local mean to determine the trend
%     .filt_len: 2 element vector of positive integers indicating the size of
%     the trend estimation filter where the first element is the boxcar
%     filter length in fast-time (row) and the second element is the boxcar
%     filter length in slow-time (columns). Elements must be odd since
%     fir_dec used. Elements can be greater than the corresponding
%     dimensions Nt or Nx in which case an ordinary mean is used in that
%     dimension.
%
%   'mean': takes the mean in the along-track and averages this in the
%   cross-track. This is the default method.
%
%   'polynomial': polynomial fit to data between two layers, outside of the
%   region between the two layers uses nearest neighborhood interpolation
%     .layer_bottom: bottom layer
%     .layer_top: top layer
%     .order: nonnegative integer scalar indicating polynomial order
%
%   'tonemap': uses Matlab's tonemap command
%
% OUTPUTS:
%
% data: detrended input
%
% Examples:
%
% fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_standard/20140512_01/Data_20140512_01_018.mat';
% mdata = load(fn);
%
% imagesc(lp(echo_detrend(mdata)));
%
% [surface,bottom] = layerdata.load_layers(mdata,'','surface','bottom');
% imagesc(lp(echo_detrend(mdata, struct('method','polynomial','layer_top',surface,'layer_bottom',bottom))));
% 
% imagesc(lp(echo_detrend(mdata, struct('method','tonemap'))));
% 
% imagesc(lp(echo_detrend(mdata, struct('method','local','filt_len',[51 101]))));
%
% Author: John Paden

if isstruct(mdata)
  data = mdata.Data;
else
  data = mdata;
end

if ~exist('param','var') || isempty(param)
  param = [];
end

if ~isfield(param,'method') || isempty(param.method)
  param.method = 'mean';
end

switch param.method
  case 'local'
    if ~isfield(param,'filt_len') || isempty(param.filt_len)
      param.filt_len = [21 51];
    end
    Nt = size(data,1);
    Nx = size(data,2);
    if param.filt_len(1) < Nt
      trend = nan_fir_dec(data.',ones(1,param.filt_len(1))/param.filt_len(1),1).';
    else
      trend = repmat(nan_mean(data,1), [Nt 1]);
    end
    if param.filt_len(2) < Nx
      trend = nan_fir_dec(trend,ones(1,param.filt_len(2))/param.filt_len(2),1);
    else
      trend = repmat(nan_mean(trend,2), [Nx 1]);
    end
    data = data ./ trend;
    
  case 'mean'
    data = bsxfun(@times,data,1./mean(data,2));
    
  case 'polynomial'
    if ~isfield(param,'order') || isempty(param.order)
      param.order = 2;
    end
    if any(param.layer_bottom<1)
      % layer is two way travel time
      if isstruct(mdata)
        param.layer_bottom = round(interp_finite(interp1(mdata.Time,1:length(mdata.Time),param.layer_bottom,'linear','extrap'),0));
      end
    end
    if any(param.layer_top<1)
      % layer is two way travel time
      if isstruct(mdata)
        param.layer_top = round(interp_finite(interp1(mdata.Time,1:length(mdata.Time),param.layer_top,'linear','extrap'),0));
      end
    end
    
    Nt = size(data,1);
    Nx = size(data,2);
    mask = false(size(data));
    x_axis = nan(size(data));
    for rline = 1:Nx
      bins = max(1,min(Nt,param.layer_top(rline))):max(1,min(Nt,param.layer_bottom(rline)));
      mask(bins,rline) = true;
      x_axis(bins,rline) = (bins - param.layer_top(rline)) / (param.layer_bottom(rline) - param.layer_top(rline));
    end
    
    % Find the polynomial coefficients
    detrend_poly = polyfit(x_axis(mask),lp(data(mask)), param.order);
    
    trend = zeros(Nt,1);
    for rline = 1:Nx
      % Section 2: Evaluate the polynomial
      bins = max(1,min(Nt,param.layer_top(rline))):max(1,min(Nt,param.layer_bottom(rline)));
      trend(bins) = polyval(detrend_poly,x_axis(bins,rline));
      
      % Section 1: Constant
      trend(1:bins(1)-1) = trend(bins(1));
      
      % Section 3: Constant
      trend(bins(end):end) = trend(bins(end));
      
      if 0
        %% Debug plots
        figure(1); clf;
        plot(lp(data(:,rline)));
        hold on;
        plot(trend);
      end
      data(:,rline) = data(:,rline) ./ 10.^(trend/10);
    end
    
  case 'tonemap'
    tmp = [];
    tmp(:,:,1) = abs(data);
    tmp(:,:,2) = abs(data);
    tmp(:,:,3) = abs(data);
    % for generating synthetic HDR images
    data = tonemap(tmp, 'AdjustLightness', [0.1 1], 'AdjustSaturation', 1.5);
    data = single(data(:,:,2));
   
  otherwise
    error('Invalid param.method %s.', param.method);
end
