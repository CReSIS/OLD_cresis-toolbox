function data = echo_norm(mdata, param)
% data = echo_norm(mdata, param)
%
% The trend of the data is estimated using various methods and this trend
% is removed from the data.
%
% INPUTS:
%
% data = 2D input data matrix (log power)
%
% param: struct controlling how the normalization is done
%
% OUTPUTS:
%
% data: normalized input (log power normalized)
%
% Examples:
%
% fn = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_standard/20140512_01/Data_20140512_01_018.mat';
% mdata = load(fn); mdata.Data = 10*log10(mdata.Data);
%
% imagesc(echo_norm(mdata)); colorbar; caxis([0 1]);
%
% imagesc(echo_norm(mdata,struct('scale',[25 255]))); colorbar; caxis([0 1]);
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

% valid_max_range_dB: 2 element numeric vector; specifies the valid range
% for the max value; default is [-inf inf] which effectively disables
% this valid max value range constraint;
if ~isfield(param,'valid_max_range_dB') || isempty(param.valid_max_range_dB)
  param.valid_max_range_dB = [30 inf];
end

if ~isfield(param,'scale') || isempty(param.scale)
  param.scale = [0.1 0.9];
end


% Estimate the noise
if isstruct(mdata)
  mdata.Data = 10.^(mdata.Data/10);
  noise = 10*log10(mean(echo_noise(mdata, param)));
else
  noise = 10*log10(mean(echo_noise(10.^(mdata/10), param)));
end

% Determine max of data
max_scalar = min(noise+param.valid_max_range_dB(2), ...
  max(noise+param.valid_max_range_dB(1), ...
  max(data(:))));

% Scale and offset data
data = param.scale(1) + (data-noise) / (max_scalar - noise) * param.scale(2);
