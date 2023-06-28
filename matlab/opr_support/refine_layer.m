function layer_bins = refine_layer(data, layer_bins, param)
% layer_bins = refine_layer(data, layer_bins, param)
%
% Refines layer_bins by interpolating data around the bins indicated in
% layer_bins.
%
% data: Nt by Nx matrix
% layer_bins: 1 by Nx vector (typically contains indices from 1 to Nt of the location of a layer for each column in data)
% param
%   .over_sample: amount to over sample the data (positive integer)
%   .method: string containing the interpolationg method to use (e.g. 'spline')
%   .bin_rng: a vector representing the range of bins to use in the
%      interpolation around the bin indicated by layer_bins, usually -N:N
%
% layer_bins: 1 by Nx vector output with updated values
%
% A = randn(100);
% A = filter2(ones(10),A);
% [tmp max_bins] = max(A);
% max_bins = max_bins + round(2*randn(size(max_bins)));
% max_bins_refine = refine_layer(A, max_bins, struct('over_sample',10,'method','spline','bin_rng',-5:5));
% figure;
% imagesc(A); colormap(1-gray(256));
% hold on;
% plot(max_bins);
% plot(max_bins_refine,'r');
% hold off;
%
% Author: John Paden

for rline = 1:size(data,2)
  bin_start = floor(layer_bins(rline) + min(param.bin_rng));
  if bin_start < 1
    bin_start = 1;
  end
  if bin_start > size(data,1)-1
    bin_start = size(data,1)-1;
  end
  bin_stop = ceil(layer_bins(rline) + max(param.bin_rng));
  if bin_stop < 2
    bin_stop = 2;
  end
  if bin_stop > size(data,1)
    bin_stop = size(data,1);
  end
  
  new_bin_axis = bin_start:1/param.over_sample:bin_stop;
  
  data_over_sampled = interp1(bin_start:bin_stop, data(bin_start:bin_stop,rline), new_bin_axis, param.method);
  
  [tmp idx] = max(data_over_sampled);
  
  layer_bins(rline) = bin_start + (idx - 1) / param.over_sample;
end

return
