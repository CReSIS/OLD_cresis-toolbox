function [data] = fir_dec(data, Bfilter, dec_factor, start_idx, Nidxs, renormalize_en, phase_weights)
% [data] = fir_dec(data, Bfilter, dec_factor, start_idx, Nidxs, renormalize_en, phase_weights)
%
% Applies a finite impulse response filter to the data and then decimates.
% The FIR only operates on the second dimension of 2D datasets.
%
% With 2 input arguments, it does simple coherent averaging (also known as
% stacking or presumming).
%
% data = 1D row or column vector or 2D array
%
% If only 2 arguments and second argument is a scalar:
%  Bfilter = scalar, amount to presum (boxcar window and then decimate)
%    boxcar window is normalized so that it is a mean (as opposed to sum)
%    If Bfilter (the amount to presum) is less than one, then nothing
%    happens.
%
% If 3 or more arguments:
%  Bfilter = FIR filter coefficients (row vector), must be even order
%    (i.e. odd length)
%  dec_factor = Decimation factor to apply to the data (default is 1)
%  start_idx = first output corresponds to this input index (default is 1)
%  Nidxs = number of outputs (default is the max supported by the input)
%  renormalize_en = renormalize edge effects, default is true
%
% Authors: John Paden
%
% See also: decimate, resample, downsample

if nargin == 2 && length(Bfilter) == 1
  coh_avgs = Bfilter;
  if coh_avgs >= 1
    siz = size(data);
    if siz(1) == 1
      new_len     = floor(siz(2)/coh_avgs);
      nearest_len = floor(siz(2)/coh_avgs)*coh_avgs;
      data        = mean(reshape(data(1:nearest_len),coh_avgs,new_len),1);
    elseif siz(2) == 1
      new_len     = floor(siz(1)/coh_avgs);
      nearest_len = floor(siz(1)/coh_avgs)*coh_avgs;
      data        = mean(reshape(data(1:nearest_len),coh_avgs,new_len),1).';
    else
      new_len     = floor(siz(2)/coh_avgs);
      nearest_len = floor(siz(2)/coh_avgs)*coh_avgs;
      data        = reshape(mean(reshape(data(:,1:nearest_len), ...
        [siz(1) coh_avgs new_len]),2),[siz(1) new_len]);
    end
  end
  
else
  if ~exist('start_idx','var') || isempty(start_idx)
    start_idx = 1;
  end
  if ~exist('dec_factor','var') || isempty(dec_factor)
    dec_factor = 1;
  end
  
  Nidxs_max = ceil((size(data,2) - start_idx + 1) / dec_factor);
  if ~exist('Nidxs','var') || isempty(Nidxs) || Nidxs > Nidxs_max
    Nidxs = Nidxs_max;
  end
  
  if ~exist('renormalize_en','var') || isempty(renormalize_en)
    renormalize_en = true;
  end
  
  if ~exist('phase_weights','var') || isempty(phase_weights)
    phase_weights = [];
  end
  
  if size(Bfilter,1) ~= 1
    error('Bfilter must be a row vector');
  end
  
  filter_order = (length(Bfilter) - 1);
  
  if mod(filter_order,2)
    error('Only even filter orders allowed (Bfilter should be odd length)');
  end
    
  if size(data,2) == 1
    data = data.';
  end
  
  out_data = zeros(size(data,1), Nidxs, class(data));
  
  for out_idx = 1:Nidxs
    renormalize = false;
    
    idx_rng(1) = start_idx + (out_idx-1)*dec_factor - filter_order/2;
    if idx_rng(1) < 1
      filter_rng(1) = 2 - idx_rng(1);
      idx_rng(1) = 1;
      renormalize = true;
    else
      filter_rng(1) = 1;
    end
    
    idx_rng(2) = start_idx + (out_idx-1)*dec_factor + filter_order/2;
    if idx_rng(2) > size(data,2)
      filter_rng(2) = size(Bfilter,2) - (idx_rng(2) - size(data,2));
      idx_rng(2) = size(data,2);
      renormalize = true;
    else
      filter_rng(2) = size(Bfilter,2);
    end
    
    if ~renormalize
      if isempty(phase_weights)
        out_data(:,out_idx) = sum(bsxfun(@times, data(:,idx_rng(1):idx_rng(2)), Bfilter),2);
      else
        dphase = phase_weights(idx_rng(1):idx_rng(2));
        dphase = conj(dphase./exp(1i * angle(mean(dphase))));
        out_data(:,out_idx) = sum(bsxfun(@times, data(:,idx_rng(1):idx_rng(2)), Bfilter .* dphase),2);
      end

    else
      if isempty(phase_weights)
        out_data(:,out_idx) = sum(bsxfun(@times, data(:,idx_rng(1):idx_rng(2)), ...
          Bfilter(:,filter_rng(1):filter_rng(2))), 2);
      else
        dphase = phase_weights(idx_rng(1):idx_rng(2));
        dphase = conj(dphase./exp(1i * angle(mean(dphase))));
        out_data(:,out_idx) = sum(bsxfun(@times, data(:,idx_rng(1):idx_rng(2)), ...
          Bfilter(:,filter_rng(1):filter_rng(2)) .* dphase), 2);
      end
      
      if renormalize_en
        renormalization_factor ...
          = sum(Bfilter(1,:)) / sum(Bfilter(1,filter_rng(1):filter_rng(2)));
        out_data(:,out_idx) = out_data(:,out_idx) * renormalization_factor;
      end
    end
  end

  data = out_data;
end

return;

