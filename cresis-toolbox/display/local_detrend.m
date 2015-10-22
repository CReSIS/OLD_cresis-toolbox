function data = local_detrend(data, B_noise, B_sig, cmd, minVal)
% data = local_detrend(data, B_noise, B_sig, cmd, minVal)
%
% data = 2D input data matrix (linear power)
% B_noise = 2x1 vector giving size of noise window
%   B_noise(1) = row/first-dim size
%   B_noise(2) = col/second-dim size
% B_sig = 2x1 vector giving size of signal window
%   B_sig(1) = row/first-dim size
%   B_sig(2) = col/second-dim size
% cmd = scalar command to combine estimates
%   1: combine with max
%   2: combine with min
%   3: use just method 1 [default]
%   4: use just method 2
%   5: use no detrend
% minVal = scalar minimum value of detrend [-inf default]
%   Units of dB (so 10*log10(data))
%
% This local detrend function combines to noise power
% estimates:
% 1. average power for each bin (mean on second dimension)
% 2. average power from noise window (mean from filter2)
% The signal power is then scaled relative to the noise power.
%
% imagesc(10*log10(local_detrend(data,[80 20], [4 2])));
%
% imagesc(10*log10(local_detrend(data,[80 20], [4 2], , 220)));
%
% Look into:
% Data = tonemap(repmat(Data,[1 1 3]), 'AdjustLightness', [0.1 1], 'AdjustSaturation', 1.5);
% Data = Data(:,:,2);

debugLevel = 1;

if ~exist('cmd','var') || isempty(cmd)
  cmd = 3;
end
if ~exist('minVal','var') || isempty(minVal)
  minVal = -inf;
end
if ~exist('B_noise','var') || isempty(B_noise)
  B_noise = [1 1];
end
if ~exist('B_sig','var') || isempty(B_sig)
  B_sig = [1 10];
end

if cmd == 5
  data = filter2(ones(B_sig(1),B_sig(2))/(B_sig(1)*B_sig(2)),data);
  return;
end

% Prevent log function from returning -inf
data(data==0) = min(data(data>0));

% Detrend method 1
detrend1 = mean(data,2);
filter_length = ceil(size(data,1)/25);
detrend1_tmp = 10*log10(filter(exp(-(linspace(-2,2,filter_length)).^2),1,detrend1));
detrend1(1:ceil(filter_length/2)) = detrend1_tmp(filter_length);
detrend1(floor(filter_length/2)+1:end-ceil(filter_length/2)) = detrend1_tmp(filter_length+1:end);
detrend1(end-ceil(filter_length/2)+1:end) = detrend1_tmp(end);

% [B,A] = butter(2,1/(0.05*size(data,1)));
% detrend1 = filtfilt(B,A,10*log10(detrend1));
if debugLevel == 2
  figure(1); clf;
  plot(detrend1);
end
detrend1 = repmat(detrend1,[1 size(data,2)]);

if cmd == 1 || cmd == 2 || cmd == 4
  % Detrend method 2
  detrend2 = filter2(ones(B_noise(1),B_noise(2))/(B_noise(1)*B_noise(2)),data);
  %detrend2 = medfilt2(data,[B_noise(1),B_noise(2)]);
  detrend2 = 10*log10(detrend2);
  if debugLevel == 2
    figure(1); clf;
    imagesc(detrend2);
    colorbar;
  end
end

% Detrend methods combined
switch cmd
  case 1
    detrendc = max(detrend1,detrend2);
  case 2
    detrendc = min(detrend1,detrend2);
  case 3
    detrendc = detrend1;
  case 4
    detrendc = detrend2;
  otherwise
    warning('Invalid command (using default of 1)');
    detrendc = max(detrend1,detrend2);
end

% The detrend normalizes the data, this limits the
% normalization (noise powers below the threshold
% won't be normalized more than the threshold)
detrendc(detrendc<minVal) = minVal;

% Convert back to linear power scale
detrendc = 10.^(detrendc/10);

% Filter to find signal power and apply detrend
data = filter2(ones(B_sig(1),B_sig(2))/(B_sig(1)*B_sig(2)),data)./detrendc;
%data = medfilt2(data,[B_sig(1),B_sig(2)])./detrendc;

return;
