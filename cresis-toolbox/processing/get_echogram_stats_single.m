function [layer_SNR] = get_echogram_stats_single(img,layer_row,noise_bin_rng)
% [layer_SNR] = get_echogram_stats_single(img,layer_row,noise_bin_rng,noise_mode)
%
% Returns signal to noise ratio for layer. Signal power is estimated by
% taking the value at img(layer_row(col),col). Noise power is found by
% taking the average of all the pixels that fall in
% img(layer_row(col)+noise_bin_rng(1):end-noise_bin_rng(end),col) and
% estimating the power of those pixels. If there are no pixels available
% for the noise power, then layer_SNR will return all NaN.
%
% INPUTS:
%
% img:
%  Assumed to be log power image of size Nt by Nx
%
% layer_row:
%  1 by Nx vector such that layer_row(col) indicates the row in img for the
%  layer for column "col".
%
% noise_bin_rng:
%  2 element vector that specifies the row range of the pixels to calculate
%  the noise power with. The first number is the number of pixels after the
%  layer to start at and the second number is the number of pixels to
%  ignore at the bottom of the image. For example, the bins that will be
%  used to calculated the noise for column, col, are:
%    noise_bins = layer_row(col) + noise_bin_rng(1) : Nt - noise_bin_rng(end)
%
% OUTPUTS:
%
% layer_SNR:
%  1 by Nx vector such that layer_SNR(col) indicates the signal to noise
%  ratio (SNR) in dB for the layer for column "col".
%
% EXAMPLE:
% 
% % Load data
% fn = '/N/dcwan/projects/cresis/output/rds/2018_Greenland_P3/CSARP_standard/20180405_01/Data_20180405_01_014.mat';
% mdata = load(fn);
% % Convert data to log scale
% img = 10*log10(mdata.Data);
% 
% % Load layer
% fn = '/N/dcwan/projects/cresis/output/rds/2018_Greenland_P3/CSARP_layerData/20180405_01/Data_20180405_01_014.mat';
% lay = load(fn);
% % Convert layer two way travel time to rows
% layer_row = interp1(mdata.Time,1:length(mdata.Time),lay.layerData{2}.value{2}.data);
% layer_row = round(interp1(lay.GPS_time,layer_row,mdata.GPS_time,'linear','extrap'));
% 
% % For noise, use bins starting 400 bins after the layer and stopping 50
% % bins from the end of each column.
% noise_bin_rng = [400 50];
% 
% % Calculate and plot SNR for the layer
% layer_SNR = get_echogram_stats_single(img,layer_row,noise_bin_rng);
% plot(layer_SNR); grid on; xlabel('Column'),ylabel('SNR (dB)');
%
% Author: John Paden

if 0
  % Example/Test
  
  % Load data
  fn = '/N/dcwan/projects/cresis/output/rds/2018_Greenland_P3/CSARP_standard/20180405_01/Data_20180405_01_014.mat';
  mdata = load(fn);
  % Convert data to log scale
  img = 10*log10(mdata.Data);
  
  % Load layer
  fn = '/N/dcwan/projects/cresis/output/rds/2018_Greenland_P3/CSARP_layerData/20180405_01/Data_20180405_01_014.mat';
  lay = load(fn);
  % Convert layer two way travel time to rows
  layer_row = interp1(mdata.Time,1:length(mdata.Time),lay.layerData{2}.value{2}.data);
  layer_row = round(interp1(lay.GPS_time,layer_row,mdata.GPS_time,'linear','extrap'));
  
  noise_bin_rng = [401 500];
  
  % Calculate and plot SNR for the layer
  layer_SNR = get_echogram_stats_single(img,layer_row,noise_bin_rng);
  plot(layer_SNR);
  
end

% Setup
Nt = size(img,1);
Nx = size(img,2);
layer_SNR = zeros(1,Nx);
signal_bin_rng = 0;
noise_bin_rng(end) = max(0,noise_bin_rng(end));
noise_mask = false(size(img)); % memory inefficient, but easy

% Check the signal and noise for each column (range line) of the img
for rline = 1:Nx
  signal_rows = max(1,layer_row(rline) + signal_bin_rng(1)) : min(Nt,layer_row(rline) + signal_bin_rng(end));
  
  % Put the signal into layer_SNR:
  if isempty(signal_rows)
    % No valid signal bins for this column
    layer_SNR(rline) = NaN;
  else
    layer_SNR(rline) = max(img(signal_rows,rline));
  end
  
  noise_rows = max(1,layer_row(rline) + noise_bin_rng(1)) : min(Nt,Nt - noise_bin_rng(end));
  
  noise_mask(noise_rows,rline) = true;
end

if any(noise_mask(:))
  % Divide signal in layer_SNR by average noise power
  layer_SNR = layer_SNR - 10*log10(mean(10.^(img(noise_mask)/10)));
else
  % No valid noise pixels found
  layer_SNR(:) = NaN;
end
