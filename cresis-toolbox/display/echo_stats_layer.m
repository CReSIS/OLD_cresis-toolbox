function [max_power,mean_power,neighbor_power,waveforms] = echo_stats_layer(data,layer,noise_bin_rng)
% [max_power,mean_power,neighbor_power,waveforms] = echo_stats_layer(data,layer,noise_bin_rng)
%
% Returns signal to noise ratio for layer. Signal power is estimated by
% taking the value at data(layer(col),col). Noise power is found by
% taking the average of all the pixels that fall in
% data(layer(col)+noise_bin_rng(1):end-noise_bin_rng(end),col) and
% estimating the power of those pixels. If there are no pixels available
% for the noise power, then layer_SNR will return all NaN.
%
% INPUTS:
%
% data:
%  Assumed to be log power image of size Nt by Nx
%
% layer:
%  1 by Nx vector such that layer(col) indicates the row in data for the
%  layer for column "col".
%
% noise_bin_rng:
%  2 element vector that specifies the row range of the pixels to calculate
%  the noise power with. The first number is the number of pixels after the
%  layer to start at and the second number is the number of pixels to
%  ignore at the bottom of the image. For example, the bins that will be
%  used to calculated the noise for column, col, are:
%    noise_bins = layer(col) + noise_bin_rng(1) : Nt - noise_bin_rng(end)
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
% data = 10*log10(mdata.Data);
% 
% % Load layer
% fn = '/N/dcwan/projects/cresis/output/rds/2018_Greenland_P3/CSARP_layerData/20180405_01/Data_20180405_01_014.mat';
% lay = load(fn);
% % Convert layer two way travel time to rows
% layer = interp1(mdata.Time,1:length(mdata.Time),lay.layerData{2}.value{2}.data);
% layer = round(interp1(lay.GPS_time,layer,mdata.GPS_time,'linear','extrap'));
% 
% % For noise, use bins starting 400 bins after the layer and stopping 50
% % bins from the end of each column.
% noise_bin_rng = [400 50];
% 
% % Calculate and plot SNR for the layer
% layer_SNR = get_echogram_stats_single(data,layer,noise_bin_rng);
% plot(layer_SNR); grid on; xlabel('Column'),ylabel('SNR (dB)');
%
% Author: John Paden

if 0
  % Example/Test
  
  % Load data
  fn = '/N/dcwan/projects/cresis/output/rds/2018_Greenland_P3/CSARP_standard/20180405_01/Data_20180405_01_014.mat';
  mdata = load(fn);
  % Convert data to log scale
  data = 10*log10(mdata.Data);
  
  % Load layer
  fn = '/N/dcwan/projects/cresis/output/rds/2018_Greenland_P3/CSARP_layerData/20180405_01/Data_20180405_01_014.mat';
  lay = load(fn);
  % Convert layer two way travel time to rows
  layer = interp1(mdata.Time,1:length(mdata.Time),lay.layerData{2}.value{2}.data);
  layer = round(interp1(lay.GPS_time,layer,mdata.GPS_time,'linear','extrap'));
  
  % For noise, use bins starting 400 bins after the layer and stopping 50
  % bins from the end of each column.
  noise_bin_rng = [400 50];
  
  % Calculate and plot SNR for the layer
  layer_SNR = get_echogram_stats_single(data,layer,noise_bin_rng);
  plot(layer_SNR);
  
end

%% Input Check
if ~exist('noise_bin_rng','var') || isempty(noise_bin_rng)
  % Default is to use all pixels below the layer
  noise_bin_rng = [0 0];
end

%% Setup
Nt = size(data,1);
Nx = size(data,2);
layer_SNR = zeros(1,Nx);
signal_bin_rng = 0;
noise_bin_rng(end) = max(0,noise_bin_rng(end));
noise_mask = false(size(data)); % memory inefficient, but easy

%% Check the signal and noise for each column (range line) of the data
for rline = 1:Nx
  signal_rows = max(1,layer(rline) + signal_bin_rng(1)) : min(Nt,layer(rline) + signal_bin_rng(end));
  
  % Put the signal into layer_SNR:
  if isempty(signal_rows)
    % No valid signal bins for this column
    layer_SNR(rline) = NaN;
  else
    layer_SNR(rline) = max(data(signal_rows,rline));
  end
  
  noise_rows = max(1,layer(rline) + noise_bin_rng(1)) : min(Nt,Nt - noise_bin_rng(end));
  
  noise_mask(noise_rows,rline) = true;
end

%% Final calculation of SNR
if any(noise_mask(:))
  % Divide signal in layer_SNR by average noise power
  layer_SNR = layer_SNR - 10*log10(mean(10.^(data(noise_mask)/10)));
else
  % No valid noise pixels found
  layer_SNR(:) = NaN;
end
