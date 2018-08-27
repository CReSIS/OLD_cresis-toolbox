function [ detect_data ] = set_detect_data( img_data, det_params, idx )
%SET_DETECT_DATA Summary of this function goes here
%   Detailed explanation goes here

  detect_data = img_data.img(:,:, idx);
  detect_data(detect_data > det_params.threshold) = det_params.threshold;
  detect_data = fir_dec(detect_data.',hanning(3).'/3,1).'; % what's this? just curious.
  
end

