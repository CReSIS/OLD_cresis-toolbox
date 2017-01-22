function layer = threshold(img_slice,params)
% layer = threshold(img_slice,params)
%
% Track the ice-surface layer in a 3D image slice.
%
% Inputs:
% img_slice: single slice from 3D imagery (i.e. img_slice is a 2D image)
% params: struct array controlling how threshold works
%  .Noise_range_bin : "Noise range bin" specified by user to elimate the 
%                     clutter part above the surface and threshold function
%                     only allow surface detection to occur at a later bin
%                     than the noise range bin.
%  FOLLOWING FIELDS ARE NOT CURRENTLY SUPPORTED:
%  .row_range   :  a vector specifying relative row range to search in.
%                  E.g. "-30:30".
%  .threshold   :  surface choosen is above the noise floor by "threshold dB"
%                  and if surface contains no values above threshold, then
%                  the value is set to NaN.
%
% Outputs:
% layer: row corresponding to the surface for each column in the slice
%
% Author: Sravya Athinarapu
%
% See also: slicetool_threshold.m

Range_bins = 1 : size(img_slice,1);
slice_columns = size(img_slice,2);
layer = zeros(1,slice_columns);
% To remove the clutter above the surface
index_noise = params.Noise_Range_bin; 
Range_bins_noise_free = Range_bins(index_noise + 1:end);

for col = 1 : slice_columns
    rline_data_noise_free = img_slice(index_noise+1:end,col);
    
    if col == 1
        [~,id_req_pixel] = max(rline_data_noise_free);
        layer(col) = Range_bins_noise_free(id_req_pixel);
        prev_layer = layer(col) - index_noise;
    else
        
        if col <= slice_columns/2 + 3  % for left part and center part of the slice
            data_col = rline_data_noise_free (1 : prev_layer);
            Range_bins_col = Range_bins_noise_free (1 : prev_layer );
            [~,id_req_pixel_both] = max(data_col);
            layer(col) = Range_bins_col(id_req_pixel_both);
            prev_layer = layer(col) - index_noise;
            
        else    % for right part of slice
            layer_allow_range_bins = 30;
            data_col = rline_data_noise_free(prev_layer : prev_layer + layer_allow_range_bins);
            Range_bins_col = Range_bins_noise_free(prev_layer : prev_layer + layer_allow_range_bins);
            [~,id_req_pixel_both] = max(data_col);
            layer(col)  = Range_bins_col(id_req_pixel_both);
            prev_layer = layer(col) - index_noise;
            
        end
        
    end
    clear data_col Range_bins_col id_req_pixel_both
end

layer = layer';
return