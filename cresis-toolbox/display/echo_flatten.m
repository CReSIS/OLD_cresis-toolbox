function [mdata] = echo_flatten(mdata,layer_params,fill_value)
% [mdata] = echo_flatten(mdata,layer_params)
%
% Shifts columns up/down in echogram in order to flatten one or more 1D
% layers in the data. Defaults to flattening the surface stored in
% CSARP_layer if no layer_params given. If a single layer is given, all
% bins are shifted up and down the same amount. If more than one layer is
% given, then bins are adjusted according to a piecewise linear polynomial
% between each of the layers and nearest neighbor for extrapolating above
% and below the top and bottom layers respectively. If a slope field is
% given, then the bin shifts are made from left to right. Interpolation is
% used for the shifts, so oversampling (e.g. 2x) in the fast-time domain is
% recommended for best results when the shifts are non-integer. To
% compensate for elevation, just pass in a custom layer_params similar to
% opsCopyLayers and set the twtt = elev*c/2.
%
% INPUTS:
%
% mdata: echogram structure from an echogram file (e.g. returned from
% echo_load)
%
% layer_params: struct array controlling how the flattening is done if a
% character, then assumes it is a slope field file path (ct_filename_out
% will be used to determine filepath)
% Should have a 'slope' field containing layer_bins to flatten.
% Mirror of layer is flattened when 'mirror' field is present and true (as for elevation).
% When double, assumed to be layer_bins.
%
% fill_value: value with which to fill extraneous pixels that are 
% introduced above and below the shifted matrix. If a string, passed to
% interp_finite as the interp_method.
%
% interp_method: string containing 'linear' or 'interpft' (default is
% 'linear'). interpft is not recommended unless the data are Nyquist
% sampled (usually this means 2x interpolation in the fast-time before
% power detection occurs).
% 
% OUTPUTS:
%
% mdata: echogram structure (t0 and dt fields added)
%
% Examples:
%
% param = read_param_xls(ct_filename_param('rds_param_2019_Antarctica_Ground.xls'),'20200107_01');
% mdata = echo_load(param,'standard',1);
%
% Author: John Paden, Reece Mathews
%
% See also: echo_detrend, echo_filt, echo_mult_suppress, echo_noise,
% echo_norm, echo_param, echo_stats, echo_stats_layer, echo_xcorr,
% echo_xcorr_profile

physical_constants;

if isstruct(mdata)
  data = mdata.Data;
else
  data = mdata;
  mdata = struct();
  mdata.Data = data;
end

Nx = size(data, 2);
mirror = false;

if ~exist('layer_params','var') || isempty(layer_params)
    if ~isstruct(mdata) || ~isfield(mdata,'Surface') || isempty(mdata.Surface)
        error('If no layer params are passed, mdata must be struct with Surface field so that Surface can be flattened by default.');
    end

    surf_bins = interp1(mdata.Time, 1:length(mdata.Time), mdata.Surface);
    slope = surf_bins;
    
elseif isa(layer_params, 'double')
    % Passed layer_bins directly
    slope = layer_params;
    
elseif isstruct(layer_params)
    if isfield(layer_params, 'mirror') && layer_params.mirror
       mirror = true;
    end
    if ~isfield(layer_params, 'slope') || isempty(layer_params.slope)
        error('Slope field must be present on layer_params struct.');
    end
    slope = layer_params.slope;
    
    % TODO[reece]: When multiple layers are given, shift bins to flatten
    %              all layers (nearest neighbor above and below top and
    %              bottom)
    % TODO[reece]: Allow passing slope field files
end

slope_range = slope - slope(1);  % Use first column as baseline
if mirror
    slope_range = -slope_range;
end
new_size = size(data, 1) + floor(range(slope_range));
new_range = (1:new_size);
mdata.Data = nan(new_size, Nx);

cushioned_data = [nan(abs(ceil(min(slope_range))), Nx); data; nan(ceil(max(slope_range)), Nx)];

for c = 1:Nx
  mdata.Data(:, c) = interp1(1:size(cushioned_data, 1), cushioned_data(:, c), new_range + slope_range(c) - 1);
  new_idxs = find(~isnan(mdata.Data(:, c)));
  mdata.Vertical_Bounds([1 2], c) = [min(new_idxs), max(new_idxs)];
end

% TODO[reece]: Why does the nearest neighbor seem to sometimes come from
%              the other end of the previous column?
% TODO[reece]: What is the interp_method arg and interpft
% TODO[reece]: Update echogram's Time to reflect shifts

if ~exist('fill_value','var') || isempty(fill_value)
   fill_value = 'nearest'; 
end

if isa(fill_value, 'char') || isa(fill_value, 'string')
  mdata.Data = interp_finite(mdata.Data, nan, fill_value);
else
  mdata.Data(isnan(mdata.Data)) = fill_value;
end