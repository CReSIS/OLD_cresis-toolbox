function RMSE = crosstrack_rmse(param,results,src_sep)
% RMSE = sim.crosstrack_rmse(param,results,src_sep)
%
% Compute root mean squared error (RMSE) for each source from
% sim.crosstrack param and results. Results are separated based on
% the source separation input argument.
%
% param and results: directly from crosstrack.m
% src_sep: array of doa separations, default is [0:180]
% RMSE: structure of results
%  .num_src(actual number of sources, estimated number of sources,
%   actual values of sources, estimated values of sources)
%  .RMSE_total: RMSE of all sources
%  .num_total: number of pixels contributing to RMSE_total
%  .RMSE{actual_num_sources}(src_separation_idx): RMSE error of sources
%    binned based on source separation (defined as the minimum distance to
%    the nearest other source) and binned based on the number of sources
%  .num{actual_num_sources}(src_separation_idx): number of pixels
%    contributing to each RMSE bin
%  .RMSE{actual_num_sources}(source_idx,:): results restructured so that
%    errors are grouped by numbers of sources and the error for each source
%
% Author: John Paden, Theresa Stumpf
%
% See also: crosstrack.m, crosstrack_example.m

physical_constants;

% Extract results from "crosstrack":

% Nt: number of range axis samples
% Nsig: number of sources (aka signals)
% Nx: number of along track samples

% range: range to target (Nt by 1 vector, units: seconds)
range = results.array_param.wfs.time(results.array_param.bins)*c/2;

% doa: direction of arrival to the target
% power: scattering power of target (higher is better)
% cost: cost function result (lower is better)
% hessian: hessian of cost function (higher is better)
%   (Nt by Nsig by Nx array, units: radians)
doa = results.tomo.doa;
power = results.tomo.power;
cost = results.tomo.cost;
hessian = results.tomo.hessian;

% surf_model: the target surface grid
surf_model = results.surf_model;
% param: parameters used to setup the simulation
param = results.param;
% array_param: array processing parameters
array_param = results.array_param;

% x: along track position of target
% y: cross track position of target
% z: elevation position of target
%   (Nt by Nsig by Nx array, units: meters)
y = bsxfun(@times,range,sin(doa));
z = bsxfun(@times,-range,cos(doa));
x = repmat(permute(surf_model.x(array_param.lines),[1 3 2]),[size(y,1) size(y,2) 1]);

% x_Axis: along track positions (Nx by 1, units: meters)
x_axis = param.monte.target_param{1}.x(array_param.lines);

%% Determine the number of sources and their DOA for each image pixel
Nt = length(results.array_param.bins);
Nx = length(results.array_param.lines);
actual_doa = cell(Nt,Nx);

dr = range(2) - range(1);
for x_idx = 1:Nx
  line = results.array_param.lines(x_idx);
  for y_idx = 1:length(surf_model.y)
    % Determine range to the target
    R = sqrt(surf_model.y(y_idx).^2 + surf_model.z(y_idx,line).^2);
    src_doa = atan2(surf_model.y(y_idx),surf_model.z(y_idx,line));
    rbins = find(abs(R-range) < dr/2*2);
    for rbin = rbins(:).'
      actual_doa{rbin,x_idx} = [actual_doa{rbin,x_idx} src_doa];
    end
  end
end

actual_doa_len = cellfun(@length,actual_doa);
max_doa_len = max(actual_doa_len(:));

figure(1); clf;
imagesc(actual_doa_len);
h_axis = gca;
h_colorbar = colorbar;
set(get(h_colorbar,'YLabel'),'String','Number of sources');
xlabel('Along-track range line');
ylabel('Range bin');
jet_colormap = jet(max_doa_len+1);
colormap([0 0 0; jet_colormap(1:end-1,:)]);

%% Group and filter DOAs that are close together

% Convert doa from theta to ky
k = 4*pi*param.src.fc/c;
actual_doa = cellfun(@(theta) k*sin(theta),actual_doa,'UniformOutput',false);

% Determine dky resolution to use
My = 4;
phase_center = param.src.lever_arm.fh(param.src.lever_arm.fh_args{:});
Y = max(phase_center(2,:)) - min(phase_center(2,:));
dky = 2*pi / Y / My;

for x_idx = 1:Nx
  for rbin = 1:Nt
    actual_doa{rbin,x_idx} = sort(actual_doa{rbin,x_idx});
    doa_idx = 1;
    while doa_idx <= length(actual_doa{rbin,x_idx})
      % Collect all doa's within the resolution cell and average into one
      % doa bin
      doa_bin_edge = actual_doa{rbin,x_idx}(doa_idx) + dky;
      doa_end_idx = doa_idx;
      while doa_end_idx+1 <= length(actual_doa{rbin,x_idx}) ...
          && actual_doa{rbin,x_idx}(doa_end_idx+1) < doa_bin_edge
        doa_end_idx = doa_end_idx + 1;
      end
      actual_doa{rbin,x_idx} ...
        = [actual_doa{rbin,x_idx}(1:doa_idx-1), ...
          mean(actual_doa{rbin,x_idx}(doa_idx:doa_end_idx)), ...
          actual_doa{rbin,x_idx}(doa_end_idx+1:end)];
        
      % Skip to the next doa index that was not binned
      doa_idx = doa_end_idx + 1;
    end
  end
end

% Convert doa from ky to theta
actual_doa = cellfun(@(ky) asin(ky/k),actual_doa,'UniformOutput',false);

actual_doa_len = cellfun(@length,actual_doa);
max_doa_len = max(actual_doa_len(:));

figure(2); clf;
imagesc(actual_doa_len);
h_axis(2) = gca;
h_colorbar = colorbar;
set(get(h_colorbar,'YLabel'),'String','Number of sources');
xlabel('Along-track range line');
ylabel('Range bin');
jet_colormap = jet(max_doa_len+1);
colormap([0 0 0; jet_colormap(1:end-1,:)]);

%% For each image pixel, compare estimated DOAs to actual DOAs
est_error = cell(Nt,Nx);

% RMSE = [];
% RMSE.RMSE_total = 0;
% RMSE.num_total = 0;
% RMSE.RMSE = cell(max_doa_len);
% RMSE.num = cell(max_doa_len);

for x_idx = 1:Nx
  for rbin = 1:Nt
    if ~isempty(actual_doa{rbin,x_idx})
      est_doa = squeeze(doa(rbin,:,x_idx));
      %       est_power = squeeze(power(rbin,:,x_idx))
      %       est_cost = cost(rbin,x_idx)
      %       est_hessian = squeeze(hessian(rbin,:,x_idx))
      
      % Match up doa's (assume best case matchup)
      for doa_idx = 1:length(actual_doa{rbin,x_idx})
        [~,best_idx] = min(abs(actual_doa{rbin,x_idx}(doa_idx) - est_doa));
        est_error{rbin,x_idx}(doa_idx) = actual_doa{rbin,x_idx}(doa_idx) - est_doa(best_idx);
      end
    end
  end
end

%% Compute the Root Mean Squared Error (RMSE)

if 0
  % Grab mean of all
  RMSE_mean = cellfun(@(bin_error) single(sqrt(mean(abs(bin_error).^2))),est_error);
  num = cellfun(@(bin_error) length(bin_error),est_error);
  
elseif 0
  % Grab smallest
  RMSE_mean = cellfun(@(bin_error) single(sum(min(abs(bin_error)))),est_error);
  num = cellfun(@(bin_error) ~isempty(bin_error),est_error);
  
elseif 1
  % Grab N smallest
  Nsig = array_param.Nsig;
  num = cellfun(@(bin_error) length(bin_error),est_error);
  est_error = cellfun(@(bin_error) sort(abs(bin_error)), est_error, 'UniformOutput', false);
  est_error = cellfun(@(bin_error) bin_error(1:min(Nsig,end)), est_error, 'UniformOutput', false);
  RMSE_mean = cellfun(@(bin_error) single(sqrt(mean(abs(bin_error).^2))),est_error);
  RMSE_mean(num>2) = NaN;
  num = cellfun(@(bin_error) length(bin_error),est_error);
  num(isnan(RMSE_mean)) = NaN;
end

% Remove outliers
threshold = 2 * nanstd(RMSE_mean(:));
RMSE_mean(RMSE_mean>threshold) = NaN;

num_total = nansum(num(:));
RMSE_total = nansum(RMSE_mean(:) .* num(:)) ./ nansum(num(:)) * 180/pi;

figure(3); clf;
colormap(jet(256));
imagesc(RMSE_mean * 180/pi);
h_axis(3) = gca;
h_colorbar = colorbar;
set(get(h_colorbar,'YLabel'),'String','RMS error (deg)');
caxis([0 5]);
xlabel('Along-track range line');
ylabel('Range bin');
linkaxes(h_axis,'xy');

RMSE.num_total = num_total;
RMSE.RMSE_total = RMSE_total;

return;
