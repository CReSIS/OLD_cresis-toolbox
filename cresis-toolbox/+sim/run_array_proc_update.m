% script run_array_proc_update
%
% Example of how to create simulated data using crosstrack_example and then
% see how the cost function changes when that result is perturbed.
% 
% Author: John Paden
%
% See also: array_proc_update.m, sim.crosstrack_example.m

% =========================================================================
%% Step 1. Run crosstrack_example.m example 1 with debug_level set to 2 or
% higher
% =========================================================================

fprintf('Ensure sim.crosstrack_example.m runs example #1!!!\n\n');
sim.crosstrack_example;

% Extract results from "crosstrack_example":

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
surf_model = results.surf_model;
param = results.param;
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

if 0
  % Debug: Basic surface extraction/gridding process
  figure(2); clf;
  imagesc(surf_model.x, surf_model.y, surf_model.z-param.monte.target_param{1}.z.mean);
  title('Ground truth');
  h_axis = gca;
  xlabel('Along-track (m)');
  ylabel('Cross-track (m)');
  cc = caxis;
  h_cb = colorbar;
  set(get(h_cb,'YLabel'),'String','WGS-84 elevation (m)');
  
  good_mask = zeros(size(tomo.doa));
  good_mask = good_mask | db(tomo.power) > 10;
  good_mask = good_mask | repmat(permute(tomo.cost,[1 3 2]),[1 size(good_mask,2) 1]) < -25;
  
  z_grid = griddata(double(x(good_mask)),double(y(good_mask)),double(z(good_mask)), ...
    surf_model.x,surf_model.y);
  
  figure(3); clf;
  imagesc(surf_model.x, surf_model.y, z_grid-param.monte.target_param{1}.z.mean);
  h_axis(2) = gca;
  title('Array processing with basic surface extraction');
  xlabel('Along-track (m)');
  ylabel('Cross-track (m)');
  caxis(cc);
  h_cb = colorbar;
  set(get(h_cb,'YLabel'),'String','WGS-84 elevation (m)');
end

% =========================================================================
%% Step 2. Perturb the output
% =========================================================================

% Step 2a: Choose which bin_idx and line_idx you want to update:
bin_idx = 210;
line_idx = 1;
% Set the new DOA for this pixel (perturb first source):
doa_new = doa(bin_idx,:,line_idx);
doa_new(1) = doa_new(1) + 0.0;

% Step 2b: Update the results with the new DOA:
[tomo_update] = array_proc_update(results.sim_data,results.array_param,bin_idx,line_idx,doa_new);

% Step 2c: Compare new results with old results:
fprintf('\n');
fprintf('Direction of arrivals (original and perturbed)\n');
doa(bin_idx,:,line_idx)
tomo_update.doa

fprintf('Source power (original and perturbed)\n');
power(bin_idx,:,line_idx)
tomo_update.power

fprintf('Local cost function (original and perturbed)\n');
cost(bin_idx,line_idx)
tomo_update.cost

return
