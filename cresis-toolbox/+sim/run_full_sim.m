% Script sim.run_full_sim
%
% Calls all necessary scipts and functions to demostrate the FullSimulator
%
% Authors: John Paden, Hara Madhav Talasila
hara;

% =========================================================================
fprintf('=====================================================================\n');
fprintf('%s: (%s)\n', mfilename, datestr(now));
fprintf('=====================================================================\n');


% param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140410/param.mat';
  param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20140502/param.mat';
%   param_fn = '/cresis/snfs1/dataproducts/ct_data/ct_tmp/sim3D/rds/2014_Greenland_P3sim/20120330/param.mat';
  
if 1 % enable this to rerun everything
  
  param_fn = sim.input_full_sim;
  
  run_load_data(7, param_fn); % example 7
  
  sim.run_qlook(1, param_fn);
  sim.run_qlook(0, param_fn);
  return;
  sim.run_sar(1, param_fn);
  sim.run_sar(0, param_fn);
  
  sim.run_array(1, param_fn);
  sim.run_array(0, param_fn);
  
elseif 0 % OLD -- enable this to rerun everything
  
  sim.input_full_sim;
  run_load_data; % example 7
  
  sim.run_qlook(1);
  sim.run_qlook(0);
  
  sim.run_sar(1);
  sim.run_sar(0);
  
  sim.run_array(1);
  sim.run_array(0);
  
else % just for plots
  
  sim.run_qlook(0);
  
  sim.run_sar(0);
  
  sim.run_array(0);
  
end
