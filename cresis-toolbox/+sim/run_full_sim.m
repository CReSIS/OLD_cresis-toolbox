% Script sim.run_full_sim
%
% Calls all necessary scipts and functions to demostrate the FullSimulator
%
% Authors: John Paden, Hara Madhav Talasila

if 1 % enable this to rerun everything
  
  sim.input_full_sim
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
