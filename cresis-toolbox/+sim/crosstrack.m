function results = crosstrack(param)
% results = sim.crosstrack(param)
%
% Function for simulating test data and testing direction of arrival
% algorithms on the test data.
%
% param: structure describing the simulation and testing
%  .src: structure describing signal source characteristics
%    .f0: low frequency (Hz)
%    .f1: high frequency (Hz)
%    .lever_arm.fh: function handle to a lever arm function that returns
%      the locations of the sensor phase centers. Phase centers returned
%      should be in a 3 by Nc array. First row is "x" or along-track
%      and is generally not used, second row is "y" , and third row is "z".
%      Nc is the number of sensors/channels of data.
%    .lever_arm.args: arguments to be passed to the lever arm function
%      handle
%    FIELDS REQUIRED BY sim.crossover_data.m
%  .monte: structure describing the monte carlo simulation setup
%    .runs: positive integer scalar representing number of simulation runs
%    .random_seed_offset: scalar double usually set to zero to allow
%                         repeatability
%
% results: structure containing results from simulation runs
%
% Author: Sean Holloway, John Paden, Theresa Stumpf
%
% See also: sim.crosstrack.m, sim.crosstrack_data.m, sim.crosstrack_example.m

%% Setup

% Load standard physical constants like c = speed of light
physical_constants;

% fc: center frequency
param.src.fc = (param.src.f0 + param.src.f1)/2;
% fs: sampling frequency
param.src.fs = abs(param.src.f1 - param.src.f0);

% phase_center: get phase center positions of antenna array
[phase_center] = param.src.lever_arm.fh(param.src.lever_arm.args{:});

% Convert phase_center's from aircraft coordinate system into flight coordinate system
% y_pc: cross track with positive pointing to the left
% z_pc: elevation with positive pointing up
param.src.y_pc   = -phase_center(2,:).';
param.src.z_pc   = -phase_center(3,:).';

%% Monte Carlo Trials
%==========================================================================

start_time = tic;
for run_idx = 1:param.monte.runs
  fprintf('Run %4d of %d (%.1f sec)\n', ...
    run_idx, param.monte.runs, toc(start_time));
  
  %% Setup random number generator
  rng_args(run_idx) = param.monte.random_seed_offset + run_idx;
  rng(rng_args(run_idx));
  
  %% Create layered surface model
  fprintf('  Create simulated data\n');
  surf_model.x = [];
  surf_model.y = [];
  surf_model.z = [];
  surf_model.rcs = [];
  surf_model.layer = [];
  if param.debug_level >= 3
    % Plot surface model
    figure(11); clf;
  end
  for surf_idx = 1:length(param.monte.target_param)
    new_surf = param.monte.target_func(param.monte.target_param{surf_idx});
    surf_model.x = [surf_model.x; new_surf.x];
    surf_model.y = [surf_model.y; new_surf.y];
    if surf_idx == 1
      surf_model.z = [surf_model.z; new_surf.z];
    else
      % Do not allow negative layer thickness
      new_surf.z(new_surf.z<0) = 0;
      surf_model.z = [surf_model.z; surf_model.z + new_surf.z];
    end
    surf_model.rcs = [surf_model.rcs; new_surf.rcs];
    surf_model.layer = [surf_model.layer; surf_idx*ones(size(new_surf.rcs(:,1)))];
    
    if param.debug_level >= 3
      % Plot surface model
      plot_style = {'b-','r:'};
      figure(11);
      plot(new_surf.y, new_surf.z(:,1), plot_style{1+mod(surf_idx-1,length(plot_style))});
      hold on;
    end
  end
  
  %% Create simulated data
  sim_data = sim.crosstrack_data(param, surf_model);
  
  if 0
    imagesc(squeeze(lp(sim_data(:,1,:))))
  end
  
  % =======================================================================
  %% Process simulated data
  % =======================================================================

  %% Setup array processing input arguments
  sim_data = {permute(sim_data,[1 3 4 5 2])};
  array_param = param.array_param;
  array_param.three_dim.en = true;
  for wf_adc_idx = 1:length(param.src.y_pc)
    array_param.fcs{1}{wf_adc_idx}.pos(2,:) = repmat(param.src.y_pc(wf_adc_idx), [1 size(sim_data{1},2)]);
    array_param.fcs{1}{wf_adc_idx}.pos(3,:) = repmat(param.src.z_pc(wf_adc_idx), [1 size(sim_data{1},2)]);
  end
  
  array_param.wfs.time = (param.src.t0:1/param.src.fs:param.src.t1).';
  array_param.wfs.fs = param.src.f1 - param.src.f0;
  array_param.wfs.fc = (param.src.f1 + param.src.f0)/2;
  array_param.imgs = {1};
  
  %% Loop through and run each array processing method
  for method_idx = 1:length(param.method.list)

    fprintf('  Array processing method %d\n', param.method.list(method_idx));
    
    %% Run array processing
    % Some notes about the outputs
    %  array_param.theta: radians, zero points toward -z_pc, increases toward positive y_pc
    array_param.method = param.method.list(method_idx);
    [array_param,tomo] = array_proc(array_param, sim_data);
    
    %% Debug Plots
    if param.debug_level >= 2
      range = array_param.wfs.time(array_param.bins)*c/2;
      
      if array_param.method < 7
        z = bsxfun(@times,-range,cos(theta));
        y = bsxfun(@times,range,sin(theta));
        y_axis = (min(y(:)):max(y(:))).';
        z_axis = min(z(:)):0.05:max(z(:));
        
        tomo.img_rect = griddata(y,z,double(tomo.img(:,:,1)),y_axis,z_axis);
        
        figure(1); clf;
        subplot(2,1,1);
        plot_style = {'b-','r:'};
        for layer = unique(surf_model.layer)
          plot(surf_model.y(surf_model.layer == layer), ...
            surf_model.z(surf_model.layer == layer,1)-param.monte.target_param{1}.z.mean, ...
            plot_style{1+mod(layer-1,length(plot_style))});
          hold on;
        end
        hold off;
        ylabel('Elevation (m)');
        xlim([-40 40]);
        ylim([-1 2.25]);
        h_axis = gca;
        grid on;
        
        subplot(2,1,2);
        imagesc(y_axis,z_axis-param.monte.target_param{1}.z.mean,lp(tomo.img_rect));
        set(gca,'YDir','normal');
        xlabel('Cross-track (m)');
        ylabel('Elevation (m)');
        xlim([-40 40]);
        ylim([-1 2.25]);
        h_axis(2) = gca;
        
        linkaxes(h_axis,'xy');
        
      else
        z = bsxfun(@times,-range,cos(tomo.doa(:,:,1)));
        y = bsxfun(@times,range,sin(tomo.doa(:,:,1)));
        
        figure(1); clf;
        subplot(2,1,1);
        plot_style = {'b-','r:'};
        for layer = unique(surf_model.layer).'
          plot(surf_model.y(surf_model.layer == layer), ...
            surf_model.z(surf_model.layer == layer,1)-param.monte.target_param{1}.z.mean, ...
            plot_style{1+mod(layer-1,length(plot_style))});
          hold on;
        end
        hold off;
        ylabel('Elevation (m)');
        xlim([-40 40]);
        ylim([-1 2.25]);
        h_axis = gca;
        grid on;
        
        subplot(2,1,2);
        plot(y,z-param.monte.target_param{1}.z.mean,'.');
        set(gca,'YDir','normal');
        xlabel('Cross-track (m)');
        ylabel('Elevation (m)');
        xlim([-40 40]);
        ylim([-1 2.25]);
        h_axis(2) = gca;
        grid on;
        linkaxes(h_axis,'xy');
        
        figure(2); clf;
        imagesc(surf_model.x, surf_model.y, surf_model.z-param.monte.target_param{1}.z.mean);
        title('Ground truth');
        h_axis = gca;
        xlabel('Along-track (m)');
        ylabel('Cross-track (m)');
        cc = caxis;
        h_cb = colorbar;
        set(get(h_cb,'YLabel'),'String','WGS-84 elevation (m)');
        
        z = bsxfun(@times,-range,cos(tomo.doa));
        y = bsxfun(@times,range,sin(tomo.doa));
        x = repmat(permute(surf_model.x(array_param.lines),[1 3 2]),[size(y,1) size(y,2) 1]);
        
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
        
        linkaxes(h_axis,'xy');
        
      end
    end
    
    if param.debug_level >= 3
      % End early for debug testing of outputs
      results.param = param;
      results.tomo = tomo;
      results.array_param = array_param;
      results.sim_data = sim_data;
      results.surf_model = surf_model;
      return
    end
    
    %% Determine number of sources errors and source DOA errors
    
    % To do...
    
  end
end

% Copy outputs into output argument structure

return;
