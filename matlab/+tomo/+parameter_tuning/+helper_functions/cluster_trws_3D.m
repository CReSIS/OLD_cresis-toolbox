function [ result ] = cluster_trws_3D( cpu_time, memory_req, detect_params_array, ...
  stats_array, num_combinations, sources, references )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

  cluster_compile;

  ctrl = cluster_new_batch;
  sparam.task_function = 'cluster_kernel_trws_3D';
  sparam.num_args_out = 1;
  sparam.notes = '';
  sparam.success = '';
  sparam.cpu_time = cpu_time;
  sparam.mem = memory_req;

  for idx = 1:num_combinations  
    dparam.argsin{1} = sources;
    dparam.argsin{2} = references;
    dparam.argsin{3} = detect_params_array(idx);
    dparam.argsin{4} = stats_array(idx); 

    [ctrl, ~] = cluster_new_task(ctrl,sparam,dparam,'dparam_save',0);
  end
  
  ctrl = cluster_save_dparam(ctrl);
  ctrl_chain = {{ctrl}};
  [~,chain_id] = cluster_save_chain(ctrl_chain);

  while true
      % These are the three lines of code that should be run to poll the job:
      fprintf('%s\nPoll the chain:\n',repmat('=',[1 40]));
      ctrl_chain = cluster_run(chain_id,false);
      cluster_save_chain(ctrl_chain,chain_id,false);

      % Check to see if any tasks are not complete and did not fail
      if ~any(isfinite(cluster_chain_stage(ctrl_chain)))
        % Done with all tasks (either completed or failed)
        break;
      end
      pause(30);
  end
  
  result = {};

  for idx = 1:num_combinations
    [~, result{idx}] = cluster_print(ctrl_chain{1}{1}.batch_id,idx,0);
  end

  cluster_cleanup(ctrl_chain);
  
end

