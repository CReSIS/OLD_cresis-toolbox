function cluster_get_chain_list(param)
% cluster_get_chain_list(param)
%
% Gets a list of chains using the specified data location.
%
% param: parameter structure, leave undefined for default (default uses
%  global gRadar variable)

if nargin == 0
  global gRadar;
  param = gRadar;
end

chain_fns = get_filenames(param.cluster.data_location,'chain_','','',struct('type','f'));

if isempty(chain_fns)
  fprintf('No chains found in %s\n', param.cluster.data_location);
end

chain_ids = zeros(size(chain_fns));
time_stamps = zeros(size(chain_fns));
for chain_idx = 1:length(chain_fns)
  [tmp chain_fn_name] = fileparts(chain_fns{chain_idx});
  chain_id = strtok(chain_fn_name(7:end),'_');
  chain_id = str2double(chain_id);
  chain_ids(chain_idx) = chain_id;
  finfo = dir(chain_fns{chain_idx});
  time_stamps(chain_idx) = finfo.datenum;
end

[~,chain_sorted_idxs] = sort(time_stamps);

for chain_idx = chain_sorted_idxs(:).'
  fprintf('%d: %s\n', chain_ids(chain_idx), chain_fns{chain_idx});
  try
    load(chain_fns{chain_idx});
  catch ME
    warning(ME.getReport());
    continue;
  end
  fprintf('  %s: %d chain(s) with batches: ', datestr(time_stamps(chain_idx)), length(ctrl_chain));
  first_batch = true;
  notes_str = '';
  for ctrl_chain_idx = 1:length(ctrl_chain)
    for batch_idx = 1:length(ctrl_chain{ctrl_chain_idx})
      if first_batch
        first_batch = false;
        fprintf('%d', ctrl_chain{ctrl_chain_idx}{batch_idx}.batch_id);
        if numel(ctrl_chain{ctrl_chain_idx}{batch_idx}.notes) >= 1
          notes_str = ctrl_chain{ctrl_chain_idx}{batch_idx}.notes{1};
        end
      else
        fprintf(', %d', ctrl_chain{ctrl_chain_idx}{batch_idx}.batch_id);
      end
    end
  end
  fprintf('\n');
  fprintf('  %s\n', notes_str);
  
end


