% Scripts nsidc_delivery_script
%
% Delivery L1B and L2 data to nsidc
% For L1B
%        Input: .mat files
%        Output: .premet, .spatial and netcdf files
% For L2
%        INput: .csv file
%        Output: .premet, .spatial and ascii files
%
% Example:
%  See run_nsidc_delivery_script.m to run.
%
% Author: Yi Zhu, John Paden
%
% See also: type "nsidc_help.m"

%% Automated Section
fprintf('===============================================\n');
fprintf('NSIDC Delivery Script\n\n');

skip_phrase = 'do not process';

% Control parameters for netcdf file
netcdf_param = [];
netcdf_param(end+1).mat_name = 'GPS_time';
netcdf_param(end).cdf_name = 'time';
netcdf_param(end).dim_name = {'time'};
% netcdf_param(end).eval = @(x) (epoch_to_sod(x - utc_leap_seconds(x(1)),param.day_seg(1:8)));
% This is done later
netcdf_param(end).attributes = {'units' 'TO BE FILLED IN LATER', ...
  'calendar', 'standard', ...
  'long_name', 'Time of day UTC', ...
  'standard_name', 'time', ...
  'axis', 'T'};

netcdf_param(end+1).mat_name = 'Time';
netcdf_param(end).cdf_name = 'fasttime';
netcdf_param(end).dim_name = {'fasttime_dim'};
netcdf_param(end).eval = @(x) (x*1e6);
netcdf_param(end).attributes = {'units' 'microseconds', ...
  'long_name', '2-way travel time in useconds', ...
  'standard_name', 'time', ...
  'positive', 'down', ...
  'axis', 'Z'};

netcdf_param(end+1).mat_name = 'Elevation';
netcdf_param(end).cdf_name = 'altitude';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'standard_name','height', ...
  'units', 'meters', ...
  'positive', 'up', ...
  'long_name', 'Altitude of antenna above nominal sea level (WGS84)'};

netcdf_param(end+1).mat_name = 'Latitude';
netcdf_param(end).cdf_name = 'lat';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'standard_name','latitude', ...
  'units', 'degrees_north', ...
  'long_name', 'Latitude of sample', ...
  'axis', 'Y'};

netcdf_param(end+1).mat_name = 'Longitude';
netcdf_param(end).cdf_name = 'lon';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'standard_name','longitude', ...
  'units', 'degrees_east', ...
  'long_name', 'Longitude of sample', ...
  'axis', 'X'};

netcdf_param(end+1).mat_name = 'Data';
netcdf_param(end).cdf_name = 'amplitude';
netcdf_param(end).dim_name = {'fasttime_dim' 'time'};
netcdf_param(end).eval = @(x) (10*log10(x));
netcdf_param(end).attributes = {'units', 'counts in dB', ...
  'long_name', 'Amplitude of low/high gain merged radar reflection after processing', ...
  'coordinates', 'time travel channel'};

netcdf_param(end+1).mat_name = 'Heading';
netcdf_param(end).cdf_name = 'heading';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).eval = @(x) (x*180/pi);
netcdf_param(end).attributes = {'standard_name','heading', ...
  'units', 'degrees', ...
  'long_name', 'Heading of the platform.  Positive is clockwise from above.  Zero is true north.'};

netcdf_param(end+1).mat_name = 'Pitch';
netcdf_param(end).cdf_name = 'pitch';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).eval = @(x) (x*180/pi);
netcdf_param(end).attributes = {'standard_name','pitch', ...
  'units', 'degrees', ...
  'long_name', 'Pitch of the platform.  Positive is nose up.  Zero is horizontal.'};

netcdf_param(end+1).mat_name = 'Roll';
netcdf_param(end).cdf_name = 'roll';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).eval = @(x) (x*180/pi);
netcdf_param(end).attributes = {'standard_name','roll', ...
  'units', 'degrees', ...
  'long_name', 'Roll of the platform.  Positive is right wing down.  Zero is horizontal.'};

netcdf_param(end+1).mat_name = 'Surface';
netcdf_param(end).cdf_name = 'Surface';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'seconds', ...
  'long_name', 'Two way travel time to surface used during processing. Not the final picked surface.'};

netcdf_param(end+1).mat_name = 'Bottom';
netcdf_param(end).cdf_name = 'Bottom';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'seconds', ...
  'long_name', 'Two way travel time to bottom used during processing. Not the final picked bottom.'};

netcdf_param(end+1).mat_name = 'Elevation_Correction';
netcdf_param(end).cdf_name = 'Elevation_Correction';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'range bins', ...
  'long_name', 'Represents the number of zeros that were inserted during elevation compensation for each range line to simulate near-level flight. These zeros are not included in the truncation noise statistics.'};

netcdf_param(end+1).mat_name = 'Truncate_Bins';
netcdf_param(end).cdf_name = 'Truncate_Bins';
netcdf_param(end).dim_name = {'fasttime_dim'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Indices into the original fasttime vector for which the data "amplitude" are available.'};

netcdf_param(end+1).mat_name = 'Truncate_Mean';
netcdf_param(end).cdf_name = 'Truncate_Mean';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Represents a mean of the noise power for the truncated range bins before the surface return. When no range bins were truncated before the surface return the value is NaN.'};

netcdf_param(end+1).mat_name = 'Truncate_Median';
netcdf_param(end).cdf_name = 'Truncate_Median';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Represents a median of the noise power for the truncated range bins before the surface return. When no range bins were truncated before the surface return the value is NaN.'};

netcdf_param(end+1).mat_name = 'Truncate_Std_Dev';
netcdf_param(end).cdf_name = 'Truncate_Std_Dev';
netcdf_param(end).dim_name = {'time'};
netcdf_param(end).attributes = {'units', 'n/a', ...
  'long_name', 'Represents a standard deviation of the noise power for the truncated range bins before the surface return. When no range bins were truncated before the surface return the value is NaN.'};

% Control parameters for supplement netcdf file
supplement_netcdf_param = [];
supplement_netcdf_param(end+1).mat_name = 'coh_noise_removal_artifact';
supplement_netcdf_param(end).cdf_name = 'coh_noise_removal_artifact';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Coherent noise removal artifacts in radar echogram'};
supplement_netcdf_param(end+1).mat_name = 'deconvolution_artifact';
supplement_netcdf_param(end).cdf_name = 'deconvolution_artifact';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Deconvolution artifacts (e.g. sidelobes) in radar echogram'};
supplement_netcdf_param(end+1).mat_name = 'vertical_stripes_artifact';
supplement_netcdf_param(end).cdf_name = 'vertical_stripes_artifact';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Vertical stripes or raised noise floor artifacts in radar echogram'};
supplement_netcdf_param(end+1).mat_name = 'missing_data';
supplement_netcdf_param(end).cdf_name = 'missing_data';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Radar echogram is missing data because radar range gate clips echogram'};
supplement_netcdf_param(end+1).mat_name = 'no_good_data';
supplement_netcdf_param(end).cdf_name = 'no_good_data';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Radar echogram has no good data'};
supplement_netcdf_param(end+1).mat_name = 'low_SNR';
supplement_netcdf_param(end).cdf_name = 'low_SNR';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Low signal to noise ratio'};
supplement_netcdf_param(end+1).mat_name = 'unclassified_artifact';
supplement_netcdf_param(end).cdf_name = 'unclassified_artifact';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Unclassified artifacts exist in radar echogram'};
supplement_netcdf_param(end+1).mat_name = 'land_ice';
supplement_netcdf_param(end).cdf_name = 'land_ice';
supplement_netcdf_param(end).dim_name = {'unit'};
supplement_netcdf_param(end).attributes = {'units' 'mask', ...
  'long_name', 'Land ice, ice shelf, or ice berg contained in echogram'};

param = param_override;

if strcmpi(param.sched.type,'custom_torque')
  global ctrl; % Make this global for convenience in debugging
  ctrl = torque_new_batch(param);
  fprintf('Torque batch: %s\n', ctrl.batch_dir);
  torque_compile('nsidc_delivery_script_task.m',ctrl.sched.hidden_depend_funs,ctrl.sched.force_compile);
end

%% Main loop for each segment
params = read_param_xls(params_fn);
for param_idx = 1:length(params)
  
  param = params(param_idx);
  param = merge_structs(param, param_override);
  if ~isempty(skip_phrase) ...
      && ~isempty(strfind(lower(param.cmd.notes),skip_phrase)) ...
      || ~params(param_idx).cmd.generic
    continue;
  end
  
  param.nsidc.USER_SPECIFIED_DIRECTORY_BASE = USER_SPECIFIED_DIRECTORY_BASE;
  param.nsidc.L1B_cmd = L1B_cmd;
  param.nsidc.data_dir_L1 = data_dir_L1;
  param.nsidc.data_dir_L1_extra = data_dir_L1_extra;
  param.nsidc.image_extra = image_extra;
  param.nsidc.L1B_supplement_cmd = L1B_supplement_cmd;
  param.nsidc.L1B_supplement_name = L1B_supplement_name;
  param.nsidc.L1B_supplement_name_extra = L1B_supplement_name_extra;
  param.nsidc.L2_cmd = L2_cmd;
  param.nsidc.data_dir_L2 = data_dir_L2;
  param.nsidc.premet_param_L1B = premet_param_L1B;
  param.nsidc.premet_param_L2 = premet_param_L2;
  param.nsidc.mcf_version_id_L1B = mcf_version_id_L1B;
  param.nsidc.mcf_version_id_L2 = mcf_version_id_L2;
  param.nsidc.frm_types = frm_types;
  
  param.nsidc.netcdf_param = netcdf_param;
  param.nsidc.supplement_netcdf_param = supplement_netcdf_param;

  % =================================================================
  % Execute tasks/jobs
  fh = @nsidc_delivery_script_task;
  arg{1} = param;
  
  if strcmp(param.sched.type,'custom_torque')
    create_task_param.conforming = true;
    create_task_param.notes = sprintf('%s', param.day_seg);
    ctrl = torque_create_task(ctrl,fh,1,arg,create_task_param);
    
  else
    fprintf('%s (%s)', param.day_seg, datestr(now));
    [success] = fh(arg{1});
  end
  
end

if strcmpi(param.sched.type,'custom_torque')
  % Wait until all submitted jobs to complete
  ctrl = torque_rerun(ctrl);
  if ~all(ctrl.error_mask == 0)
    if ctrl.sched.stop_on_fail
      torque_cleanup(ctrl);
      error('Not all jobs completed, but out of retries (%s)', datestr(now));
    else
      warning('Not all jobs completed, but out of retries (%s)', datestr(now));
      keyboard;
    end
  else
    fprintf('Jobs completed (%s)\n\n', datestr(now));
  end
  torque_cleanup(ctrl);
end

return;
