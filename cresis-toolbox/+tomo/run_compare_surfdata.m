% script compare_surfData
%
% Compares layer A (reference) and layer B
%   from a surfData file.
%
% Authors: Victor Berger, Mohanad Al-Ibadi, John Paden
%
% See also: tomo.run_compare_surfdata.m, tomo.compare_surfdata,
%   tomo.surfdata.m

%% Setup

params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'');
params = ct_set_params(params,'cmd.generic',0);
params = ct_set_params(params,'cmd.generic',1,'20140401_03|20140506_01|20140325_05|20140325_06|20140325_07');
params = ct_set_params(params,'cmd.frms',[]);
% params = ct_set_params(params,'cmd.generic',1,'20140401_03');
% params = ct_set_params(params,'cmd.frms',[37 39 43:47]);
% params = ct_set_params(params,'cmd.generic',0);

% surfdata_ref: Directory where reference (ideal) layer exists
param_override.compare.surfdata_ref = 'surfData_IGARSS';
% layer_name_ref: Name of layer in reference file
% param_override.compare.layer_name_ref = 'bottom';
% param_override.compare.cdf_title = 'TRW-S mod';
% param_override.compare.layer_name_ref = 'bottom extract';
% param_override.compare.cdf_title = 'TRW-S';
% param_override.compare.layer_name_ref = 'bottom detect';
% param_override.compare.cdf_title = 'Viterbi';
param_override.compare.layer_name_ref = 'bottom viterbi';
param_override.compare.cdf_title = 'Viterbi mod';

% surfdata_ref: Directory where reference (ideal) layer exists
param_override.compare.surfdata_quality = 'surfData_IGARSS';
% layer_name_ref: Name of layer in reference file
param_override.compare.layer_name_quality = 'bottom quality';

% surfdata_other: Directory where layer to compare exists
param_override.compare.surfdata_other = 'surfData';
% layer_name_other: Name of layer comparison file
param_override.compare.layer_name_other = 'bottom';

% compare.cutoffs: Cutoff interval in % for statistics;
param_override.compare.cutoffs = [1 5 10 15 20 25]; 

% compare.DOA_trim: DOA bins to ignore in the comparison. First number is
%   number of columns to trim from start and second number is number of
%   columns to trim from the end of the DOA bin list.
param_override.compare.DOA_trim = [3 4];

% compare.out_dir: Output directory to store CDF image (cumulative
%   distribution function).
param_override.compare.out_dir = '~/'; 
param_override.compare.save_cdf_flag = true;

% compare.display_status_flag: logical, set to true to display results for
%   each frame as the comparisons are done
param_override.compare.display_status_flag = false;


%% Automated section
global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

collate_results = false;
frm_id_f = {};
rmse_f = [];
diff_f = [];
med_f = [];
min_f = [];
max_f = [];
total_diff = [];
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  [frm_id_s,rmse_s,diff_s,med_s,min_s,max_s,total_diff_s] = tomo.compare_surfdata(param,param_override);
  frm_id_f = cat(2,frm_id_f,frm_id_s);
  rmse_f = cat(2,rmse_f,rmse_s);
  diff_f = cat(2,diff_f,diff_s);
  med_f = cat(2,med_f,med_s);
  min_f = cat(2,min_f,min_s);
  max_f = cat(2,max_f,max_s);
  total_diff = cat(2,total_diff,total_diff_s);
  collate_results = true;
end

%% Collate results
if collate_results
  % STATISTICS: between frames
  fprintf('\n\nSTATISTICS: over all frames \n');
  fprintf('Mean RMSE\t%.4f\n', nanmean(rmse_f));
  fprintf('Median RMSE\t%.4f\n', nanmedian(rmse_f));
  fprintf('Mean\t%.4f\n', mean(diff_f));
  fprintf('Median\t%.4f\n', median(med_f));
  fprintf('Minimum\t%.4f\n', min(min_f));
  fprintf('Maximum\t%.4f\n', max(max_f));
  
  % STATISTICS, make sure 'cutoffs' vector is properly set
  dl = numel(total_diff);
  fprintf('\n\nSTATISTICS: from the CDF\n');
  for i = param_override.compare.cutoffs
    hits = length(find(abs(total_diff) <= i));
    fprintf('  Within %d bins of the reference: %d out of %d, or %.2f%%\n', ...
      i, hits, dl, ((100 * hits) / dl));
  end
  
  % Plot the CDF
  [cdf_vals,cdf_pts] = ecdf(sort(total_diff(:)));
  desired_idx = find(cdf_pts==25,1);

  if ~exist('h_comp_fig','var') || ~isvalid(h_comp_fig)
    h_comp_fig = figure;
    h_comp_leg = {};
    h_axes = axes('parent',h_comp_fig);
    pos = get(h_comp_fig,'Position');
    set(h_comp_fig,'Position',[pos(1:2)   471   230])
  end
  hold(h_axes,'on');
  plot(cdf_pts(1:desired_idx),cdf_vals(1:desired_idx),'LineWidth',1.5,'Parent',h_axes)
  xlabel('|Error| (range-bins)','Parent',h_axes)
  ylabel('F(|Error|)','Parent',h_axes)
  %title('cdf(|Error|)','Parent',h_axes)
  ylim(h_axes,[0 1])
  grid(h_axes,'on');
  h_comp_leg{end+1} = param_override.compare.cdf_title;
  legend(h_axes,h_comp_leg,'Location','southeast');
  
  % Save
  if param_override.compare.save_cdf_flag
    set(h_comp_fig,'PaperPositionMode','auto');% This property prevents saved figure from being distorted
    
    out_fn_name = sprintf('cdf_layer_errors');
    saveas(h_comp_fig,[fullfile(param_override.compare.out_dir,out_fn_name),'.fig']);
    saveas(h_comp_fig,[fullfile(param_override.compare.out_dir,out_fn_name),'.jpg']);
  end
end
