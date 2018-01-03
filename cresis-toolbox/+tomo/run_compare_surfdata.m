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
ct_set_params(params,'cmd.generic',0);
ct_set_params(params,'cmd.generic',1,'20140401_03|20140506_01|20140325_05|20140325_06|20140325_07');

% surfdata_ref: Directory where reference (ideal) layer exists
param_override.compare.surfdata_ref = 'surfData_v2';
% layer_name_ref: Name of layer in reference file
param_override.compare.layer_name_ref = 'bottom';

% surfdata_other: Directory where layer to compare exists
param_override.compare.surfdata_other = 'surfData_v2_no_MC';
% layer_name_other: Name of layer comparison file
param_override.compare.layer_name_other = 'bottom';

% compare.cutoffs: Cutoff interval in % for statistics;
param_override.compare.cutoffs = [0 5 10 15 20 25]; 

% compare.DOA_trim: DOA bins to ignore in the comparison. First number is
%   number of columns to trim from start and second number is number of
%   columns to trim from the end of the DOA bin list.
param_override.compare.DOA_trim = [3 3];

% compare.out_dir: Output directory to store CDF image (cumulative
%   distribution function).
param_override.compare.out_dir = ''; 
param_override.compare.save_cdf_flag = false;

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

for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  [rmse_f,diff_f,med_f,min_f,max_f,total_diff] = tomo.compare_surfdata(param,param_override);
  
end

%% Collate results

% STATISTICS: between frames
fprintf('\n\nSTATISTICS: over all frames \n');
fprintf('  Mean RMSE over all frames: %.4f\n', nanmean(rmse_f));
fprintf('  Median RMSE over all frames: %.4f\n', nanmedian(rmse_f));
fprintf('  Average mean difference: %.4f\n', mean(diff_f));
fprintf('  Average median difference: %.4f\n', mean(med_f));
fprintf('  Minimum mean difference: %.4f\n', min(min_f));
fprintf('  Maximum mean difference: %.4f', max(max_f));

% STATISTICS, make sure 'cutoffs' vector is properly set
dl = length(total_diff);
fprintf('\n\nSTATISTICS: from the CDF\n');
for i = cutoffs
  hits = length(find(abs(total_diff) <= i));
  fprintf('  Within %d bins of the reference: %d out of %d, or %.2f%%\n', ...
    i, hits, dl, ((100 * hits) / dl));
end

% Plot the CDF
total_diff = sort(total_diff);
[cdf_vals,cdf_pts] = ecdf(total_diff);
desired_idx = find(cdf_pts==25,1);

f = figure(11);
plot(cdf_pts(1:desired_idx),cdf_vals(1:desired_idx),'LineWidth',1.5)
xlabel('|Error|(range-bins)')
ylabel('F(|Error|)')
title('cdf(|Error|)')
ylim([0 1])
grid on

% Save
if param_override.compare.save_cdf_flag
  fig = gcf;
  fig.PaperPositionMode = 'auto';% This property prevents saved figure from being distorted
  
  out_fn_name = sprintf('cdf_layer_errors');
  saveas(f,[fullfile(param_override.compare.out_dir,out_fn_name),'.fig']);
  saveas(f,[fullfile(param_override.compare.out_dir,out_fn_name),'.jpg']);
end
