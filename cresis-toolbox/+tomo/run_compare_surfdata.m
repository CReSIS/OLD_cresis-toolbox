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

% params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'');
% params = ct_set_params(params,'cmd.generic',0);
% params = ct_set_params(params,'cmd.generic',1,'day_seg','20140401_03|20140506_01|20140325_05|20140325_06|20140325_07');
% % params = ct_set_params(params,'cmd.generic',1,'day_seg','20140401_03');
% params = ct_set_params(params,'cmd.frms',[]);
% % params = ct_set_params(params,'cmd.generic',1,'20140401_03');
% % params = ct_set_params(params,'cmd.frms',[37]);
% % params = ct_set_params(params,'cmd.generic',0);



  params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'),'','post');
%   params = ct_set_params(params,'cmd.generic',0);
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20140325_05|20140325_06|20140325_07|20140401_03|20140506_01');
%   params = ct_set_params(params,'cmd.generic',1,'day_seg','20140401_03');
%   params = ct_set_params(params,'cmd.frms',[20]);











% surfdata_ref: Directory where reference (ideal) layer exists
% compare_params.surfdata_ref = 'old/surfData';
compare_params.surfdata_ref = 'old/surfData';
% surf_name_ref: Name of layer in reference file
compare_params.surf_name_ref = 'bottom';

% surfdata_ref: Directory where reference (ideal) layer exists
compare_params.surfdata_quality = 'old/surfData';
% surf_name_ref: Name of layer in reference file
compare_params.layer_name_quality = 'bottom quality';

% surfdata_other: Directory where layer to compare exists
compare_params.surfdata_other = {};
% surf_name_other: Name of layer comparison file
compare_params.surf_name_other = {};
% cdf_title: Title or label to use in plots
compare_params.cdf_title = {};

% Victor's IGARSS18 paper
% compare_params.surfdata_other{end+1} = 'shane_music_surfData';
% compare_params.surf_name_other{end+1} = 'bottom detect';
% compare_params.cdf_title{end+1} = 'Viterbi';
% compare_params.surfdata_other{end+1} = 'shane_music_surfData';
% compare_params.surf_name_other{end+1} = 'bottom extract';
% compare_params.cdf_title{end+1} = 'TRW-S';
compare_params.surfdata_other{end+1} = 'surfdata_Victor_TEST';
compare_params.surf_name_other{end+1} = 'bottom viterbi';
compare_params.cdf_title{end+1} = 'Viterbi Mod';
compare_params.surfdata_other{end+1} = 'surfdata_Victor_TEST';
compare_params.surf_name_other{end+1} = 'bottom trws';
compare_params.cdf_title{end+1} = 'TRW-S Mod';

% Mohanad's RadConf18 paper
% compare_params.surfdata_other{end+1} = 'surfData_no_MC_test';
% compare_params.surfdata_other{end+1} = 'surfData_no_MC';
% compare_params.surf_name_other{end+1} = 'bottom';
% compare_params.cdf_title{end+1} = 'TRW-S/Vitirbi';

% compare.cutoffs: Cutoff interval in % for statistics;
compare_params.cutoffs = [1 5 10 15 20 25];

% compare.DOA_trim: DOA bins to ignore in the comparison. First number is
%   number of columns to trim from start and second number is number of
%   columns to trim from the end of the DOA bin list.
compare_params.DOA_trim = [3 4];

% compare.out_dir: Output directory to store CDF image (cumulative
%   distribution function).
compare_params.out_dir = '~/';
compare_params.save_cdf_flag = true;

% compare.display_status_flag: logical, set to true to display results for
%   each frame as the comparisons are done
compare_params.display_status_flag = false;


%% Automated section
global gRadar;

param_override = struct('compare',compare_params);

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

collate_results = false;
for param_idx = 1:length(params)
  param = params(param_idx);
  if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
    continue;
  end
  
  [frm_id_s,rmse_s,mean_s,median_s,min_s,max_s,total_diff_s] = tomo.compare_surfdata(param,param_override);
  if ~collate_results
    collate_results = true;
    frm_id_f = frm_id_s;
    rmse_f = rmse_s;
    mean_f = mean_s;
    median_f = median_s;
    min_f = min_s;
    max_f = max_s;
    total_diff = total_diff_s;
  else
    frm_id_f = cat(2,frm_id_f,frm_id_s);
    for surf_idx = 1:length(rmse_s)
      rmse_f{surf_idx} = cat(2,rmse_f{surf_idx},rmse_s{surf_idx});
      mean_f{surf_idx} = cat(2,mean_f{surf_idx},mean_s{surf_idx});
      median_f{surf_idx} = cat(2,median_f{surf_idx},median_s{surf_idx});
      min_f{surf_idx} = cat(2,min_f{surf_idx},min_s{surf_idx});
      max_f{surf_idx} = cat(2,max_f{surf_idx},max_s{surf_idx});
      total_diff{surf_idx} = cat(2,total_diff{surf_idx},total_diff_s{surf_idx});
    end
  end
end

%% Collate results
if collate_results
  % STATISTICS: between frames
  
  fprintf('\n\nSTATISTICS: over all segments and frames \n');
  
  rmse_f_all = [];
  for surf_idx = 1:length(rmse_f)
    rmse_f_all(surf_idx) = sqrt(nanmean(rmse_f{surf_idx}.^2));
    fprintf('RMSE\t%s\t%.1f\n', compare_params.cdf_title{surf_idx}, rmse_f_all(surf_idx));
  end
  
  mean_f_all = [];
  for surf_idx = 1:length(mean_f)
    mean_f_all(surf_idx) = nanmean(mean_f{surf_idx});
    fprintf('Mean\t%s\t%.1f\n', compare_params.cdf_title{surf_idx}, mean_f_all(surf_idx));
  end
  
  median_f_all = [];
  for surf_idx = 1:length(median_f)
    median_f_all(surf_idx) = nanmedian(median_f{surf_idx});
    fprintf('Median\t%s\t%.1f\n', compare_params.cdf_title{surf_idx}, median_f_all(surf_idx));
  end
  
  min_f_all = [];
  for surf_idx = 1:length(min_f)
    min_f_all(surf_idx) = nanmin(min_f{surf_idx});
    fprintf('Min\t%s\t%.1f\n', compare_params.cdf_title{surf_idx}, min_f_all(surf_idx));
  end
  
  max_f_all = [];
  for surf_idx = 1:length(max_f)
    max_f_all(surf_idx) = nanmax(max_f{surf_idx});
    fprintf('Max\t%s\t%.1f\n', compare_params.cdf_title{surf_idx}, max_f_all(surf_idx));
  end
  
  if ~exist('h_comp_fig','var') || ~isvalid(h_comp_fig)
    h_comp_fig = figure;
    h_comp_leg = {};
    h_axes = axes('parent',h_comp_fig);
    pos = get(h_comp_fig,'Position');
    set(h_comp_fig,'Position',[pos(1:2)   471   230])
  end
  fprintf('\n\nSTATISTICS: from the CDF\n');
  h_plot = [];
  legend_str = {};
  for surf_idx = 1:length(total_diff)
    % Concatenate total_diff
    total_diff_all = [];
    for frm_idx = 1:length(total_diff{surf_idx})
      total_diff_all = cat(2,total_diff_all,total_diff{surf_idx}{frm_idx});
    end
    
    % STATISTICS, make sure 'cutoffs' vector is properly set
    total_diff_all = total_diff_all(:);
    dl = numel(total_diff_all);
    fprintf('Method\t');
    for cutoff = compare_params.cutoffs
      fprintf('%.0f%%\t', cutoff);
    end
    fprintf('\n');
    fprintf('%s\t', compare_params.cdf_title{surf_idx});
    for cutoff = compare_params.cutoffs
      hits = sum(abs(total_diff_all) <= cutoff);
      fprintf('%.2f%%\t', 100*hits/ dl);
    end
    fprintf('\n');
    
    % Plot the CDF
    [cdf_vals,cdf_pts] = ecdf(sort(total_diff_all));
    desired_idx = find(cdf_pts==25,1);
    
    h_plot(surf_idx) = plot(cdf_pts(1:desired_idx),cdf_vals(1:desired_idx),'LineWidth',1.5,'Parent',h_axes);
    hold(h_axes,'on');
    
    legend_str{surf_idx} = compare_params.cdf_title{surf_idx};
  end
  xlabel('|Error| (range-bins)','Parent',h_axes)
  ylabel('F(|Error|)','Parent',h_axes)
  %title('cdf(|Error|)','Parent',h_axes)
  ylim(h_axes,[0 1])
  grid(h_axes,'on');
  legend(h_axes,h_plot,legend_str,'Location','southeast');
  
%   end
  % Save
  if compare_params.save_cdf_flag
    set(h_comp_fig,'PaperPositionMode','auto');% This property prevents saved figure from being distorted
    
    out_fn_name = fullfile(compare_params.out_dir,sprintf('cdf_layer_errors'));
    if ~exist(compare_params.out_dir,'dir')
      mkdir(compare_params.out_dir);
    end
    saveas(h_comp_fig,[out_fn_name,'.fig']);
    saveas(h_comp_fig,[out_fn_name,'.jpg']);
  end
end
