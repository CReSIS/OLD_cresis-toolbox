% script compare_surfData
%
% Compares layer A (reference) and layer B
%   from a surfData file.
%
% Authors: Victor Berger, Mohanad Al-Ibadi

%%%%%%%% LAYER A (REFERENCE):
surfdata_ref       = 'CSA_music_surfData';
layer_name_ref     = 'bottom';
%%%%%%%% LAYER B (OTHER):
surfdata_other     = 'CSA_music_surfData_no_QC';
layer_name_other   = 'bottom';
%%%%%%%% BOTH lAYERS:
param.radar_name   = 'rds';
param.season_name  = '2014_Greenland_P3';
segs               = {'20140401_03','20140506_01','20140325_05','20140325_06','20140325_07'};
param.slices       = []; % Leave empty to run all slices
cutoffs            = [0 5 10 15 20 25]; % Cutoff interval for statistics
DOA_trim           = 3;

% Output directory of CDF image (cumulative distribution function)
out_dir = '';

display_status_flag = false;
save_cdf_flag       = false;


%%%%%%%%%%%%%%%%%%%%%
%%% Automated section
rmse_f  = [];
diff_f = [];
med_f  = [];
max_df = [];
min_df = [];
total_diff = [];
counter_f  = 1;

if isempty(out_dir) && save_cdf_flag
  error('Please enter output directory for CDF image.');
end

for seg_idx = 1:length(segs)
  param.day_seg = segs{seg_idx};
  fn = fullfile(ct_filename_out(param,surfdata_ref,''));
  d = dir(fn);
  num_frms = length(d)-2; % Number of frames
  
  frames = [];
  for idx = 1:num_frms
    tmp = d(idx+2).name;
    frames(idx) = str2num(tmp(18:end-4)); % Exact frame numbers
  end
  
  for frame_idx = frames
    fprintf('\nWorking on frame %s_%03d. ',param.day_seg, frame_idx);
    fn_A = fullfile(ct_filename_out(param,surfdata_ref,''),sprintf('Data_%s_%03d.mat',param.day_seg,frame_idx));
    fn_B = fullfile(ct_filename_out(param,surfdata_other,''),sprintf('Data_%s_%03d.mat',param.day_seg,frame_idx));
    
    if exist(fn_A, 'file') ~= 2
      error('surfData file for LAYER A not found, check parameters.');
    elseif exist(fn_B, 'file') ~= 2
      error('surfData file for LAYER B not found, check parameters.');
    end
    
    fprintf('Loading data...');
    data_A = load(fn_A);
    data_B = load(fn_B);
    fprintf(' done.');
    
    idx_A = strcmp({data_A.surf.name},layer_name_ref);
    idx_B = strcmp({data_B.surf.name},layer_name_other);
    
    if ~any(idx_A)
      error('Could not find desired layer name in LAYER A');
    elseif ~any(idx_B)
      error('Could not find desired layer name in LAYER B');
    end
    
    idx_A = find(idx_A);
    idx_B = find(idx_B);
    
    if(isempty(param.slices))
      slices = [1 : size(data_A.surf(idx_A).y(:,:),2)];
    else
      slices = param.slices;
    end
    
    layer_diff  = abs(data_B.surf(idx_B).y(1+DOA_trim:end-DOA_trim+1,slices) - data_A.surf(idx_A).y(1+DOA_trim:end-DOA_trim+1,slices));
    rmse        = sqrt(mean(abs(layer_diff(:)).^2));
    mean_diff   = nanmean(layer_diff(:));
    median_diff = nanmedian(layer_diff(:));
    min_diff    = nanmin(layer_diff(:));
    max_diff    = nanmax(layer_diff(:));
    total_diff  = cat(2, total_diff, layer_diff(:)');
    
    rmse_f(counter_f) = rmse;
    diff_f(counter_f) = mean_diff;
    med_f(counter_f)  = median_diff;
    min_df(counter_f) = min_diff;
    max_df(counter_f) = max_diff;
    counter_f         = counter_f + 1;
    
    if(~any(layer_diff))
      fprintf('\n Both layers are exactly equal, no error.\n')
      continue;
    end
    
    if display_status_flag
      fprintf('\n\nBetween slices %d:%d of frame %d:',slices(1),slices(end),frame_idx);
      fprintf('\n Root mean squared error: %.4f', rmse);
      fprintf('\n Mean difference: %.4f', mean_diff);
      fprintf('\n Median difference: %.4f', median_diff);
      fprintf('\n Minimum difference: %.4f', min_diff);
      fprintf('\n Maximum difference: %.4f\n', max_diff);
      imagesc(layer_diff); colorbar;
    end
  end
end

% STATISTICS: between frames
fprintf('\n\nSTATISTICS: over all frames \n');
fprintf('\n  Mean RMSE over all frames: %.4f\n', nanmean(rmse_f));
fprintf('\n  Median RMSE over all frames: %.4f\n', nanmedian(rmse_f));
fprintf('\n  Average mean difference: %.4f\n', mean(diff_f));
fprintf('\n  Average median difference: %.4f\n', mean(med_f));
fprintf('\n  Minimum mean difference: %.4f\n', min(min_df));
fprintf('\n  Maximum mean difference: %.4f', max(max_df));

% STATISTICS, make sure 'cutoffs' vector is properly set
dl = length(total_diff);
fprintf('\n\nSTATISTICS: from the CDF\n');
for i = cutoffs
  hits = length(find(abs(total_diff) <= i));
  fprintf('\n  Within %d bins of the reference: %d out of %d, or %.2f%%\n', ...
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
if save_cdf_flag
  fig = gcf;
  fig.PaperPositionMode = 'auto';% This property prevents saved figure from being distorted
  out_fn_name = sprintf('cdf_layer_errors');
  saveas(f,[fullfile(out_dir,out_fn_name),'.fig']);
  saveas(f,[fullfile(out_dir,out_fn_name),'.jpg']);
end