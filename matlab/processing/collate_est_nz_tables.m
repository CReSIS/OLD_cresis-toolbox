% function collate_est_nz_tables(params_en,param_override)
%
% function collate_est_nz
%
% Collates results from run_collate_est_nz to create a ldap and table of
% unique coherent noise waveforms for entire season.
%
% General order for processing:
% run_collate_est_nz, run_collate_est_nz_tables, run_collate_est_nz_records
%
% Authors: John Paden, Hara Madhav Talasila

%% General Setup
% =========================================================================

param = param_override;
%
% fprintf('=============================================================\n');
% fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
% fprintf('=============================================================\n');

%% Input Checks
% =========================================================================
if ~isfield(param,'tables') || isempty(param.tables)
  param.tables = [];
end

if ~isfield(param.tables, 'enable_visible_plot')
  param.tables.enable_visible_plot = 0;
end

param.tables.in_dir = fileparts(ct_filename_ct_tmp(params_en{1},'','collate_est_nz',''));
if ~exist(param.tables.in_dir, 'dir')
  fprintf('Empty directory --> run_analysis results (%s)\n',param.tables.in_dir);
  return;
end

if ~isfield(param.(mfilename),'debug_out_dir') || isempty(param.(mfilename).debug_out_dir)
  param.tables.debug_out_dir = mfilename;
end

param.tables.out_dir = fileparts(ct_filename_ct_tmp(params_en{1},'',param.tables.debug_out_dir,''));
if ~exist(param.tables.out_dir, 'dir')
  mkdir(param.tables.out_dir);
end

%% Collate the results
% =========================================================================
%%% Let this code thinks there is only one image [1 1];
wf = 1;
adc = 1;

ldap = struct();
ldap.day_seg = {};
ldap.set = [];
ldap.coh_wf = [];
ldap.twtt = [];
ldap.table_idx = [];
ldap.nyquistzone_hw=[];
ldap_idx = 0;

for idx = 1:length(params_en)
  fprintf('%d >> %s\n',idx,params_en{idx}.day_seg);
  %% Load the nz_est and records file
  % =====================================================================
  %
  try
    fn = fullfile(param.tables.in_dir, sprintf('collate_est_nz_wf_%d_adc_%d_%s.mat',wf,adc,params_en{idx}.day_seg));
    tmp = load(fn);
    for set_idx = 1:tmp.coh_wf.sets
      ldap_idx = ldap_idx+1;
      ldap(ldap_idx).day_seg = params_en{idx}.day_seg;
      ldap(ldap_idx).set = set_idx;
      ldap(ldap_idx).coh_wf = tmp.coh_wf.coh_noise(:,set_idx);
      ldap(ldap_idx).twtt = tmp.coh_wf.twtt(:,set_idx);
      ldap(ldap_idx).nyquistzone_hw = NaN;
    end
    clear tmp;
  catch
    fprintf('Missing nz_est file: %s\n',fn);
    continue;
  end
end

%%
% ldap(1).table_idx = [];

table = struct();
% table.xcorr = [];
% table.coh_wf = [];
% % table.nyquistzone = [];
% table.old_day_seg = [];
% table.cur_day_seg = [];
% table.twtt = [];
% table.old_set = [];
% table.cur_set = [];
% table.ldap_idx = [];

table(1).xcorr = 1;
table(1).coh_wf = ldap(1).coh_wf;
% table.nyquistzone(1) = [];
table(1).old_day_seg = ldap(1).day_seg;
table(1).cur_day_seg = ldap(1).day_seg;
table(1).twtt = ldap(1).twtt;
table(1).old_set = 1;
table(1).cur_set = 1;
table(1).ldap_idx = 1;
table(1).nyquistzone_hw=NaN;

xcorr_norm=[];
for idx = 1:length(ldap)
  
  for t_idx = 1:length(table)
    xcorr_norm(idx,t_idx) = abs( sum( (ldap(idx).coh_wf-mean(ldap(idx).coh_wf)) .* conj( (table(t_idx).coh_wf-mean(table(t_idx).coh_wf)) ) ) / ( std(ldap(idx).coh_wf)*std(table(t_idx).coh_wf) ) /length(table(t_idx).coh_wf) ); %Calculate xcorr_norm
  end
  [max_val, max_idx] = max(xcorr_norm(idx,:));
  if abs(max_val - table(max_idx).xcorr) < 0.15
    if max_val > table(max_idx).xcorr || table(max_idx).xcorr ==1
      table(max_idx).xcorr = max_val;
      table(max_idx).coh_wf = ldap(idx).coh_wf;
      table(max_idx).old_day_seg = table(max_idx).cur_day_seg;
      table(max_idx).cur_day_seg = ldap(idx).day_seg;
      table(max_idx).twtt = ldap(idx).twtt;
      table(max_idx).old_set = table(max_idx).cur_set;
      table(max_idx).cur_set = ldap(idx).set;
      table(max_idx).ldap_idx = idx;
      table(max_idx).nyquistzone_hw=NaN;
    end
    ldap(idx).table_idx = max_idx;
  else
    t_idx = t_idx+1;
    table(t_idx).xcorr = 1;
    table(t_idx).coh_wf = ldap(idx).coh_wf;
    table(t_idx).old_day_seg = ldap(idx).day_seg;
    table(t_idx).cur_day_seg = ldap(idx).day_seg;
    table(t_idx).twtt = ldap(idx).twtt;
    table(t_idx).old_set = table(t_idx-1).cur_set;
    table(t_idx).cur_set = 1;
    table(t_idx).ldap_idx = idx;
    table(max_idx).nyquistzone_hw=NaN;
    ldap(idx).table_idx = t_idx;
  end
end


figure(99);
clf(99);
leg=[];
for idx = 1:length(table)
  plot(lp(table(idx).coh_wf));hold on;
  leg = [leg {sprintf('twtt %.2f us xcorr(%.4f) for %s [%d] and %s [%d] ',table(idx).twtt/1e-6, table(idx).xcorr, strrep(table(idx).cur_day_seg,'_','-'), table(idx).cur_set, strrep(table(idx).old_day_seg,'_','-'),table(idx).old_set) } ];
end
legend(leg);
hold off;


user_input = input(sprintf('Enter %d nyquist zone estimates example: [1 2 1 0 3]\n',length(table)));



if length(user_input) == length(table)
  for idx=1:length(table)
    table(idx).nyquistzone_hw = user_input(idx);
    [ldap([ldap.table_idx]==idx).nyquistzone_hw] = deal(user_input(idx));
  end
end


if 0 %display ldap
  fprintf('day_seg\t\t[nz_hw');
  prev = [];
  
  for idx = 1:length(ldap)
    if ~strcmpi(ldap(idx).day_seg,prev)
      fprintf(']\n%s\t[%d',ldap(idx).day_seg,ldap(idx).nyquistzone_hw);
    else
      fprintf(' %d',ldap(idx).nyquistzone_hw);
    end
    prev = ldap(idx).day_seg;
  end
  fprintf(']\n');
end

for idx = 1:length(params)
  fprintf('%s\t%d\t',params(idx).day_seg(1:8),str2double(params(idx).day_seg(10:11)));
  [indices] = find(strcmp({ldap.day_seg},params(idx).day_seg)==1);
  if ~isempty(indices)
    fprintf('[ ');
    for idxp = 1:length(indices)
      fprintf('%d ',ldap(indices(idxp)).nyquistzone_hw);
    end
    fprintf(']\t new\n');
  else
    fprintf('[ ');
    for idxp = 1:length(params(idx).collate_nz_est.nz_table)
      fprintf('%d ',params(idx).collate_nz_est.nz_table(idxp));
    end
    fprintf(']\t old\n');
  end
end
  
  
  %     if ~param.collate_est_nz.enable_visible_plot
  %       try
  %         delete(cur_fig);
  %       catch
  %       end
  %     end
  
  
  
