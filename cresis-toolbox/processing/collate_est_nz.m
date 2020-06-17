function collate_est_nz(param,param_override)
%
% function collate_est_nz
%
% Classifies a day_segment into sets based on normalized cross correlation
% between blocks of coherent noise.
%
% General order for processing:
% run_collate_est_nz, run_collate_est_nz_tables, run_collate_est_nz_records
%
% Authors: John Paden, Hara Madhav Talasila

%% General Setup
% =========================================================================

param = merge_structs(param, param_override);

fprintf('=============================================================\n');
fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
fprintf('=============================================================\n');

%% Input Checks
% =========================================================================
if ~isfield(param,'collate_est_nz') || isempty(param.collate_est_nz)
  param.collate_est_nz = [];
end

if ~isfield(param.collate_est_nz, 'enable_visible_plot')
  param.collate_est_nz.enable_visible_plot = 0;
end

if ~isfield(param.collate_est_nz, 'enable_extralines_on_plot')
  param.collate_est_nz.enable_extralines_on_plot = 0;
end

if ~isfield(param.collate_est_nz, 'enable_verbose')
  param.collate_est_nz.enable_verbose = 0;
end

if ~isfield(param.collate_est_nz, 'reuse_files')
  param.collate_est_nz.reuse_files = 0;
end

if ~isfield(param.collate_est_nz, 'disconnect_mul_factor')
  param.collate_est_nz.disconnect_mul_factor = 2;
end

if ~isfield(param.collate_est_nz, 'disconnect_threshold')
  param.collate_est_nz.disconnect_threshold = 10; %dBm
end

if ~isfield(param.collate_est_nz, 'xcorr_tol')
  param.collate_est_nz.xcorr_tol = 0.01;
end

if ~isfield(param.collate_est_nz, 'xcorr_region_tol')
  param.collate_est_nz.xcorr_region_tol = 0.02;
end

stat_max_dir = ct_filename_out(param,'analysis_max','',1);
if ~exist(stat_max_dir, 'dir')
  fprintf('Empty directory --> run_analysis results (%s)\n',stat_max_dir);
  return;
end

coh_noise_dir = ct_filename_out(param,'analysis','',1);
if ~exist(coh_noise_dir, 'dir')
  fprintf('Empty directory --> run_analysis results (%s)\n',coh_noise_dir);
  return;
end

if ~isfield(param.(mfilename),'debug_out_dir') || isempty(param.(mfilename).debug_out_dir)
  param.collate_est_nz.debug_out_dir = mfilename;
end

param.collate_est_nz.out_dir = fileparts(ct_filename_ct_tmp(param,'',param.collate_est_nz.debug_out_dir,''));
if ~exist(param.collate_est_nz.out_dir, 'dir')
  mkdir(param.collate_est_nz.out_dir);
end

%% Collate the results
% =========================================================================

% For figure handles
fig_count = 0;
for idx =1:length(param.analysis.imgs)
  fig_count = fig_count + size(param.analysis.imgs{idx},1);
end
h_fig = get_figures(fig_count*2,param.collate_est_nz.enable_visible_plot);
fig_idx=0;

% The big loop for each img-wf_adc pair
for img = 1:length(param.analysis.imgs)
  for wf_adc = 1:size(param.analysis.imgs{img},1)
    wf = param.analysis.imgs{img}(wf_adc,1);
    adc = param.analysis.imgs{img}(wf_adc,2);
    fprintf('[[ For wf-adc %d-%d]] (%s)\n',wf,adc,datestr(now));
    %% Load the coh_noise file
    % =====================================================================
    %
    try
      fn = fullfile(coh_noise_dir, sprintf('coh_noise_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc));
      simp = load(fn);
    catch
      fprintf('Missing file: %s\n',fn);
      return;
    end
    
    reuse_success = 0;
    if param.collate_est_nz.reuse_files
      try
        fn = ct_filename_ct_tmp(param,'',param.collate_est_nz.debug_out_dir,sprintf('%s_wf_%d_adc_%d.mat',param.collate_est_nz.debug_out_dir,wf,adc));
        reuse = load(fn);
        reuse_success = 1;
      catch
        fprintf('Missing reuse_file: %s\n',fn);
        reuse_success = 0;
      end
    end
    
    if ~reuse_success
      % Initialization
      simp.coh_noise = [];
      if param.collate_est_nz.enable_verbose; fprintf('[coh_noise] Forming a huge Nt x Nx matrix\n'); end;
      for idx = 1:length(simp.coh_ave)
        simp.coh_noise =[simp.coh_noise, simp.coh_ave{idx}]; % Form a huge Nt x Nx matrix
      end
      [Nt, Nx] = size(simp.coh_noise);
      block_upper_lim = param.analysis.block_size*(1:Nx);
      
      xc=struct(); % xcorr -> Cross convolution
      xc.std = std(simp.coh_noise,1,1,'omitnan');
      xc.mean = mean(simp.coh_noise,1,'omitnan');
      
      if param.collate_est_nz.enable_verbose; fprintf('[xcorr calculation] %d iterations - %s\n',Nx,char(datetime)); end;
      
      if 0 % very slow
        xc.norm = NaN(Nx,Nx);
        for idx = 1:Nx
          for idx2 = idx:Nx
            xc.norm(idx,idx2) = abs( sum( (simp.coh_noise(:,idx)-xc.mean(idx)) .* conj( (simp.coh_noise(:,idx2)-xc.mean(idx2)) ) ) / ( xc.std(idx)*xc.std(idx2) ) /Nt ); %Calculate xcorr_norm
            xc.norm(idx2,idx) = xc.norm(idx,idx2); % form lower triangle
          end
        end
      end
      
      xc.norm = NaN(Nx,Nx);
      olaaa = ( simp.coh_noise - repmat(xc.mean,Nt,1) ) ./ repmat(xc.std,Nt,1);
      for idx = 1:Nx
        for idx2 = idx:Nx
          xc.norm(idx,idx2) = abs( sum( olaaa(:,idx) .* conj(olaaa(:,idx2)) )/Nt);
          xc.norm(idx2,idx) = xc.norm(idx,idx2); % form lower triangle
        end
      end
      
      xc.max = -Inf*ones(1,Nx);
      xc.max_idx = NaN(1,Nx);
      [xc.max, xc.max_idx] = max(xc.norm-eye(Nx));
      if 0
        figure(111);subplot(211);plot(xc.max_idx);title(sprintf('xc.row max idx Mode=%d',mode(xc.max_idx)));
        subplot(212); plot(xc.max); title(sprintf('xc.max'));
      end
      
    else % reuse
      
      simp.coh_noise = reuse.simp.coh_noise;
      [Nt, Nx] = size(simp.coh_noise);
      block_upper_lim = param.analysis.block_size*(1:Nx);
      xc = reuse.xc;
      
    end
    
    %% Detect Disconnect
    % =====================================================================
    % Load the statistics_max file
    try
      fn = fullfile(stat_max_dir, sprintf('stats_%s_wf_%d_adc_%d.mat',param.day_seg,wf,adc));
      dd = load(fn);
    catch
      fprintf('Missing file: %s\n',fn);
      return;
    end
    % Initialization
    rline_pow = [];
    num_blocks=0;
    
    % Collate
    for b_size = 1:length(dd.stats)
      rline_pow = [rline_pow dd.stats{b_size}{1}(1,:)];
      num_blocks = num_blocks + size(dd.stats{b_size}{1}(1,:),2);
    end
    
    % Calculate the Max power for all the blocks
    rline_max = rline_pow.*conj(rline_pow) / 50; %Check Z0=50 Ohm in param sheet
    rline_max_dB(wf,adc,:) = 10*log10(rline_max) + 30; % in dBm
    temp = squeeze(rline_max_dB(wf,adc,:));
    temp2 = temp(~isinf(temp));
    rline_max_dB_mean(wf,adc) = mean(temp2,1,'omitnan');
    rline_max_dB_std(wf,adc) = std(temp2,1,'omitnan');
    rline_max_dB_median(wf,adc) = median(temp2,1,'omitnan');
    
    
    % Main structure elements
    disconnect={};
    disconnect.detected = 0;
    disconnect.records_mask=zeros(1,num_blocks);
    disconnect.blocks = [];
    
    disconnect.limit = rline_max_dB_median(wf,adc) - max( rline_max_dB_std(wf,adc) * param.collate_est_nz.disconnect_mul_factor , param.collate_est_nz.disconnect_threshold );
    disconnect.records = find(rline_max_dB(wf,adc,:)<disconnect.limit).';
    BAD_from_analysis =[];
    
    % %if DISCONNECT detected
    if ~isempty(disconnect.records)
      
      disconnect.records_mask(disconnect.records) = 1;
      disconnect.blocks = find(fir_dec(disconnect.records_mask,param.analysis.block_size))-1;
      
      % Group the continuous blocks
      disconnect.block_sets={};
      disconnect.set_idx = 1; %current set index
      disconnect.prev_set_idx = 1;
      disconnect.block_sets{disconnect.set_idx} = [disconnect.blocks(1)]; %First record is part of the set
      
      for idx=2:length(disconnect.blocks)
        if disconnect.blocks(idx) == disconnect.blocks(idx-1)+1
          if disconnect.prev_set_idx == disconnect.set_idx
            disconnect.block_sets{disconnect.set_idx} = [disconnect.block_sets{disconnect.set_idx} disconnect.blocks(idx)];
          else
            disconnect.set_idx = disconnect.set_idx + 1;
            disconnect.block_sets{disconnect.set_idx} = [disconnect.blocks(idx)];
            disconnect.prev_set_idx = disconnect.set_idx;
          end
        else
          disconnect.set_idx = disconnect.set_idx + 1;
          disconnect.block_sets{disconnect.set_idx} = [disconnect.blocks(idx)];
        end
      end
      
      % Number of disconnects after grouping the blocks
      disconnect.detected = length(disconnect.block_sets); % change the flag
      
      %Check blocks
      disconnect.block_disconn_check = zeros(1,Nx);
      for idx = 1:length(disconnect.records)
        disconnect.blk_idx = find(disconnect.records(idx)<=block_upper_lim,1);
        if ~(disconnect.blk_idx==1)
          disconnect.block_disconn_check(disconnect.blk_idx-1) = disconnect.block_disconn_check(disconnect.blk_idx-1)+1;
        end
      end
      disconnect.block_disconn = find(disconnect.block_disconn_check);
      
      for idx = 1:length(disconnect.block_disconn)
        if param.collate_est_nz.enable_verbose %Check--
          fprintf('[analysis_max] Blocks w/ Disconnect: %d -- Records %d to %d\n',disconnect.block_disconn(idx),block_upper_lim(disconnect.block_disconn(idx)),block_upper_lim(disconnect.block_disconn(idx)+1));
        end
        BAD_from_analysis = [BAD_from_analysis disconnect.block_disconn(idx)];
      end
      
    end %if DISCONNECT detected
    
    % Save disconnect plots
    fig_idx = fig_idx +1;
    cur_fig = h_fig(fig_idx);
    clf(cur_fig);
    if param.collate_est_nz.enable_visible_plot
      figure(cur_fig); % Brings the current figure to the top
    end
    h_axes = axes('parent',cur_fig);
    
    subplot(2,1,1,h_axes)
    plot(h_axes, squeeze(rline_max_dB(wf,adc,:)));
    hold(h_axes, 'on');
    axis(h_axes,'tight');
    xlabel(h_axes,'Records');
    ylabel(h_axes,'Max Amplitude, (dBm)');
    grid(h_axes, 'on');
    if ~isempty(disconnect.records)
      plot(h_axes,disconnect.records,squeeze(rline_max_dB(wf,adc,disconnect.records)),'x');
      legend(h_axes,{'Max','Disconn'}, 'Location', 'NorthWest','FontSize',12,'FontWeight', 'bold');
    end
    hold(h_axes, 'off');
    title(h_axes,sprintf('%s : For [ %d-%d ] Max Power: Disconnect Limit = %.1f dBm',...
      param.day_seg,wf,adc,disconnect.limit), 'Interpreter', 'none');
    
    h_axes = axes('parent',cur_fig);
    subplot(2,1,2,h_axes)
    plot(h_axes,disconnect.records_mask);
    hold(h_axes,'on');
    axis(h_axes,'tight');
    ylim(h_axes,[-0.2 1.2]);
    xlabel(h_axes,'Records');
    ylabel(h_axes,sprintf('True for disconnect')); % (Block Size = %d) ,param.analysis.block_size
    
    if ~isempty(disconnect.records)
      xlim(h_axes,[disconnect.blocks(1)*param.analysis.block_size-round(param.analysis.block_size*0.2) (disconnect.blocks(end)+1)*param.analysis.block_size+round(param.analysis.block_size*0.2)]);
      if param.collate_est_nz.enable_extralines_on_plot
        for idx = 1:length(disconnect.blocks)
          try
            line(h_axes,[(disconnect.blocks(idx)+0)*param.analysis.block_size (disconnect.blocks(idx)+0)*param.analysis.block_size],[-1 2],'Color','g','LineStyle','--','LineWidth',2);
            line(h_axes,[(disconnect.blocks(idx)+1)*param.analysis.block_size (disconnect.blocks(idx)+1)*param.analysis.block_size],[-1 2],'Color','k','LineWidth',2);
          catch
            line([(disconnect.blocks(idx)+0)*param.analysis.block_size (disconnect.blocks(idx)+0)*param.analysis.block_size],[-1 2],'Color','g','LineStyle','--','LineWidth',2);
            line([(disconnect.blocks(idx)+1)*param.analysis.block_size (disconnect.blocks(idx)+1)*param.analysis.block_size],[-1 2],'Color','k','LineWidth',2);
          end
        end
        leg = {'Records Mask','Start of block','End of block'};
        legend(h_axes,leg, 'Location', 'South','FontSize',12,'FontWeight', 'bold');
      end
    end
    title(h_axes,sprintf('%s : For [ %d-%d ] Detected Disconnects = %d',param.day_seg,wf,adc,disconnect.detected), 'Interpreter', 'none');
    
    fig_fn = [ct_filename_ct_tmp(param,'',param.collate_est_nz.debug_out_dir,sprintf('%s_disconnect_wf_%d_adc_%d',param.collate_est_nz.debug_out_dir,wf,adc)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(cur_fig,fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'',param.collate_est_nz.debug_out_dir,sprintf('%s_disconnect_wf_%d_adc_%d',param.collate_est_nz.debug_out_dir,wf,adc)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(cur_fig,fig_fn);
    
    % Closing the disconnect section
    disconnect.Nregions = disconnect.detected + 1;
    if param.collate_est_nz.enable_verbose;
      fprintf('[analysis_max] Disconnects detected = %d ==> Regions in coh_noise Normalized xcorr = %d\n',disconnect.detected,disconnect.Nregions);
    end
    
    if ~param.collate_est_nz.enable_visible_plot
      try
        delete(cur_fig);
      catch
      end
    end
    clear dd rline_pow num_blocks b_size rline_max rline_max_dB temp temp2 rline_max_dB_mean rline_max_dB_std rline_max_dB_median
    % carry forward: disconnect
    
    %% Categorize blocks into sets
    
    set={};
    set_idx = 1; %current set index
    prev_set_idx = 1;
    set{set_idx} = [1]; %First record is part of the set
    
    %     xcorr_tol = 0.01; %<<<<<<<<<<<<<<<<<*####################################(LEGACY)
    xcorr_tol = param.collate_est_nz.xcorr_tol;
    
    idx=1;
    while idx<Nx %creating sets with continous blocks and transitions (for loop is not possible here)
      
      for idx2=idx+1:Nx
        if xc.norm(idx,idx2-1)-xc.norm(idx,idx2) < xcorr_tol
          if prev_set_idx==set_idx
            set{set_idx} = [set{set_idx} idx2];
          else
            set_idx=set_idx+1;
            set{set_idx} = [idx2];
            prev_set_idx = set_idx;
          end
        else
          set_idx = set_idx + 1;
          set{set_idx} = [idx2];
          break;
        end
      end
      idx = idx2;
      
    end
    
    clear prev_set_idx set_idx% add others
    
    %% Regions and Boundaries
    
    BAD_from_categories = [];
    region=struct();
    region.set={};
    region.boundary = [];
    region.bad_set={};
    blk_mask = NaN(1,Nx);
    
    for idx = 1: length(set)
      if length(set{idx}) > 4 % if the categorized set has more than 4 blocks, it is a region
        region.set = [region.set set{idx}];
        region.boundary = [region.boundary; set{idx}(1), set{idx}(end)];
        blk_mask(set{idx}(1):set{idx}(end)) = 1; % GOOD BLOCKS
      else
        region.bad_set = [region.bad_set set{idx}];
        blk_mask(set{idx}) = 0; % BAD BLOCKS %Assigning for checks
        if param.collate_est_nz.enable_verbose;
          fprintf('[categories] Blocks : %d -- Records %d to %d\n' , set{idx} , block_upper_lim(set{idx}) , block_upper_lim(set{idx}+1) ); %Bug:if 2 continuous blocks in bad_set
        end
        BAD_from_categories = [BAD_from_categories set{idx}];
      end
    end
    region.Nr = length(region.set);
    region.N_bad = length(region.bad_set);
    
    
    %% Calculate mean of regions
    %         fprintf('>> Mean xcorr in region ');
    region.id = 1:region.Nr;
    region.mean = NaN(region.Nr);
    
    for idx = 1:region.Nr
      for idx2 = 1:region.Nr
        temp = xc.norm(region.boundary(idx,1):region.boundary(idx,2),region.boundary(idx2,1):region.boundary(idx2,2));
        temp_mean = mean(temp,2);
        region.mean(idx,idx2) = mean(temp_mean);
      end
    end
    
    clear temp temp_mean
    %% Assign unique region ID based on mean
    %     fprintf('>> Assigning IDs\n');
    %     xcorr_region_tol = 0.02; %<<<<<<<<<<<<<<<<<*####################################(LEGACY)
    
    region.id_new = region.id;
    for idx=1:region.Nr
      for idx2=idx+1:region.Nr
        if abs(region.mean(idx,idx2)-region.mean(idx,idx)) < param.collate_est_nz.xcorr_region_tol %%% tol and region_mean
          region.id_new(idx2) = region.id_new(idx);
        elseif region.id_new(idx2) > region.id_new(1,idx2-1) && region.id_new(idx2)~=max(region.id_new(1:idx2-1))
          region.id_new(idx2) = max(region.id_new(1:idx2-1))+1;
        end
      end
    end
    
    region.unique_id = unique(region.id_new);
    region.unique_N = length(region.unique_id);
    
    %% Pick wf based on max(mean of xcorr in unique region)
    
    region.mask_id = zeros(Nx);
    region.mask = zeros(Nx);
    
    for idx = 1:region.Nr
      for idx2 = 1:region.Nr
        if region.id_new(idx2) == region.id_new(idx)
          region.mask_id(region.boundary(idx,1):region.boundary(idx,2),region.boundary(idx2,1):region.boundary(idx2,2)) = region.id_new(idx2);
          region.mask(region.boundary(idx,1):region.boundary(idx,2),region.boundary(idx2,1):region.boundary(idx2,2)) = 1;
        end
      end
    end
    
    %%
    region.unique_max=[];
    region.unique_max_idx=[];
    
    for idx=1:region.unique_N
      olala = region.mask_id==region.unique_id(idx);
      tmp_row_mean = sum(xc.norm.*olala ,2) ./ sum(olala ,2);
      [region.unique_max(idx), region.unique_max_idx(idx)] = max(tmp_row_mean);
      clear olala tmp_row_mean;
    end
    
    %% output structure
    layer = opsLoadLayers(param,struct('name','surface','source','layerData','layerdata_source','layerData_koenig'));
    
    coh_wf=struct();
    coh_wf.sets = region.unique_N;
    coh_wf.block_idxs = region.unique_max_idx;
    coh_wf.gps_time = simp.gps_time(coh_wf.block_idxs);
    coh_wf.elev = simp.elev(coh_wf.block_idxs);
    coh_wf.coh_noise = simp.coh_noise(:,coh_wf.block_idxs);
    
    leg={};
    for idx=1:region.unique_N
      coh_wf.twtt(idx) = interp1(layer.gps_time,layer.twtt,coh_wf.gps_time(idx));
      leg=[leg sprintf('%d: block %d %0.3f us)',region.unique_id(idx),coh_wf.block_idxs(idx),coh_wf.twtt(idx)/1e-6)]; %#ok<AGROW>
    end
    
    coh_wf.block_mask_set = zeros(1,Nx);
    coh_wf.records_mask_set = NaN(1,length(disconnect.records_mask));
    for idx=1:length(region.set)
      coh_wf.block_mask_set(region.boundary(idx,1):region.boundary(idx,2)) = region.id_new(idx);
      low_lim = block_upper_lim(region.boundary(idx,1))-param.analysis.block_size+1;
      up_lim = block_upper_lim(region.boundary(idx,2));
      if up_lim > length(disconnect.records_mask)
        up_lim = length(disconnect.records_mask);
      end
      coh_wf.records_mask_set(low_lim:up_lim) = region.id_new(idx);
    end
    
    mat_fn = ct_filename_ct_tmp(param,'',param.collate_est_nz.debug_out_dir,sprintf('%s_wf_%d_adc_%d.mat',param.collate_est_nz.debug_out_dir,wf,adc));
    fprintf('Saving %s\n', mat_fn);
    save(mat_fn, 'coh_wf', 'simp', 'xc');% , 'param',  'region', 'disconnect'
    %%
    % Save OUTPUT plots
    fig_idx = fig_idx +1;
    cur_fig = h_fig(fig_idx);
    clf(cur_fig);
    if param.collate_est_nz.enable_visible_plot
      figure(cur_fig); % Brings the current figure to the top
    end
    h_axes = axes('parent',cur_fig);
    
    subplot(4,4,[1,2,5,6],h_axes);
    imagesc(xc.norm);
    hold(h_axes,'on');
    colorbar(h_axes);
    
    for idx=1:length(region.set)
      try
        line(h_axes,[region.boundary(idx,1) region.boundary(idx,1)],[region.boundary(idx,1) region.boundary(idx,2)],'Color','g','LineWidth',2); %Low (vertical)
        line(h_axes,[region.boundary(idx,2) region.boundary(idx,2)],[region.boundary(idx,1) region.boundary(idx,2)],'Color','g','LineWidth',2); %Up (vertical)
        line(h_axes,[region.boundary(idx,1) region.boundary(idx,2)],[region.boundary(idx,1) region.boundary(idx,1)],'Color','g','LineWidth',2); %Low (horizontal)
        line(h_axes,[region.boundary(idx,1) region.boundary(idx,2)],[region.boundary(idx,2) region.boundary(idx,2)],'Color','g','LineWidth',2); %Up (horizontal)
      catch
        line([region.boundary(idx,1) region.boundary(idx,1)],[region.boundary(idx,1) region.boundary(idx,2)],'Color','g','LineWidth',2); %Low (vertical)
        line([region.boundary(idx,2) region.boundary(idx,2)],[region.boundary(idx,1) region.boundary(idx,2)],'Color','g','LineWidth',2); %Up (vertical)
        line([region.boundary(idx,1) region.boundary(idx,2)],[region.boundary(idx,1) region.boundary(idx,1)],'Color','g','LineWidth',2); %Low (horizontal)
        line([region.boundary(idx,1) region.boundary(idx,2)],[region.boundary(idx,2) region.boundary(idx,2)],'Color','g','LineWidth',2); %Up (horizontal)
      end
    end
    for idx=1:region.unique_N
      try
        line(h_axes,[1 Nx],[region.unique_max_idx(idx) region.unique_max_idx(idx)],'Color','w','LineWidth',2,'LineStyle',':'); %Mean Line
      catch
        line([1 Nx],[region.unique_max_idx(idx) region.unique_max_idx(idx)],'Color','w','LineWidth',2,'LineStyle',':'); %Mean Line
      end
    end
    
    if param.collate_est_nz.enable_extralines_on_plot
      for idx=1:length(region.bad_set)
        try
          line(h_axes,[region.bad_set{idx} region.bad_set{idx}],[1 Nx],'Color','r','LineWidth',2); %BAD from sets
        catch
          line([region.bad_set{idx} region.bad_set{idx}],[1 Nx],'Color','r','LineWidth',2); %BAD from sets
        end
      end
      for idx=1:length(BAD_from_analysis)
        try
          line(h_axes,[BAD_from_analysis(idx) BAD_from_analysis(idx)],[1 Nx],'Color','k','LineWidth',3,'LineStyle','--'); %BAD because disconnect
        catch
          line([BAD_from_analysis(idx) BAD_from_analysis(idx)],[1 Nx],'Color','k','LineWidth',3,'LineStyle','--'); %BAD because disconnect
        end
      end
    end
    
    grid(h_axes,'on');
    xlabel(h_axes,sprintf('Nx (%d) -->',Nx));
    ylabel(h_axes,sprintf('Nx (%d) -->',Nx));
    title(h_axes,sprintf('Norm xcorr (%.2f,%.2f)',param.collate_est_nz.xcorr_tol,param.collate_est_nz.xcorr_region_tol));
    hold(h_axes,'off');
    
    h_axes = axes('parent',cur_fig);
    subplot(4,4,[9,10],h_axes);
    plot(h_axes,coh_wf.block_mask_set,'LineWidth',2);
    hold(h_axes,'on');
    xlabel(h_axes,sprintf('Blocks Nx (%d) -->',Nx));
    ylabel(h_axes,'Set ID');
    grid(h_axes,'on');
    axis(h_axes,'tight');
    ylim(h_axes,[0 max(region.id_new)+1]);
%     title(h_axes,'Blocks classified into sets');
    plot(h_axes,coh_wf.block_idxs,coh_wf.block_mask_set(coh_wf.block_idxs),'x','LineWidth',4);
    hold(h_axes,'off');
    
    h_axes = axes('parent',cur_fig);
    subplot(4,4,[13,14],h_axes);
    plot(h_axes,coh_wf.records_mask_set,'LineWidth',2);
    hold(h_axes,'on');
    xlabel(h_axes,sprintf('Records (%d) -->',length(coh_wf.records_mask_set)));
    ylabel(h_axes,'Set ID');
    grid(h_axes,'on');
    axis(h_axes,'tight');
    ylim(h_axes,[0 max(region.id_new)+1]);
%     title(h_axes,'Blocks classified into sets');
%     plot(h_axes,coh_wf.block_idxs,coh_wf.block_mask_set(coh_wf.block_idxs),'x','LineWidth',4);
    hold(h_axes,'off');
    
    h_axes = axes('parent',cur_fig);
    subplot(4,4,[3,4,7,8,11,12,15,16],h_axes);
    plot(h_axes,lp(coh_wf.coh_noise));
    grid(h_axes,'on');
    axis(h_axes,'tight');
    ylim(h_axes,[-100 50]);
    legend(h_axes,leg);
    title(h_axes,'Coh noise for block in each set');
    
    try
      sgtitle(cur_fig,sprintf('%s : For [ %d-%d ] collate_est_nz',param.day_seg,wf,adc), 'Interpreter', 'none');
    catch
      suptitle(sprintf('%s : For [ %d-%d ] collate est nz',strrep(param.day_seg,'_','-'),wf,adc));
    end
    
    fig_fn = [ct_filename_ct_tmp(param,'',param.collate_est_nz.debug_out_dir,sprintf('%s_wf_%d_adc_%d',param.collate_est_nz.debug_out_dir,wf,adc)) '.fig'];
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(cur_fig,fig_fn);
    fig_fn = [ct_filename_ct_tmp(param,'',param.collate_est_nz.debug_out_dir,sprintf('%s_wf_%d_adc_%d',param.collate_est_nz.debug_out_dir,wf,adc)) '.jpg'];
    fprintf('Saving %s\n', fig_fn);
    ct_saveas(cur_fig,fig_fn);
    
    if ~param.collate_est_nz.enable_visible_plot
      try
        delete(cur_fig);
      catch
      end
    end
    
  end
end

end