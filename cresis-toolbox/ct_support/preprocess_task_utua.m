function success = preprocess_task_utua(param)
% success = preprocess_task_utua(param)
%
% Support function for preprocess.m
%
% Example:
% Called from preprocess_task.m
%
% Author: John Paden
%
% See also: run_preprocess.m, preprocess.m, preprocess_task.m,
% preprocess_task_arena.m, preprocess_task_cresis.m

% 1. Use what we already have
% 2. Use file version instead of radar name
% 3. Use last read to allow header wraps... or not.

%% Input checks
% =========================================================================

if ~isfield(param.config,'field_time_gap') || isempty(param.config.field_time_gap)
  param.config.field_time_gap = 'utc_time_sod';
end

if ~isfield(param.config,'plots_visible') || isempty(param.config.plots_visible)
  param.config.plots_visible = 1;
end

%% Read Headers
% =========================================================================

% Save concatenated temporary file
fn_board_hdrs = ct_filename_ct_tmp(param,'','headers', fullfile(param.config.date_str,'board_hdrs.mat'));
num_board_to_load = numel(param.config.board_map);
board_hdrs = cell(1,num_board_to_load);
failed_load = cell(1,num_board_to_load);
fns_list = cell(1,num_board_to_load);
if exist(fn_board_hdrs,'file')
  try
    fprintf('Found %s\n  Trying to load...\n', fn_board_hdrs);
    load(fn_board_hdrs,'board_hdrs','fns_list','failed_load');
    num_board_to_load = 0;
  catch ME
    ME.getReport
  end
end

for board_idx = 1:num_board_to_load
  %% Read Headers: Filenames
  board = param.config.board_map{board_idx};
  board_folder_name = param.config.board_folder_name;
  board_folder_name = regexprep(board_folder_name,'%b',board);
  
  get_filenames_param = struct('regexp',param.config.file.regexp);
  fns = get_filenames(fullfile(param.config.base_dir,board_folder_name), ...
    param.config.file.prefix, param.config.file.midfix, ...
    param.config.file.suffix, get_filenames_param);
  
  if param.config.online_mode == 2
    fns = fns(end);
  end
  fns_list{board_idx} = fns;
  
  if isempty(fns)
    error('No files found matching %s*%s*%s', ...
      fullfile(param.config.base_dir,board_folder_name,param.config.file.prefix), ...
      param.config.file.midfix, param.config.file.suffix);
  end
  
  % Assumption is that fns is in chronological order. Most radar systems
  % have filenames that are in chronological order with a simple
  % alphabetical sort.
  
  %% Read Headers: File Loop
  failed_load{board_idx} = false(size(fns));
  board_hdrs{board_idx}.unknown = [];
  board_hdrs{board_idx}.radar_time = [];
  board_hdrs{board_idx}.radar_time_1pps = [];
  board_hdrs{board_idx}.epri = [];
  board_hdrs{board_idx}.seconds = [];
  board_hdrs{board_idx}.fraction = [];
  board_hdrs{board_idx}.waveform_ID = [];
  board_hdrs{board_idx}.counter = [];
  board_hdrs{board_idx}.file_idxs = [];
  for fn_idx = 1:length(fns)
    
    % Create temporary filename that will store the header information for
    % this file.
    fn = fns{fn_idx};
    if any(param.config.file.version == [405 406])
      [~,fn_name,ext] = fileparts(fn);
      fn_name = [fn_name,ext];
    else
      [~,fn_name] = fileparts(fn);
    end
    tmp_hdr_fn = ct_filename_ct_tmp(param,'','headers', ...
      fullfile(board_folder_name, [fn_name '.mat']));
    tmp_hdr_fn_dir = fileparts(tmp_hdr_fn);
    if ~exist(tmp_hdr_fn_dir,'dir')
      mkdir(tmp_hdr_fn_dir);
    end
    
    keyboard
    block = utua_rds.tdms2block_loop(fn,tmp_hdr_fn);
    
    if any(param.config.file.version == [1])
      [file_size offset epri seconds fraction] = basic_load_hdr_mex(fn,hdr_param.frame_sync,hdr_param.field_offsets,hdr_param.field_types,hdr_param.file_mode);
      
      % Find bad records by checking their size (i.e. the distance between
      % frame syncs which should be constant).
      expected_rec_size = median(diff(offset));
      meas_rec_size = diff(offset);
      bad_mask = meas_rec_size ~= expected_rec_size;
      bad_mask(end+1) = file_size < offset(end) + expected_rec_size;
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      epri = double(epri(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','wfs');
      
      
      % Remove bad records (i.e. ones with sizes that are not expected
      offset = double(offset(~bad_mask));
      epri = double(epri(~bad_mask));
      seconds = double(seconds(~bad_mask));
      fraction = double(fraction(~bad_mask));
      counter = double(counter(~bad_mask));
      
      save(tmp_hdr_fn,'offset','epri','seconds','fraction','counter','wfs');
      
    end
    
    % Load and concatenate temporary file
    if any(param.config.file.version == [101])
      hdr = load(tmp_hdr_fn);
      board_hdrs{board_idx}.unknown ...
        = cat(2,board_hdrs{board_idx}.unknown,reshape(hdr.unknown,[1 length(hdr.unknown)]));
      board_hdrs{board_idx}.seconds ...
        = cat(2,board_hdrs{board_idx}.seconds,reshape(hdr.seconds,[1 length(hdr.seconds)]));
      board_hdrs{board_idx}.fraction ...
        = cat(2,board_hdrs{board_idx}.fraction,reshape(hdr.fraction,[1 length(hdr.fraction)]));
    end
    board_hdrs{board_idx}.file_idxs = cat(2,board_hdrs{board_idx}.file_idxs,fn_idx*ones([1 length(hdr.offset)]));
    
  end
end

% Save concatenated temporary file
fn_board_hdrs_dir = fileparts(fn_board_hdrs);
if ~exist(fn_board_hdrs_dir,'dir')
  mkdir(fn_board_hdrs_dir);
end
save(fn_board_hdrs,'-v7.3','board_hdrs','fns_list','failed_load');

%% List bad files
% =========================================================================
for board_idx = 1:numel(param.config.board_map)
  if any(failed_load{board_idx})
    warning('Some files failed to load, consider deleting these to avoid problems.');
    for fn_idx = find(failed_load{board_idx})
      fprintf('  %s\n', fns_list{board_idx}{fn_idx});
    end
  end
end
  
%% Correct time variable
% =========================================================================
for board_idx = 1:numel(param.config.board_map)
  
  if any(param.config.file.version == [1 2 3 4 5 6 7 8 101 403 404 407 408])
    utc_time_sod = double(board_hdrs{board_idx}.seconds) + double(board_hdrs{board_idx}.fraction) / param.config.cresis.clk;
    epri = double(board_hdrs{board_idx}.epri);
    
    if 0
      % Test sequences
      utc_time_sod = [0 1 2 3 10000 5 6 7 8 9 10 11 12 13 24 25 26 27 28 19 20 21 22 19 20 21 22]
      utc_time_sod = utc_time_sod + 0.0001*randn(size(utc_time_sod))
      epri = 100 + [1:23, 20:23]
      epri(15) = 5000;
    end
    
    % Estimate the pulse repetition interval, PRI
    PRI = medfilt1(diff(utc_time_sod),51);
    
    % Create an EPRI sequence from the time record
    time_epri = utc_time_sod(1)/PRI(1) + [0 cumsum(diff(utc_time_sod) ./ PRI)];
    [~,good_time_idx] = min(abs(utc_time_sod - median(utc_time_sod)));
    time_epri = time_epri - time_epri(good_time_idx);
    
    % Find the difference of the time-generated epri and the recorded epri
    dtime_epri = diff(time_epri);
    depri = diff(epri);
    
    % Find good/bad differences. Mask values are:
    %  0: both differences are bad
    %  1: EPRI good
    %  2: Time-generated EPRI good
    %  3: EPRI and time-generated EPRI good
    dtime_epri_threshold = 0.1; % Allow for 10% PRI error
    mask = (depri == 1) + (2*(abs(dtime_epri-1) < dtime_epri_threshold));
    % If the EPRI's both indicate the same number of skipped records,
    % consider it a good difference.
    mask(mask ~= 3 & depri == round(dtime_epri)) = 3;
    
    % Fix differenced time-generated EPRIs using differenced EPRIs
    dtime_epri(mask==1) = depri(mask==1);
    % Fix differenced EPRIs using differenced time-generated EPRIs
    depri(mask==2) = round(dtime_epri(mask==2));
    
    % Find sequences of good records (where mask > 0) and deal with each
    % segment separately.
    good_out_mask = false(size(utc_time_sod));
    start_idx = find(mask ~= 0,1);
    while ~isempty(start_idx)
      stop_idx = start_idx-1 + find(mask(start_idx+1:end) == 0,1);
      if isempty(stop_idx)
        stop_idx = numel(mask);
      end
      
      % Find a median point in each segment and assume this value is good
      [~,good_time_idx] = min(abs(utc_time_sod(start_idx:stop_idx+1) - median(utc_time_sod(start_idx:stop_idx+1))));
      [~,good_epri_idx] = min(abs(epri(start_idx:stop_idx+1) - median(epri(start_idx:stop_idx+1))));
      
      % Reconstruct epri
      tmp = [0 cumsum(depri(start_idx:stop_idx))];
      tmp = tmp - tmp(good_epri_idx) + epri(start_idx-1+good_epri_idx);
      
      % Reconstruct time from time-generated EPRIs
      tmp = [0 cumsum(dtime_epri(start_idx:stop_idx).*PRI(start_idx:stop_idx))];
      tmp = tmp - tmp(good_time_idx) + utc_time_sod(start_idx-1+good_time_idx);
      utc_time_sod_new(start_idx:stop_idx+1) = tmp;
      
      % Mark these records as good outputs
      good_out_mask(start_idx:stop_idx+1) = true;
      
      % Find the next sequence
      start_idx = stop_idx + find(mask(stop_idx+1:end) ~= 0,1);
    end
    
    % START IMPORTANT TIME CORRECTION
    
    % NOTE: If custom per-day utc_time_sod corrections are required,
    % include that code here to change the utc_time_sod_new.m. This
    % function should be able to identify the radar and date of the data
    % and apply the custom correction based on that.
    
    % utc_time_sod_new = preprocess_task_utua_custom(); % <-- CREATE
    
    % END IMPORTANT TIME CORRECTION
    
    h_fig = get_figures(3,param.config.plots_visible,mfilename);
    
    clf(h_fig(1)); h_axes = axes('parent',h_fig(1));
    plot(h_axes,utc_time_sod);
    hold(h_axes,'on');
    plot(h_axes,utc_time_sod_new,'r');
    hold(h_axes,'off');
    xlabel(h_axes,'Record number');
    ylabel(h_axes,'UTC Time SOD (sec)');
    legend(h_axes,'Original','Corrected','location','best');
    title(h_axes,sprintf('%s: UTC time original and corrected should\nmatch except at outliers',param.config.date_str),'fontsize',10);
    
    UTC_MAX_ERROR = 0.1;
    clf(h_fig(2)); h_axes = axes('parent',h_fig(2));
    plot(h_axes,utc_time_sod - utc_time_sod_new);
    xlabel(h_axes,'Record number');
    ylabel(h_axes,'Time correction (sec)');
    ylim(h_axes,[-UTC_MAX_ERROR UTC_MAX_ERROR]);
    title(h_axes,sprintf('%s: Time correction should be within limits\nexcept for a few outliers.',param.config.date_str),'fontsize',10);
    
    clf(h_fig(3)); h_axes = axes('parent',h_fig(3));
    h_axes = subplot(2,1,1);
    plot(h_axes,diff(epri),'.');
    ylabel(h_axes,'Diff EPRI');
    h_axes = subplot(2,1,2);
    plot(h_axes,diff(epri),'.');
    ylim(h_axes,[-3 5]);
    xlabel(h_axes,'Record number');
    ylabel(h_axes,'Diff EPRI');
    title(h_axes,sprintf('%s: Should be 1 except occasional record drops.',param.config.date_str),'fontsize',10);
    
    fn_fig = ct_filename_ct_tmp(param,'','headers', sprintf('preprocess_%s_utc_time_sod_board_%d.fig',param.config.date_str,board_idx));
    fn_fig_dir = fileparts(fn_fig);
    if ~exist(fn_fig_dir,'dir')
      mkdir(fn_fig_dir);
    end
    fprintf('Saving %s\n', fn_fig);
    saveas(h_fig(1),fn_fig);
    fn_fig = [fn_fig(1:end-3), 'jpg'];
    fprintf('Saving %s\n', fn_fig);
    saveas(h_fig(1),fn_fig);
    fn_fig = ct_filename_ct_tmp(param,'','headers', sprintf('preprocess_%s_utc_time_sod_error_board_%d.fig',param.config.date_str,board_idx));
    fprintf('Saving %s\n', fn_fig);
    saveas(h_fig(2),fn_fig);
    fn_fig = [fn_fig(1:end-3), 'jpg'];
    fprintf('Saving %s\n', fn_fig);
    saveas(h_fig(2),fn_fig);
    fn_fig = ct_filename_ct_tmp(param,'','headers', sprintf('preprocess_%s_epri_board_%d.fig',param.config.date_str,board_idx));
    fprintf('Saving %s\n', fn_fig);
    saveas(h_fig(3),fn_fig);
    fn_fig = [fn_fig(1:end-3), 'jpg'];
    fprintf('Saving %s\n', fn_fig);
    saveas(h_fig(3),fn_fig);
    
    utc_time_sod = utc_time_sod_new;
    
    % Check for day wraps in the UTC time seconds of day
    day_wrap_idxs = find(diff(utc_time_sod) < -50000);
    board_hdrs{board_idx}.day_wrap_offset = zeros(size(utc_time_sod));
    if ~isempty(day_wrap_idxs)
      fprintf('Found %d potential day wraps in board %d. Unwrapping times.', numel(day_wrap_idxs), board_idx);
      board_hdrs{board_idx}.day_wrap_offset = zeros(size(utc_time_sod));
      for day_wrap_idx = day_wrap_idxs
        board_hdrs{board_idx}.day_wrap_offset(day_wrap_idx+1:end) = board_hdrs{board_idx}.day_wrap_offset(day_wrap_idx+1:end) + 86400;
      end
      utc_time_sod = utc_time_sod + board_hdrs{board_idx}.day_wrap_offset;
    end
    
    clear depri dtime_epri dtime_epri_threshold time_epri utc_time_sod_new day_wrap_idxs;
    board_hdrs{board_idx}.utc_time_sod = utc_time_sod;
    board_hdrs{board_idx}.epri = epri;
    
  elseif any(param.config.file.version == [405 406])
    utc_time_sod = board_hdrs{board_idx}.seconds; % this is actually comp_time but doesn't need to
    % be converted to actual utc_time_sod since it's only looking at gaps in
    % the data
    
    day_wrap_idxs = find(diff(utc_time_sod) < -50000);
    day_wrap_offset = zeros(size(utc_time_sod));
    for day_wrap_idx = day_wrap_idxs
      day_wrap_offset(day_wrap_idx+1:end) = day_wrap_offset(day_wrap_idx+1:end) + 86400;
    end
    utc_time_sod = utc_time_sod + day_wrap_offset;
    
    time_gaps = find(abs(diff(utc_time_sod)) > MAX_TIME_GAP);
    
    reason = {};
    
    figure(1); clf;
    plot(utc_time_sod);
    ylabel('UTC time seconds of day');
    xlabel('Record');
    grid on;
    hold on;
    plot(time_gaps, utc_time_sod(time_gaps),'ro');
    hold off;
    
    names = fieldnames(hdr_log(1));
    if param.config.file.version == 406
      change_fields = [2 5 6 7 8 9 10 11 12 17 18 19 20 21];
    elseif param.config.file.version == 405
      change_fields = [2 5 6 7 8 9 10 11 12];
    end
    
    % Detect changes in header file information that would necessitate
    % creating a new segment.
    hdr_gaps = [];
    change_log = {};
    for n_idx = change_fields
      %     fprintf('Checking %s...\n',names{n_idx})
      hdr_changes = eval(sprintf('find(diff(squeeze(horzcat(hdr_log(:).%s))) ~= 0)',names{n_idx}));
      hdr_gaps = [hdr_gaps, hdr_changes];
      if ~isempty(hdr_changes)
        reason = [reason; repmat({names{n_idx}},eval(sprintf('length(find(diff(squeeze(horzcat(hdr_log(:).%s))) ~= 0))',names{n_idx})),1)];
        for cl_idx=1:length(hdr_changes)
          change_log{end+1} = eval(sprintf('[hdr_log(hdr_changes(cl_idx)).%s hdr_log(hdr_changes(cl_idx)+1).%s]',names{n_idx},names{n_idx}));
        end
      end
    end
    
    [hdr_gaps I J] = unique(hdr_gaps);
    comb_reason = cell(1,length(I));
    comb_change = cell(1,length(I));
    for uidx=1:length(J)
      comb_reason{J(uidx)} = [comb_reason{J(uidx)} reason{uidx}];
      comb_change{J(uidx)} = [comb_change{J(uidx)} change_log{uidx}];
    end
    
    % Convert from file index to file offset
    hdr_idxs = [];
    time_idxs = [];
    for gaps_idx = 1:length(hdr_gaps)
      hdr_idxs(gaps_idx) = find(utc_time_sod == htime(hdr_gaps(gaps_idx)+1),1,'last')-1;
    end
    
    if ~isempty(hdr_gaps)
      figure(1);
      hold on;
      plot(hdr_idxs, utc_time_sod(hdr_idxs),'go');
      hold off;
    end
    
    [time_gaps I J] = unique([time_gaps hdr_idxs]);
    if ~isempty(reason)
      extend_unique = [repmat({'time'},length(time_gaps),1); comb_reason.'];
      reason_unique = extend_unique(I);
      extend_change_log = [repmat({'N/A'},length(time_gaps),1); comb_change.'];
      change_log_unique = extend_change_log(I);
    else
      reason_unique = repmat({''},length(I),1);
      change_log_unique = repmat({''},length(I),1);
    end
  else
    error('Not supported');
  end

  board_hdrs{board_idx}.utc_time_sod = utc_time_sod;
  board_hdrs{board_idx}.epri = epri;
end

% if any(strcmpi(radar_name,{'accum'}))
%   %% Accum filter of bad records
%   % Repeated data often occurs in accum during a FIFO overflow. Remove
%   % these blocks.
%   bad_mask = zeros(size(utc_time-sod));
%   start_segment_time = utc_time_sod(1);
%   cur_time = utc_time_sod(1);
%   for idx = 2:utc_time_sod
%     diff_time = utc_time_sod(idx) - utc_time_sod(idx-1);
%     if diff_time > MAX_TIME_GAP
%       start_segment_time = utc_time_sod(idx);
%       cur_time = utc_time_sod(idx);
%     else
%       if utc_time_sod(idx) < start_segment_time
%       end
%       if utc_time_sod(idx) < cur_time
%         bad_mask(idx) = 1;
%       end
%     end
%   end
% end

%% Create Segments
% =========================================================================

if any(param.config.file.version == [403 404 407 408])
  %% Create Segments: Read XML settings
  % NI XML settings files available, break segments based on settings files
  % and header information
  
  xml_version = param.config.daq.xml_version;
  cresis_xml_mapping;
  
  settings_fn_dir = fullfile(param.config.base_dir,param.config.config_folder_name);
  fprintf('\nSettings Directory: %s\n\n', settings_fn_dir);
  
  % Read XML files in this directory
  [settings,settings_enc] = read_ni_xml_directory(settings_fn_dir,xml_file_prefix,false);
  
  % Get the date information out of the filename
  fn_datenums = {};
  for board_idx = 1:numel(param.config.board_map)
    fn_datenums{board_idx} = [];
    for data_fn_idx = 1:length(fns_list{board_idx})
      fname = fname_info_mcords2(fns_list{board_idx}{data_fn_idx});
      fn_datenums{board_idx}(end+1) = fname.datenum;
    end
  end
  
  %% Create Segments: Print settings
  oparams = {};
  for set_idx = 1:length(settings)
    % Print out settings
    [~,settings_fn_name] = fileparts(settings(set_idx).fn);
    fprintf('===================== Setting %d =================\n', set_idx);
    fprintf('%s: %d waveforms\n', settings_fn_name, length(settings(set_idx).(config_var).Waveforms));
    if isfield(settings(set_idx),'XML_File_Path')
      fprintf('  %s\n', settings(set_idx).XML_File_Path{1}.values{1});
    end
    fprintf('   PRF:'); fprintf(' %g', settings(set_idx).(config_var).(prf_var)); fprintf('\n');
    fprintf('   Amp:'); fprintf(' %g', settings(set_idx).(config_var).(ram_amp_var)); fprintf('\n');
    fprintf('   Tukey:'); fprintf(' %g', settings(set_idx).(config_var).RAM_Taper); fprintf('\n');
    Tpd = double(settings(set_idx).(config_var).Waveforms(1).Len_Mult)*settings(set_idx).(config_var).Base_Len;
    fprintf('   f0-f1:'); fprintf(' %g-%g MHz %g us', settings(set_idx).(config_var).Waveforms(1).Start_Freq(1)/1e6, ...
      settings(set_idx).(config_var).Waveforms(1).Stop_Freq(1)/1e6, Tpd*1e6); fprintf('\n');
    fprintf('   Tx Mask:'); fprintf(' %g', settings(set_idx).(config_var).Waveforms(1).TX_Mask); fprintf('\n');
    for wf = 1:length(settings(set_idx).(config_var).Waveforms)
      fprintf('    WF %d Atten:', wf); fprintf(' %g', settings(set_idx).(config_var).Waveforms(wf).Attenuator_2); fprintf('\n');
      fprintf('    WF %d Len:', wf); fprintf(' %.1f us', 1e6*settings(set_idx).(config_var).Base_Len*settings(set_idx).(config_var).Waveforms(wf).Len_Mult); fprintf('\n');
    end
    
    for board_idx = 1:numel(param.config.board_map)
      if set_idx < length(settings)
        settings(set_idx).file_matches{board_idx} = find(fn_datenums{board_idx} >= settings(set_idx).datenum & fn_datenums{board_idx} < settings(set_idx+1).datenum);
      else
        settings(set_idx).file_matches{board_idx} = find(fn_datenums{board_idx} >= settings(set_idx).datenum);
      end
    end
    
    % Associate default parameters with each settings
    default = default_radar_params_settings_match(param.config.defaults,settings(set_idx));
    oparams{end+1} = default;
    oparams{end} = rmfield(oparams{end},'config_regexp');
    oparams{end} = rmfield(oparams{end},'name');
    
    % Parameter spreadsheet
    % =======================================================================
    oparams{end}.cmd.notes = default.name;
    
    oparams{end}.records.file.base_dir = param.config.base_dir;
    oparams{end}.records.file.board_folder_name = param.config.board_folder_name;
    if ~isempty(oparams{end}.records.file.board_folder_name) ...
        && oparams{end}.records.file.board_folder_name(1) ~= filesep
      % Ensures that board_folder_name is not a text number which Excel
      % will misinterpret as a numeric type
      oparams{end}.records.file.board_folder_name = ['/' oparams{end}.records.file.board_folder_name];
    end
    oparams{end}.records.file.boards = param.config.board_map;
    oparams{end}.records.file.version = param.config.file.version;
    oparams{end}.records.file.prefix = param.config.file.prefix;
    oparams{end}.records.file.clk = param.config.cresis.clk;
    oparams{end}.radar.prf = settings(set_idx).(config_var).(prf_var);

    % Usually the default.radar.wfs structure only has one waveform
    % entry which is to be copied to all the waveforms.
    if numel(oparams{end}.radar.wfs) == 1
      oparams{end}.radar.wfs = repmat(oparams{end}.radar.wfs,[1 numel(settings(set_idx).(config_var).Waveforms)]);
    end
      
    for wf = 1:numel(settings(set_idx).(config_var).Waveforms)
      oparams{end}.radar.wfs(wf).Tpd = double(settings(set_idx).(config_var).Waveforms(wf).Len_Mult)*settings(set_idx).(config_var).Base_Len;
      oparams{end}.radar.wfs(wf).f0 = settings(set_idx).(config_var).Waveforms(wf).Start_Freq(1);
      oparams{end}.radar.wfs(wf).f1 = settings(set_idx).(config_var).Waveforms(wf).Stop_Freq(1);
      oparams{end}.radar.wfs(wf).tukey = settings(set_idx).(config_var).RAM_Taper;
      % Transmit weights
      if any(param.config.file.version == [403 407 408])
        tx_mask_inv = fliplr(~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0'));
        tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv ./ param.config.max_tx.*param.config.max_tx_voltage;
      else
        tx_mask_inv = ~(dec2bin(double(settings(set_idx).(config_var).Waveforms(wf).TX_Mask),8) - '0');
        tx_weights = double(settings(set_idx).(config_var).(ram_var)) .* tx_mask_inv ./ param.config.max_tx.*param.config.max_tx_voltage;
      end

      tx_weights = tx_weights(logical(param.config.tx_enable));
      oparams{end}.radar.wfs(wf).tx_weights = tx_weights;
      
      % ADC Gains
      atten = double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_1(1)) ...
        + double(settings(set_idx).(config_var).Waveforms(wf).Attenuator_2(1));
      oparams{end}.radar.wfs(wf).adc_gains_dB = param.config.cresis.rx_gain_dB - atten(1)*ones(1,length(oparams{end}.radar.wfs(wf).rx_paths));
      
      % DDC mode and frequency
      if isfield(settings(set_idx), 'DDC_Ctrl')
        oparams{end}.radar.wfs(wf).DDC_dec = 2^(2+settings(set_idx).DDC_Ctrl.DDC_sel.Val);
        oparams{end}.radar.wfs(wf).DDC_freq = settings(set_idx).DDC_Ctrl.(NCO_freq)*1e6;
      else
        oparams{end}.radar.wfs(wf).DDC_dec = 1;
        oparams{end}.radar.wfs(wf).DDC_freq = 0;
      end
    end
    
    counters = {};
    file_idxs = {};
    for board_idx = 1:numel(param.config.board_map)
      counters{board_idx} = double(board_hdrs{board_idx}.(param.config.field_time_gap));
      file_idxs{board_idx} = board_hdrs{board_idx}.file_idxs;
      day_wrap_offset{board_idx} = board_hdrs{board_idx}.day_wrap_offset;
      % Restrict to just these settings
      if isempty(settings(set_idx).file_matches{board_idx})
        counters{board_idx} = [];
        file_idxs{board_idx} = [];
      else
        mask = file_idxs{board_idx} >= settings(set_idx).file_matches{board_idx}(1) ...
          & file_idxs{board_idx} <= settings(set_idx).file_matches{board_idx}(end);
        counters{board_idx} = counters{board_idx}(mask);
        file_idxs{board_idx} = file_idxs{board_idx}(mask);
        day_wrap_offset{board_idx} = day_wrap_offset{board_idx}(mask);
      end
    end
    [segs,stats] = create_segments(counters,file_idxs,day_wrap_offset,param.config.max_time_gap);
    
    if 1
      % Debug: Test Code
      for seg_idx = 1:length(segs)
        fprintf('Segment %d\n', seg_idx);
        disp(segs(seg_idx))
      end
      
      fprintf('On time: %g\n', sum(stats.on_time));
      fprintf('Seg\tOn%%\tOn');
      for board_idx = 1:size(stats.board_time,2)
        fprintf('\t%d%%\t%d', board_idx, board_idx);
      end
      fprintf('\n');
      
      for seg_idx = 1:length(segs)
        fprintf('%d\t%.0f%%\t%.1g', seg_idx, stats.on_time(seg_idx)/sum(stats.on_time)*100, stats.on_time(seg_idx));
        for board_idx = 1:size(stats.board_time,2)
          fprintf('\t%.0f%%\t%.1g', stats.board_time(seg_idx,board_idx)/stats.on_time(seg_idx)*100, stats.board_time(seg_idx,board_idx));
        end
        fprintf('\n');
      end
    end
    
    if isempty(segs)
      oparams = oparams(1:end-1);
    else
      for seg_idx = 1:length(segs)
        segment = segs(seg_idx);
        if seg_idx > 1
          oparams{end+1} = oparams{end};
        end
        oparams{end}.day_seg = sprintf('%s_%02d',param.config.date_str,length(oparams));
        oparams{end}.records.file.start_idx = segment.start_idxs;
        oparams{end}.records.file.stop_idx = segment.stop_idxs;
        oparams{end}.records.gps.time_offset = default.records.gps.time_offset + segment.day_wrap_offset;
      end
    end
  end
  
elseif any(param.config.file.version == [410])
  % MCRDS headers available, break segments based on filenames
  
else
  % No heading information, break segments based on time, epri, or radar
  % counter information (param.config.field_time_gap and
  % param.config.max_time_gap determine which field and gap size to use).
  counters = {};
  file_idxs = {};
  for board_idx = 1:numel(param.config.board_map)
    counters{board_idx} = double(board_hdrs{board_idx}.(param.config.field_time_gap));
    file_idxs{board_idx} = board_hdrs{board_idx}.file_idxs;
    day_wrap_offset{board_idx} = board_hdrs{board_idx}.day_wrap_offset;
  end
  [segs,stats] = create_segments(counters,file_idxs,day_wrap_offset,param.config.max_time_gap);
  
  if 1
    % Debug: Test Code
    for seg_idx = 1:length(segs)
      fprintf('Segment %d\n', seg_idx);
      disp(segs(seg_idx))
    end
    
    fprintf('On time: %g\n', sum(stats.on_time));
    fprintf('Seg\tOn%%\tOn');
    for board_idx = 1:size(stats.board_time,2)
      fprintf('\t%d%%\t%d', board_idx, board_idx);
    end
    fprintf('\n');
    
    for seg_idx = 1:length(segs)
      fprintf('%d\t%.0f%%\t%.1g', seg_idx, stats.on_time(seg_idx)/sum(stats.on_time)*100, stats.on_time(seg_idx));
      for board_idx = 1:size(stats.board_time,2)
        fprintf('\t%.0f%%\t%.1g', stats.board_time(seg_idx,board_idx)/stats.on_time(seg_idx)*100, stats.board_time(seg_idx,board_idx));
      end
      fprintf('\n');
    end
  end

  % Create the parameters to output
  oparams = {};
  for segment_idx = 1:length(segs)
    segment = segs(segment_idx);
    
    % Determine which default parameters to use
    % =======================================================================
    match_idx = 1;
    oparams{end+1} = param.config.defaults{match_idx};
    oparams{end} = rmfield(oparams{end},'config_regexp');
    oparams{end} = rmfield(oparams{end},'name');
    
    % Parameter spreadsheet
    % =======================================================================
    oparams{end}.day_seg = sprintf('%s_%02d',param.config.date_str,segment_idx);
    oparams{end}.cmd.notes = param.config.defaults{match_idx}.name;
    
    oparams{end}.records.file.start_idx = segment.start_idxs;
    oparams{end}.records.file.stop_idx = segment.stop_idxs;
    oparams{end}.records.gps.time_offset = param.config.defaults{match_idx}.records.gps.time_offset + segment.day_wrap_offset;
    
    oparams{end}.records.file.base_dir = param.config.base_dir;
    oparams{end}.records.file.board_folder_name = param.config.board_folder_name;
    if ~isempty(oparams{end}.records.file.board_folder_name) ...
        && oparams{end}.records.file.board_folder_name(1) ~= filesep
      % Ensures that board_folder_name is not a text number which Excel
      % will misinterpret as a numeric type
      oparams{end}.records.file.board_folder_name = ['/' oparams{end}.records.file.board_folder_name];
    end
    if ~isnan(str2double(oparams{end}.records.file.board_folder_name))
      oparams{end}.records.file.board_folder_name = ['/' oparams{end}.records.file.board_folder_name];
    end
    oparams{end}.records.file.boards = param.config.board_map;
    oparams{end}.records.file.version = param.config.file.version;
    oparams{end}.records.file.prefix = param.config.file.prefix;
    oparams{end}.records.file.clk = param.config.cresis.clk;
  end
  
end

%% Print out segments
% =========================================================================
if ~isempty(param.config.param_fn)
  % Print parameter spreadsheet values to stdout and param_txt_fn
  % =========================================================================
  fid = 1;
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  cmd\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'cmd',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  records\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'records',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  qlook\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'qlook',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  sar\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'sar',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  array\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'array',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  radar\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'radar',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  post\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'post',oparams,fid);
  fprintf(fid,'<strong>%s\n','='*ones(1,80)); fprintf(fid,'  analysis\n'); fprintf(fid,'%s</strong>\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'analysis',oparams,fid);
  fprintf(fid,'\n');
  
  param_txt_fn = ct_filename_ct_tmp(param,'','param', [param.config.date_str,'.txt']);
  fprintf('Writing %s\n\n', param_txt_fn);
  param_txt_fn_dir = fileparts(param_txt_fn);
  if ~exist(param_txt_fn_dir,'dir')
    mkdir(param_txt_fn_dir);
  end
  [fid,msg] = fopen(param_txt_fn,'wb');
  if fid<0
    error('Could not write to %s: %s\n', param_txt_fn, msg);
  end
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  cmd\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'cmd',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  records\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'records',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  qlook\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'qlook',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  sar\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'sar',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  array\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'array',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  radar\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'radar',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  post\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'post',oparams,fid);
  fprintf(fid,'\n');
  fprintf(fid,'%s\n','='*ones(1,80)); fprintf(fid,'  analysis\n'); fprintf(fid,'%s\n','='*ones(1,80));
  read_param_xls_print(param.config.param_fn,'analysis',oparams,fid);
  fprintf(fid,'\n');
  fclose(fid);
end

%% Exit task
% =========================================================================
fprintf('%s done %s\n', mfilename, datestr(now));

success = true;
