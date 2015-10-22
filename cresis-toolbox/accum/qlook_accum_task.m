function [output_data,output_time,output_gps] = qlook_accum_task(fn,param)
% [output_data,output_time,output_gps] = qlook_accum_task(fn,param)
%
% Qlook processing task called by qlook_accum, but may also be called
% directly from the command line.
%
% One output argument (for calling from create_task or qlook_accum)
% output_data
%   Status boolean (true == success)
%
% Two or more output arguments (for calling from the command line)
% output_data
%   2-D array of single, first dim is fast-time, second dim is slow-time
% output_time
%   Fast-time vector
% output_gps
%   GPS data synchronized to data
%
% Author: John Paden, Sam Buchanan
%
% See also: master, qlook_accum, run_basic_load_accum

qlook_accum_task_tstart = tic;
output_data = false;
physical_constants;

% Setup GPS filename and load
if isfield(param.qlook, 'records_en') && param.qlook.records_en
  % GPS data can be found in records file, so load the whole thing up front
  records_fn = ct_filename_support(param,'','records');
  records = read_records_aux_files(records_fn, param.qlook.load_idxs, 1);
  trajectory_param = struct('rx_path', 1, ...
    'tx_weights', 1, 'lever_arm_fh', param.qlook.lever_arm_fh);
  records = trajectory_with_leverarm(records,trajectory_param);
else
  % Slow mode (i.e. records file/aux files do not exist)
  % NOTE: this mode uses uncorrected UTC times from the data header
  if isfield(param,'season_name')
    gps_fn = ct_filename_support(param,param.qlook.gps.fn,'gps',true);
  else
    gps_fn = '';
  end
  if exist(gps_fn,'file') && ~exist(gps_fn,'dir')
    gps = load(gps_fn);
  end
end

% Software presums (rename for convenience)
presums = param.qlook.sw_presums;

% =======================================================================
% Load and process blocks of data at a time to avoid memory overflow
%  - Make the block size a multiple of the presums
%  - Keep loading blocks until the loaded block is smaller than that
%    requested
% =======================================================================
block_size = round(4200/presums)*presums;
data_size_read = block_size;
data_idx = 1;
Data = [];
Latitude = [];
Longitude = [];
Elevation = [];
GPS_time = [];
while data_size_read == block_size
  
  % =======================================================================
  % Load data
  % =======================================================================
  fprintf('  Loading data %s (%.1f sec)\n', fn, toc(qlook_accum_task_tstart));
  finfo = fname_info_accum(fn);
  if finfo.file_idx == 0
    [hdr,data] = basic_load_accum(fn, struct('clk',param.radar.fs,'first_byte',2^26,'recs',[data_idx block_size]));
  else
    [hdr,data] = basic_load_accum(fn, struct('clk',param.radar.fs,'recs',[data_idx block_size]));
  end
  data_size_read = size(data,2);
  
  % Just get waveforms of interest
  data = data(:,:,param.radar.wfs);
  chan_weights = param.qlook.band_window_func(length(param.radar.wfs)).'./param.radar.chan_equal(param.radar.wfs);
  
  % Remove digital error bursts
  %  Probably should make these values programmable from spreadsheet, not
  %  sure if they are always 44047 and 3840
  test_idxs = find(data(:) == 44047).';
  for test_idx = test_idxs
    if test_idx >= 1 && data(test_idx-1) < data(test_idx) ...
        && test_idx < numel(data) && data(test_idx+1) < data(test_idx)
      data(max(1,test_idx-2):min(test_idx+1,numel(data))) ...
        = 2^(param.radar.adc_bits-1)*hdr.wfs(1).presums/2^hdr.wfs(1).bit_shifts;
    end
  end
  test_idxs = find(data(:) == 3840).';
  for test_idx = test_idxs
    if test_idx >= 1 && data(test_idx-1) < data(test_idx) ...
        && test_idx < numel(data) && data(test_idx+1) < data(test_idx)
      data(max(1,test_idx-2):min(test_idx+1,numel(data))) ...
        = 2^(param.radar.adc_bits-1)*hdr.wfs(1).presums/2^hdr.wfs(1).bit_shifts;
    end
  end
  
  % =====================================================================
  % Load the GPS data
  % =====================================================================
  fprintf('  Syncing GPS data (%.1f sec)\n', toc(qlook_accum_task_tstart));
  
  if isfield(param.qlook, 'records_en') && param.qlook.records_en
    % Use records file to load data if it has been created
    % Records data is presynced, so we can skip the sync steps
    
    Latitude = cat(2,Latitude,records.lat(data_idx:data_idx+data_size_read-1));
    Longitude = cat(2,Longitude,records.lon(data_idx:data_idx+data_size_read-1));
    Elevation = cat(2,Elevation,records.elev(data_idx:data_idx+data_size_read-1));
    GPS_time = cat(2,GPS_time,records.gps_time(data_idx:data_idx+data_size_read-1));
    param.qlook.gps.gps_source = records.gps_source;
    
  else
    % Slow mode (i.e. records file/aux files do not exist)
    % NOTE: this mode uses uncorrected UTC times from the data header
    if exist(gps_fn,'file') && ~exist(gps_fn,'dir')
      
      % Get the GPS seconds of day to sync to radar
      UTC_sod = epoch_to_sod(gps.gps_time - utc_leap_seconds(gps.gps_time(1)), param.day_seg(1:8));
      
      % Remove repeat values in time
      [UTC_sod sort_idxs] = unique(UTC_sod);
      gps.lat = gps.lat(sort_idxs);
      gps.lon = gps.lon(sort_idxs);
      gps.elev = gps.elev(sort_idxs);
      gps.gps_time = gps.gps_time(sort_idxs);
      
      
      % Check for seconds of day roll over and unwrap (assume jump backward
      % of more than 23 hours is a roll over)
      wrap_idxs = find(diff(hdr.utc_time_sod) < -(86400-3600));
      for wrap_idx = wrap_idxs
        hdr.utc_time_sod(wrap_idx+1:end) = hdr.utc_time_sod(wrap_idx+1:end) + 86400;
      end
      hdr.utc_time_sod = hdr.utc_time_sod + param.qlook.gps.time_offset;
      
      % When we roll over midnight on a previous block or file, we must
      % correct the time for all files after it.  We will assume that if
      % the current hdr.utc_time_sod is less than the first value of
      % UTC_sod than we need to add 24*3600 (one day) to the
      % hdr.utc_time_sod values
      if min(hdr.utc_time_sod) < min(UTC_sod)
          hdr.utc_time_sod = hdr.utc_time_sod + (24*3600);
      end
      
      % Correct occasional header time error
      epri = median(diff(hdr.utc_time_sod));
      medfilt1_utc_time_sod = medfilt1(hdr.utc_time_sod,5);
      if any(abs(medfilt1_utc_time_sod - hdr.utc_time_sod > epri))
        fprintf('Correcting header times\n');
        bad_idxs = logical(abs(medfilt1(hdr.utc_time_sod,5) - hdr.utc_time_sod > epri));
        hdr.utc_time_sod(bad_idxs) = medfilt1_utc_time_sod(bad_idxs);
      end
      
      Latitude = cat(2,Latitude,interp1(UTC_sod,gps.lat,hdr.utc_time_sod));
      Longitude = cat(2,Longitude,interp1(UTC_sod,gps.lon,hdr.utc_time_sod));
      Elevation = cat(2,Elevation,interp1(UTC_sod,gps.elev,hdr.utc_time_sod));
      GPS_time = cat(2,GPS_time,interp1(UTC_sod,gps.gps_time,hdr.utc_time_sod));
      param.qlook.gps.gps_source = gps.gps_source;
    else
      % GPS file does not exist, so fill in with zeros
      fprintf('    GPS file %s does not exist (writing radar UTC time as GPS time!)\n', gps_fn);
      Latitude = cat(2,Latitude,zeros(size(hdr.utc_time_sod)));
      Longitude = cat(2,Longitude,zeros(size(hdr.utc_time_sod)));
      Elevation = cat(2,Elevation,zeros(size(hdr.utc_time_sod)));
      GPS_time = cat(2,GPS_time,hdr.utc_time_sod);
      param.qlook.gps.gps_source = '';
    end
  end
  
  data_idx = data_idx + data_size_read;
  
  % =======================================================================
  % Convert from quantization to voltage @ ADC
  % =======================================================================
  fprintf('  Conversion from quantization to voltage (%.1f sec)\n', toc(qlook_accum_task_tstart));
  adc_bits = param.radar.adc_bits;
  Vpp_scale = param.radar.Vpp_scale;
  data = (data-mean(data(:))) * Vpp_scale/2^adc_bits * 2^hdr.wfs(1).bit_shifts / hdr.wfs(1).presums;
  
  % =======================================================================
  % Pulse compress parameters
  % =======================================================================
  pc_param.f0 = abs(param.radar.f0 - param.radar.fLO);
  pc_param.f1 = pc_param.f0 + param.radar.BW;
  pc_param.Tpd = param.radar.Tpd;
  Nt_ref = floor(pc_param.Tpd*param.radar.fs) + 1;
  dt = 1/param.radar.fs;
  pc_param.time = hdr.wfs(1).t0 - param.radar.td - (Nt_ref-1)*dt + dt*(0:(size(data,1)+Nt_ref-1)-1).';
  pc_param.tukey = param.qlook.tukey;
  pc_param.zero_pad = false;
  pc_param.td_window_func = param.qlook.td_window_func;
  pc_param.window_func = param.qlook.window_func;
  pc_param.decimate = false;
  
  % =======================================================================
  % Elevation compensation
  % =======================================================================
  if param.qlook.elev_comp.en
    fprintf('  Elevation compensation (%.1f sec)\n', toc(qlook_accum_task_tstart));
    new_elev = polyval(polyfit([1 length(Elevation)],Elevation([1 end]),1),1:length(Elevation));
    data = fft(permute(data,[1 3 2]),[],1);
    % Create frequency axis for time shifting
    for wf_idx = 1:length(param.radar.wfs)
      wf = param.radar.wfs(wf_idx);
      dt = 1/param.radar.fs;
      Nt = size(data,1);
      df = 1/(Nt*dt);
      freq_start = param.radar.fLO + (wf-1)*param.radar.f_step;
      freq(:,wf) = linspace(freq_start, freq_start-param.radar.fs+df, Nt);
    end
    
    % Apply time shift
    drange = new_elev - Elevation;
    for rline = 1:size(data,3)
      data(:,:,rline) = data(:,:,rline) .* exp(j*2*pi*freq/(c/2)*drange(rline));
    end
    data = ifft(permute(data,[1 3 2]),[],1);
    Elevation = new_elev;
  end
  
  % =======================================================================
  % Presum data
  % =======================================================================
  fprintf('  Presums/coherent averages (%.1f sec)\n', toc(qlook_accum_task_tstart));
  num_rec = floor(size(data,2)/presums);
  for wf=1:size(data,3)
    data(:,1:num_rec,wf) = fir_dec(data(:,:,wf),presums);
  end
  data = data(:,1:num_rec,:);
  
  % =======================================================================
  % Pulse compress
  % =======================================================================
  fprintf('  Pulse compressing (%.1f sec)\n', toc(qlook_accum_task_tstart));
  
  % zero-pad for pulse compression and elevation compensation
  data = cat(1,zeros(Nt_ref-1,size(data,2),size(data,3)),data);
  
  clear pc_data;
  for wf_idx = 1:size(data,3)
    [pc_data(:,:,wf_idx),pc_time] ...
      = pulse_compress(data(:,:,wf_idx),pc_param);
  end
  clear data;
  
  % =======================================================================
  % Combine waveforms
  % =======================================================================
  fprintf('  Combining waveforms (%.1f sec)\n', toc(qlook_accum_task_tstart));
  
  dt_out = 1/(size(pc_data,3)*abs(param.radar.f_step));
  Time = pc_time(1):dt_out:pc_time(end).';
  Depth = Time*c/2;
  
  if 0
    % Code for finding and testing chan_weights
    psd = squeeze(mean(abs(fft(pc_data)).^2,2));
    plot(lp(psd))
    chan_weights = 1 ./ sqrt(mean(psd(170:440,:)));
    chan_weights = chan_weights ./ max(chan_weights);
    chan_weights
  end
  
  % Apply channel coefficients for 16 waveforms
  for wf_idx = 1:size(pc_data,3)
    pc_data(:,:,wf_idx) = pc_data(:,:,wf_idx) * chan_weights(wf_idx);
  end
  
  % 20 MHz spaced channels are like network analyzer with unambiguous range
  % of 1/20 MHz and num_wf samples in between, taking fft is like
  if param.radar.f_step < 0
    pc_data = fft(pc_data(:,:,end:-1:1),[],3);
  else
    pc_data = fft(pc_data,[],3);
  end
  
  % Transform to fast-time frequency domain
  pc_data = fft(pc_data,[],1);
  
  % Complex baseband data by removing conjugate symmetric upper half
  pc_data = pc_data(1:size(pc_data,1)/2,:,:);
  
  % Apply fast-time hanning window to data
  % pc_data = pc_data .* repmat(hanning(size(pc_data,1)), [1 size(pc_data,2) size(pc_data,3)]);
  
  % Fourier interpolate data to output time axis
  pc_data = ifft([zeros(length(Time)-size(pc_data,1),size(pc_data,2),size(pc_data,3)); pc_data],[],1);
  
  offset = 6;
  % DEBUG CODE FOR CHECKING offset
  % figure(1); clf;
  % for offset=0:size(data,3)-1
  
  new_data = zeros(size(pc_data,1),size(pc_data,2));
  for bin = 1:size(pc_data,1)
    new_data(bin,:) = pc_data(bin,:,mod(offset+bin-1,size(pc_data,3))+1);
  end
  Data = cat(2,Data,new_data);
  clear pc_data;
end

% % DEBUG CODE FOR CHECKING offset
% plot(lp(Data(:,400)))
% hold on;
% offset
% pause
% end

% =======================================================================
% Decimate GPS data
% =======================================================================

Latitude = fir_dec(Latitude,param.qlook.sw_presums);
Latitude = filter2(ones(1,param.qlook.incoh_ave(2))/param.qlook.incoh_ave(2), Latitude,'valid');
Latitude = Latitude(1:param.qlook.decimate:end);

Longitude = fir_dec(Longitude,param.qlook.sw_presums);
Longitude = filter2(ones(1,param.qlook.incoh_ave(2))/param.qlook.incoh_ave(2), Longitude,'valid');
Longitude = Longitude(1:param.qlook.decimate:end);

Elevation = fir_dec(Elevation,param.qlook.sw_presums);
Elevation = filter2(ones(1,param.qlook.incoh_ave(2))/param.qlook.incoh_ave(2), Elevation,'valid');
Elevation = Elevation(1:param.qlook.decimate:end);

GPS_time = fir_dec(GPS_time,param.qlook.sw_presums);
GPS_time = filter2(ones(1,param.qlook.incoh_ave(2))/param.qlook.incoh_ave(2), GPS_time,'valid');
GPS_time = GPS_time(1:param.qlook.decimate:end);

% =======================================================================
% Plot results
% =======================================================================
fprintf('  Create quick look outputs (%.1f sec)\n', toc(qlook_accum_task_tstart));

% Incoherent averaging and decimation
if mod(param.qlook.incoh_ave(1),2) == 0
  error('Number of incoherent averages in fast-time must be odd, %d', param.qlook.incoh_ave(1));
end
Data = filter2(ones(param.qlook.incoh_ave)/prod(param.qlook.incoh_ave), abs(Data).^2,'valid');
Data = Data(:,1:param.qlook.decimate:end);
clear data;
valid_removed = (param.qlook.incoh_ave-1)/2;
Time = Time(1+valid_removed : end-valid_removed);
Depth = Depth(1+valid_removed : end-valid_removed);
if ~isempty(param.qlook.detrend_poly_order)
  median_Data = median(lp(Data,2),2);
  p_median_Data = polyfit((1:size(Data,1)).', median_Data, 1);
  pf_median_Data = polyval(p_median_Data, (1:size(Data,1)).');
  replace_idxs = find(median_Data > pf_median_Data+3);
  median_Data(replace_idxs) = pf_median_Data(replace_idxs)+3;
  replace_idxs = find(median_Data < pf_median_Data-3);
  median_Data(replace_idxs) = pf_median_Data(replace_idxs)-3;
  p_median_Data = polyfit((1:size(Data,1)).', median_Data, param.qlook.detrend_poly_order);
  pf_median_Data = polyval(p_median_Data, (1:size(Data,1)).');
  if 0 % Debug for detrending data
    figure(1); clf;
    plot(median_Data)
    hold on;
    plot(pf_median_Data,'r')
    hold off;
  end
  
  Data = Data./repmat(10.^(pf_median_Data/20),[1 size(Data,2)]);
end

% Convert time min_bin into range bins
min_bin = find(Time > param.qlook.surf.min_bin, 1);

% Set initial point
%  - Look at the maximum along Ninit_pnts equally spaced range lines
%    and choose the sort_ind bin. The range line and corresponding
%    bin become the initial point.
Ninit_pnts = 15;
sort_ind = 3;
startInds = unique(round(size(Data,2) * linspace(0.2,0.8,Ninit_pnts)));
surfBins = zeros(1,size(Data,2));
[tmp surfBins_init] = max(Data(min_bin:end,startInds));
surfBins_init = surfBins_init + min_bin - 1;
[tmp surfBins_init_sort_ind] = sort(surfBins_init);
startInd = startInds(surfBins_init_sort_ind(sort_ind));
surfBins(startInd) = surfBins_init(surfBins_init_sort_ind(sort_ind));
pnt(1).col = startInd;
pnt(1).row = surfBins(startInd);
pnt(1).method = 's';

% Automated: find remaining points (snake method)
searchRng = param.qlook.surf.search_rng; % Must be centered on zero, as in -M:M
done = zeros(1,size(Data,2));
for pntInd = 1:length(pnt)
  done(pnt(pntInd).col) = 1;
  surfBins(pnt(pntInd).col) = pnt(pntInd).row;
end
while sum(done) ~= length(done)
  for line = 1:length(done)
    if done(line) == 0
      if line < length(done) && done(line+1) == 1
        done(line) = 1;
        tmpSearchRng = searchRng(1+max(0,1-(surfBins(line+1)+searchRng(1))) : ...
          length(searchRng)-max(0,(surfBins(line+1)-searchRng(1))-size(Data,1)));
        [tmp newBin] = max(Data(surfBins(line+1)+tmpSearchRng,line));
        surfBins(line) = surfBins(line+1)+tmpSearchRng(1)-1 + newBin;
      end
      if line > 1 && done(line-1) == 1
        done(line) = 1;
        tmpSearchRng = searchRng(1+max(0,1-(surfBins(line-1)+searchRng(1))) : ...
          length(searchRng)-max(0,(surfBins(line-1)-searchRng(1))-size(Data,1)));
        [tmp newBin] = max(Data(surfBins(line-1)+tmpSearchRng,line));
        surfBins(line) = surfBins(line-1)+tmpSearchRng(1)-1 + newBin;
      end
    end
  end
end

Surface = interp1(1:length(Time),Time,surfBins);

Data = Data;

frm_id = sprintf('%s_%03.0f', param.day_seg, finfo.file_idx);

if param.qlook.out.en
  out_name = sprintf('Data_%s.mat', frm_id);
  out_dir = ct_filename_out(param, ...
    param.qlook.out.dir, 'CSARP_qlook');
  out_fn = fullfile(out_dir,out_name);
  if ~exist(out_dir,'dir')
    mkdir(out_dir)
  end
  
  fprintf('    %s\n', out_fn);
  param_qlook = param;
  save('-v6',out_fn,'Data','Time','Surface','Latitude','Longitude','Elevation','GPS_time','param_qlook');
end

if param.qlook.plot.en == 1 || param.qlook.plot.en == 2
  figure(300); clf;
  if param.qlook.plot.en == 2
    start_bin = 1;
    stop_bin = size(Data,1);
  else
    start_bin = max(min(surfBins)-100,10);
    stop_bin = min(min(surfBins)+900,size(Data,1));
  end
  imagesc(lp(abs(Data(start_bin:stop_bin,:)).^2));
  colormap(1-gray(256));
  hold on;
  plot(surfBins-start_bin+1,'b');
  hold off;
  drawnow;
elseif param.qlook.plot.en == 3
  if strcmp(param.post.img_type,'jpg')
    print_device = '-djpeg';
  elseif strcmp(param.post.img_type,'png')
    print_device = '-dpng';
  else
    error('Unsupported image type %s', param.post.img_type);
  end
  print_dpi = sprintf('-r%d', param.post.img_dpi);  % Create output directory
  image_dir = fullfile(ct_filename_out(param, ...
    param.post.out_dir, 'CSARP_post', true),'images',param.day_seg);
  if ~exist(image_dir,'dir')
    mkdir(image_dir)
  end
  time_stamp_str = datestr(epoch_to_datenum(GPS_time(1)),'HHMMSS');
  frm = finfo.file_idx;
  frm_concat_list = frm;
  
  mdata.Data = Data;
  mdata.Surface = Surface;
  mdata.Latitude = Latitude;
  mdata.Longitude = Longitude;
  mdata.Elevation = Elevation;
  mdata.GPS_time = GPS_time;
  mdata.Time = Time;
  
  echo_param.fig_hand = 2;
  echo_param.num_x_tics = 6;
  echo_param.frm_id = frm_id;
  echo_param.depth = param.post.depth_rng;
  echo_param.elev_comp = true;
  echo_param.er_ice = param.post.er_ice;
  
  lay.GPS_time = mdata.GPS_time;
  lay.Elevation = mdata.Elevation;
  lay.layerData{1}.value{1}.data = NaN*zeros(size(mdata.Surface));
  lay.layerData{1}.value{2}.data = mdata.Surface;
  lay.layerData{2}.value{1}.data = NaN*zeros(size(mdata.Surface));
  lay.layerData{2}.value{2}.data = NaN*zeros(size(mdata.Surface));
  echo_info = publish_echogram(echo_param,mdata,lay);
  set(echo_info.h_surf,'Visible','off');
  set(echo_info.h_bot,'Visible','off');
  if length(frm_concat_list) > 1
    title(sprintf('"%s" %s Frame IDs: %s', param.radar_name, param.post.mission_name, frm_id),'Interpreter','none');
  else
    title(sprintf('"%s" %s Frame ID: %s', param.radar_name, param.post.mission_name, frm_id),'Interpreter','none');
  end
  
  % Remove current file
  echo_fn = sprintf('%s*_1echo.%s',frm_id,param.post.img_type);
  echo_fn = fullfile(image_dir,echo_fn);
  delete(echo_fn);
  
  echo_fn = sprintf('%s_%s_1echo.%s',frm_id,time_stamp_str,param.post.img_type);
  echo_fn = fullfile(image_dir,echo_fn);
  set(echo_param.fig_hand,'PaperOrientation','Portrait');
  set(echo_param.fig_hand,'PaperPosition',[0.5 0.5 10 7.5]);
  print(echo_param.fig_hand,print_device,print_dpi,echo_fn);
end

fprintf('  Done (%.1f sec)\n', toc(qlook_accum_task_tstart));

if nargout >= 2
  % Called from the command line
  output_data = Data;
  output_time = Time;
  if nargout == 3
    output_gps.lat = Latitude;
    output_gps.lon = Longitude;
    output_gps.elev = Elevation;
    output_gps.gps_time = GPS_time;
  end
else
  % First argument is a status variable for create_task to use
  output_data = true;
  output_time = [];
end

return;
