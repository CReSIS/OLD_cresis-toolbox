% function airice_surface_stats(param)
% Function airice_surface_stats(param)
%
% This function is used to verify the motion compensation and to produce
% data dependent phase weights for the cross-track array.
%
% The primary spectral component from the air/ice surface is estimated
% using the root music algorithm. The root music algorithm is valid
% for linear uniform arrays. The current program also assumes half
% wavelength spacing when converting the root music spectral component
% into direction of arrival (doa) angle:
%   doa_angle = asin(freq / pi)
% If the root music algorithm is not used, the generalized music algorithm
% must be used (not currently implemented in this script). Note that
% the root music algorithm is faster and performs better...
%
% The phase and amplitude coefficients are estimated for the array...
% not completed...
%
% Authors: Logan Smith, John Paden
%

clear;
tic;
physical_constants;

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Data loaded
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------

if ~exist('data_loaded','var') || isempty(data_loaded) || ~data_loaded
  if 1
    % 2009 Antartica DC-8: 2009/10/31 seg 2 frame 46 block 1
    
    % The base path for all the data
    base_path = '/cresis/scratch2/mdce/mcords2/2011_Greenland_P3/CSARP_out';
    
    % Information for building the path to the file
    param.year = 2010;
    param.month = 05;
    param.day  = 12;
    param.seg  = 01;
    param.frame = 7;
    
    % Start and stop block to load (inf for second element loads to the end)
    param.chunk = [1 3];
    
    % Offset in the filename where the date is located (once filenames
    % are standardized this should be removed)
    param.date_offset = 1;

    % Reference channel for amplitude/phase comparisons
    param.ref_chan = 3;
    
    % caxis for echogram (iamgesc) plots (leave empty for default)
    param.caxis = [-270  -80];

    % (wf,adc) pairs to load
%     param.imgs = {[1 9; 1 10; 1 11; 1 12]};
%     param.imgs = {[1 13; 1 14; 1 15; 1 16]};
%     param.imgs = {[1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8]};
%     param.imgs = {[1 5; 1 6; 1 7; 1 8]};
%     param.imgs = {[2 1; 2 2; 2 3; 2 4]};
%     param.imgs = {[1 5; 1 6; 1 7; 1 8]};
%     param.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8; 2 9; 2 10; 2 11; 2 12; 2 13; 2 14; 2 15; 2 16]};
    param.imgs = {[2 2; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8]};
%     param.imgs = {[2 9; 2 10; 2 11; 2 12]};
%     param.imgs = {[2 13; 2 14; 2 15; 2 16]};
%     param.imgs = {[3 1; 3 2; 3 3; 3 4; 3 5; 3 6; 3 7]};
%     param.imgs = {[4 1; 4 2; 2 3; 4 3; 4 4; 4 5; 4 6; 4 7]};
%     param.imgs = {[2 1; 2 2; 3 3; 2 3; 2 4; 2 5; 2 6; 2 7; 2 8]};
%     param.imgs = {[1 1; 1 2; 1 3; 1 4; 1 5; 1 6; 1 7; 1 8; 3 1; 3 2; 3 3; 3 4; 3 5; 3 6; 3 7]};
%     param.imgs = {[3 1; 3 2; 3 3; 3 4; 3 5; 3 6; 3 7; 3 8; 2 1; 2 2; 2 3; 2 4; 2 5; 2 6; 2 7]};
    
    % lever arm function handle
    lever_arm_fh = @lever_arm_mcords_2011_greenland_P3_gravity;
    
    % Combine waveforms parameters
    param.wf_comb           = 10e-6;
  end
  
  % Debug level (1 = default)
  param.debug_level = 1;
  
  % Combine receive channels
  param.combine_channels = 0;
  
  % Take abs()^2 of the data (only runs if combine_channels runs)
  param.incoherent = 1;
  
  % Combine waveforms (only runs if incoherent runs)
  param.combine_waveforms = 0;
  
  % Parameters for local_detrend (cmd == 5 disables, only runs if incoherent runs)
  param.detrend.cmd = 5;
  param.detrend.B_noise = [100 200];
  param.detrend.B_sig = [10 20];
  param.detrend.minVal = -inf;
  
  % =======================================================================
  % Automated
  % =======================================================================
  
  % Path to the input data
  param.in_path = fullfile(base_path, ...
    sprintf('%04d%02d%02d_%02d',param.year,param.month,param.day,param.seg), ...
    sprintf('fk_data_%03d_01_01/',param.frame));
  
  [data,pos,data_param,wfs] = load_fk_data(param);
  
  % =======================================================================
  % airice_surface_stats load specific work
  % =======================================================================
  
%   eqs.ang = deg2rad(1*[145 -127 0 60]);
%   eqs.pow_dB = [1.7194    1.0352         0    5.6808];
%   eqs.ang = deg2rad(1*[151 0 -49 -189]);
%   eqs.pow_dB = [-0.3035   -1.4765         0    3.3551];
%   eqs.ang = deg2rad(1*[145 -127 0 60 151 0 -49 -189]);
%   eqs.pow_dB = [1.7194    1.0352         0    5.6808 0.3474   -0.8256    0.6509    4.0060];
%   eqs.ang = deg2rad(1*[10 0 -20 -40]);
%   eqs.pow_dB = [1.8899    0.7598         0    0.3021];
  %eqs.ang = deg2rad(1*[-60 20 -10 0 -80 -10 -10]);
  eqs.ang = deg2rad(1*[20 -10 0 -80 -10 -60 -10]);
  eqs.pow_dB = [0 0 0 0 0 0 0];
  
  equalize_data = false;
  if equalize_data
    for img_idx=1:length(data)
      for wf_adc_idx=1:size(data{1},3)
        data{img_idx}(:,:,wf_adc_idx) ...
          = data{img_idx}(:,:,wf_adc_idx)/[10.^(eqs.pow_dB(wf_adc_idx)/20).*exp(1i*eqs.ang(wf_adc_idx))].';
      end
    end
  end
  
  roll_only_compensation = false;
  if roll_only_compensation
    for img_idx = 1:length(data)
      for wf_adc_idx = 1:size(data{img_idx},3)
        wf = param.imgs{img_idx}(wf_adc_idx,1);
        adc = param.imgs{img_idx}(wf_adc_idx,2);
        
        freq = wfs(wf).freq;
        
        roll_filt = pos.roll;
        
        lev_arm = lever_arm_fh(wfs(wf).tx_weights,wfs(wf).rx_paths(adc));
        
        drange = lev_arm(2).*tan(roll_filt);
        
        % Convert to time (in air)
        dtime = 2*drange/3e8;
        
        Nt = size(data{img_idx},1);
        Nx = size(data{img_idx},2);
        % Time shift data in the frequency domain
        data{img_idx}(:,:,wf_adc_idx) = ifft(fft(data{img_idx}(:,:,wf_adc_idx)) ...
          .*exp(-1i*2*pi*repmat(freq,1,Nx).*repmat(dtime,Nt,1)));
      end
    end
  end
  
  if 0
    % View data
    view_data{1} = abs(mean(data{1},3)).^2;
    view_data{2} = abs(mean(data{2},3)).^2;

    % Copied waveform combining code from load_fk_data:
    param.wf_bins = [320 370];
    param.wf_bin_comp = [400:415];
    difference = mean(mean(view_data{1}(param.wf_bin_comp,:))) ...
      ./ mean(mean(view_data{2}(param.wf_bin_comp,:)));

    trans_bins = param.wf_bins(1)+1:param.wf_bins(2);
    weights = 0.5+0.5*cos(pi*linspace(0,1,length(trans_bins)).');
    view_data = [view_data{1}(1:param.wf_bins(1),:); ...
      repmat(weights,[1 size(view_data{1},2)]).*view_data{1}(trans_bins,:) ...
      + difference*repmat(1-weights,[1 size(view_data{2},2)]).*view_data{2}(trans_bins,:); ...
      difference*view_data{2}(param.wf_bins(2)+1:end,:)];
  end
  
  data_loaded = 1;
end

% ---------------------------------------------------------------------
% Create Depth and Time axes for the data
Depth = cell(1,length(data));
Time  = cell(1,length(data));
for img_idx = 1:length(param.imgs)
  Time{img_idx} = wfs(param.imgs{1}(1,1)).time;
  Depth{img_idx} = c/2 * Time{img_idx};
end

if 0
  % Use the current surface variable
  % 	Convert to range bin units instead of fast-time
  pos.surface_bin = interp1(Time{1}, 1:length(Time{1}), pos.surface);

else
  % Track a layer and use that as the surface variable
  tmp_data = lp(filter2(ones(1,10)/10,mean(abs(data{1}(:,:,:)).^2,3)));
  search_bins = 670:800;
%   search_bins = 800:1000;
  [max_vals,max_bins] = max(tmp_data(search_bins,:));
  max_bins = max_bins + search_bins(1)-1;
  max_bins = medfilt1(max_bins,21);
  imagesc(tmp_data);
  hold on;
  plot(max_bins,'-k');
  hold off;
  
  pos.surface_bin = max_bins;
  pos.surface = interp1(1:length(Time{1}), Time{1}, pos.surface_bin);
end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Plot the data
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
if 1
  for wf_adc_idx = 1:size(data{1},3)
    figure(wf_adc_idx); clf;
    imagesc(lp(mean(abs(data{1}(:,:,wf_adc_idx)).^2,3)));
    hold on;
    plot(pos.surface_bin,'k')
    grid on;
    hold off;
    if ~isempty(param.caxis)
      caxis(param.caxis);
    end
  end
end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Cut out surface data blocks from waveform 1 and 2
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
if param.debug_level >= 1
  fprintf('Extracting blocks (%.2f sec)\n', toc);
end

% Each output of the processing uses a block of data around the input
% pixel (neighboring pixels) to compute the angle of arrival and other
% statistics. This block is defined by binRng and lineRng
param.binRng = -5:5;
param.lineRng = -5:5;
% Spacing of output pixels is determined by dbin and dline, these
% are set so that there is no overlap in the blocks of data
param.dbin = length(param.binRng);
param.dline = length(param.lineRng);

surf = cell(2,1);

height_idxs = zeros(1,length(pos.surface));
pos.ice_height = pos.surface*3e8/2/sqrt(3.15);

for img_idx = 1:length(param.imgs)
  % Normalize data to avoid precision limits of single
  data{img_idx} = data{img_idx} ./ max(abs(data{img_idx}(:)));
  % Determine which bins/lines will be processed
  Nt = size(data{img_idx},1);
  Nx = size(data{img_idx},2);
  bins = numel(param.binRng)/2+0.5 : param.dbin ...
    : Nt-(numel(param.binRng)/2-0.5);
  lines = numel(param.lineRng)/2+0.5 : param.dline ...
    : Nx-(numel(param.lineRng)/2-0.5);
  surf{img_idx} = zeros([length(param.binRng) Nx size(data{img_idx},3)],'single');
  for line = lines
    height_idxs(line) = find(pos.surface(line)<Time{img_idx},1)-1;
    surf{img_idx}(:,line+param.lineRng,:) ...
      = data{img_idx}(height_idxs(line)+param.binRng, line+param.lineRng,:);
  end
end

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Determine primary spectral component from surface return
%  - Primary spectral component is uniquely related to angle of
%    arrival by asin(normalized_freq / pi).
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
if param.debug_level >= 1
  fprintf('Spectral analysis (%.2f sec)\n', toc);
end
f = cell(length(param.imgs),1);
for img_idx=1:length(param.imgs)
  f{img_idx} = zeros(1,length(lines));
  for lineInd = 1:length(lines)
    line = lines(lineInd);
    fblock = surf{img_idx}(:,line+param.lineRng,:);
    % Align each snapshot (SAR pixel) to be a row into the spectral estimator
    fblock = reshape(fblock,[size(fblock,1)*size(fblock,2) size(fblock,3)]);
    % Track power level from each channel
    pow_levels{img_idx}(:,lineInd) = mean(abs(fblock).^2);
    % Remove the values with weak zero-frequency (~nadir when roll = 0) return
    powVals = abs(mean(fblock,2)).^2;
    [powVals sortInds] = sort(powVals);
    peakVal{img_idx}(lineInd) = max(powVals(1)); % Save this value for later
    % Grab the strongest returns
    fblockGood = fblock(sortInds(end-3*size(fblock,2)+1:end),:);
    % Rootmusic to estimate primary spectral component
    f{img_idx}(lineInd) = rootmusic(double(fblockGood),1);
    % For some reason
    % there is a sign change (have verified this, by plugging in pure
    % frequencies directly into rootmusic):
    %   rootmusic(exp(j*[0.35 0.05 0 -0.25; 0.35 0.05 0 -0.25]),1)
    % When only a single snapshot is given, the sign is correct:
    %   rootmusic(exp(j*[0.35 0.05 0 -0.25]),1)
    if size(fblockGood,1) > 1
      f{img_idx}(lineInd) = -f{img_idx}(lineInd);
    end
    % Reference all good blocks to reference channel
    fblockGood = fblockGood ./ repmat(fblockGood(:,param.ref_chan),[1 size(fblockGood,2)]);
    peakAngles{img_idx}(lineInd,:) = mean(fblockGood,1);
  end
end

mean_pow_levels = 10*log10(mean(pow_levels{1}./repmat(pow_levels{1}(param.ref_chan,:),[size(pow_levels{1},1) 1]),2)).';
  fprintf('Mean Power Levels\n');
  fprintf('  %.2f\t', mean_pow_levels);
  fprintf('\n');

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Plot results
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
if param.debug_level >= 1
  fprintf('Plotting results (%.2f sec)\n', toc);
end

% Assumes that antennas are ordered left to right
lambda = 3e8/wfs(param.imgs{1}(1,1)).fc;
d = 15.5*0.0254; % DC-8
d = 29.5333*0.0254; % P-3 wing
threshold = -inf;
% rollEst = asind(2*f{1}/pi) - 0;
rollEst = asind(lambda*f{1}/(2*pi*d)) - 0;
[Biir,Aiir] = butter(4,0.2);
rollEst = filtfilt(Biir,Aiir,rollEst);
rollEst(10*log10(peakVal{1}) < threshold) = NaN;
rollError = rollEst(10*log10(peakVal{1}) >= threshold) ...
  - pos.roll(10*log10(peakVal{1}) >= threshold);
fprintf('Roll error mean %.2f, std %.2f\n', mean(rollError), std(rollError));
if length(data) == 2
  rollEst2 = asind(2*f{2}/pi) - 0;
  rollEst2 = filtfilt(Biir,Aiir,rollEst2);
  rollEst2(10*log10(peakVal{2}) < threshold) = NaN;
end

figure(21); clf;
if isfield(pos,'roll')
  plot(rad2deg(pos.roll),'k');
else
  warning('No roll data to plot.');
  plot(0,'k');
end
hold on;
plot(lines,rollEst,'b-x');
plot(rad2deg(pos.roll),'k');
hold off;
axis tight;
grid;
if roll_only_compensation
  title('Roll with roll-only mocomp');
else
  title('Roll without mocomp');
end
ylabel('Roll (deg)')
ylim([-20 20]);
legend('INS Roll','Est Roll');
%print(gcf,'-dpng',sprintf('/cresis/scratch1/lsmith/Figures/MoComp/MCoRDS_%02d%02d_seg%02d_%02d_%s_roll.png',...
%  param.month,param.day,param.seg,param.frame,param.proc_dir))

if 1
  fprintf('Mean Phase Offset\n');
  fprintf('  %.2f\t', angle(mean(peakAngles{1}))*180/pi);
  fprintf('\n');
  figure(22); clf;
  %   plot(medfilt1(angle(peakAngles{1}.*exp(1j.*repmat(f{1}.',1,5).*repmat(0:size(fblockGood,2)-1,length(f{1}),1))),11)*180/pi);
  hphase = plot(medfilt1(angle(double(peakAngles{1})),11)*180/pi);
  hold on
  hroll = plot(rad2deg(pos.roll(lines)),'.k')
  hold off
  axis tight;
  grid on;
  if roll_only_compensation
    title('Angles/Roll with roll-only mocomp');
  else
    title('Angles/Roll without mocomp');
  end
  label_strs = {};
  for adc_wf_idx = 1:size(param.imgs{1},1)
    rx = wfs(param.imgs{1}(adc_wf_idx,1)).rx_paths(param.imgs{1}(adc_wf_idx,2));
    label_strs{end+1} = sprintf('Adc %d/Rx %d', ...
      param.imgs{1}(adc_wf_idx,2), rx);
  end
  label_strs{end+1} = 'Roll';
  legend([hphase; hroll], label_strs);
  %   ylim([-50 50])
  %print(gcf,'-dpng',sprintf('/cresis/scratch1/lsmith/Figures/MoComp/MCoRDS_%02d%02d_seg%02d_%02d_%s_peak_angles_wf_1.png',...
  %  param.month,param.day,param.seg,param.frame,param.proc_dir))
end
if 0
  figure(4); clf;
  plot(medfilt1(angle(peakAngles{2}),11)*180/pi);
  hold on
  plot(rad2deg(pos.roll(lines)),'.k')
  hold off
  axis tight;
  grid on;
  title('Peak Angles Waveform 2');
  legend('Rx1','Rx2','Rx3','Rx6','Rx7','Roll')
  %   ylim([-50 50])
  print(gcf,'-dpng',sprintf('/cresis/scratch1/lsmith/Figures/MoComp/MCoRDS_%02d%02d_seg%02d_%02d_%s_peak_angles_wf_2.png',...
    param.month,param.day,param.seg,param.frame,param.proc_dir))
end
% round(median(angle(peakAngles{1}))*180/pi)
% round(median(angle(peakAngles{2}))*180/pi)


if 0
  figure(6); clf;
  %imagesc([1 size(view_data,2)], Depth, 10*log10(local_detrend(view_data,[50 300],[4 50],4)),[-3 1]);
  imagesc([1 size(view_data,2)], Depth{1}, 10*log10(view_data))
  hold on;
  plot(pos.ice_height-mean(pos.ice_height));
  hold off;
  print(gcf,'-dpng',sprintf('/cresis/scratch1/lsmith/Figures/MoComp/MCoRDS_%02d%02d_seg%02d_%02d_%s_echo_wf_1.png',...
    param.month,param.day,param.seg,param.frame,param.proc_dir))
  figure(7); imagesc([],Depth{2},lp(data{2}(:,:,1)))
  hold on;
  plot(pos.ice_height-mean(pos.ice_height));
  hold off;
  print(gcf,'-dpng',sprintf('/cresis/scratch1/lsmith/Figures/MoComp/MCoRDS_%02d%02d_seg%02d_%02d_%s_echo_wf_2.png',...
    param.month,param.day,param.seg,param.frame,param.proc_dir))
end

return;
