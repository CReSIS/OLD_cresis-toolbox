% Script channel_delay_analysis.m
%
% This script is used to determine the delay between channels.
%
%
% Authors: Logan Smith, John Paden
%

tic;
physical_constants;

% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Data loaded
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
fig_offset = 0;
% 2009 Antartica DC-8: 2009/10/31 seg 2 frame 46 chunk 1

% The base path for all the data
base_path = '/cresis/scratch2/mdce/MCoRDS/2010_Greenland_DC8/CSARP_out/';

% Information for building the path to the file
param.year     = 2010;
param.month    = 04;
param.day      = 20;
param.seg      = 3;
param.frame    = 3;
param.proc_dir = '';

% Start and stop chunk to load (inf for second element loads to the end)
param.chunk = [1 inf];

% Offset in the filename where the date is located (once filenames
% are standardized this should be removed)
param.date_offset = 1;

% Waveforms to load
param.wfs       = 1;
% Receiver channels to load
param.rxs       = [1 2 3 4 5 6 7 8];

sig_rxs = [1 2 3 4 5];

param.wfs       = 1;

% Combine waveforms parameters
param.wf_comb = 10e-6;


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
  sprintf('%04d%02d%02d_seg%d',param.year,param.month,param.day,param.seg), ...
  param.proc_dir, ...
  sprintf('fk_data_%02d_01/',param.frame));

if ~exist('channel_delay_analysis_loaded','var') || isempty(channel_delay_analysis_loaded) || ~channel_delay_analysis_loaded ...
    || ~strcmp(param.proc_dir,current_param.proc_dir) || param.month ~= current_param.month || param.day ~= current_param.day ...
    || param.seg ~= current_param.seg || param.frame ~= current_param.frame
  
  [data,pos,data_param] = load_fk_data(param);
  
  channel_delay_analysis_loaded = true;
end

current_param = param;
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
% Cut out surface data chunks from waveform 1 and 2
% ---------------------------------------------------------------------
% ---------------------------------------------------------------------
if param.debug_level >= 1
  fprintf('Extracting chunks (%.2f sec)\n', toc);
end
param.binRng = -10:10;
param.lineRng = [-5:5];
param.dbin = length(param.binRng);
param.dline = length(param.lineRng);

surf = cell(2,1);

Depth = cell(1,length(data));
Time  = cell(1,length(data));
for wf_idx = 1:length(param.wfs)
  wf = param.wfs(wf_idx);
  Time{wf_idx} = data_param.wfs(wf).time;
end
height_idxs = zeros(1,length(pos.surface));
pos.ice_height = pos.surface*3e8/2/sqrt(3.15);

for wf = param.wfs
  % Normalize data to avoid precision limits of single
  data{wf} = data{wf} ./ max(abs(data{wf}(:)));
  % Determine which bins/lines will be processed
  Nt = size(data{wf},1);
  Nx = size(data{wf},2);
  bins = numel(param.binRng)/2+0.5 : param.dbin ...
    : Nt-(numel(param.binRng)/2-0.5);
  lines = numel(param.lineRng)/2+0.5 : param.dline ...
    : Nx-(numel(param.lineRng)/2-0.5);
  surf{wf} = zeros([length(param.binRng) Nx size(data{wf},3)],'single');
  for line = lines
    height_idxs(line) = find(pos.surface(line)<Time{1},1)-1;
    surf{wf}(:,line+param.lineRng,:) ...
      = data{wf}(height_idxs(line)+param.binRng, line+param.lineRng,:);
  end
end

interp_factor = 10;
dt = 1/(data_param.wfs(1).fs*interp_factor);
start_bin = 1*interp_factor;
stop_bin = size(data{1},1)*interp_factor;
rlines = 11200:12200;
interp_data = zeros([size(surf{1},1)*interp_factor length(rlines) size(surf{1},3)]);
figure(4+fig_offset); clf; hold on
c = {'.b','.r','.g','.c','.m','.y','.k','xb'};
for rx_idx = 1:length(param.rxs)
  rx = param.rxs(rx_idx);
  interp_data(:,:,rx_idx) = interpft(surf{1}(:,rlines,rx_idx),interp_factor*size(surf{1},1));
  plot(lp(interp_data(:,1,rx_idx)),c{rx})
  [max_vals max_idxs] = max(interp_data(:,:,rx_idx));
  max_idx(rx) = mean(max_idxs);
  quality(rx) = std(max_idxs)/interp_factor/sqrt(length(rlines));
end
mean_sig_chan = mean(max_idx(sig_rxs));
for rx_idx = 1:length(param.rxs)
  fprintf('Time delay %d: %f ns std %f\n', param.rxs(rx_idx), (max_idx(rx_idx) - mean_sig_chan)*dt/1e-9, quality(rx_idx)*dt/1e-9)
end
hold off

figure(1+fig_offset); clf; hold on
for rx=1:8
  plot(rx,max_idx(rx),c{rx},'MarkerSize',20)
end
hold off