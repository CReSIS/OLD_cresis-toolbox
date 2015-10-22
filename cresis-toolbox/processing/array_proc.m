function [array_param,dout] = array_proc(array_param,din,dout)
% [array_param,dout] = array_proc(array_param,din,dout)
%
% Performs array processing on the input data, din, and writes the result
% to dout. Designed specifically for cresis toolbox radar data.
%
% Requires the optimization toolbox for fmincon (constrained optimization)
% when using methods 7, 8, and 9.
%
% TODO: Support steering vector tables.
%
% =========================================================================
% INPUTS
% din is complex Nt by Nx by Na by Nb by Nc matrix
%   Nt is fast-time
%   Nx is slow-time
%   Nc is channels in array
%   Na is subaperture
%   Nb is subband
% array_param
%  .method
%    0: Periodogram (aka Welch or DFT) method
%    1: Minimum Variance Distortionless Response (MVDR)
%    2: Multiple Signal Classification (MUSIC)
%    3: Eigenvector method based on Matlab's peig NOT WORKING
%    4. Re-iterative super resolution, Blunt (RISR)
%    6. Minimum Variance Distortionless Response (Robust MVDR)
%    7. Maximum Likelihood Estimator (MLE)
%    8. DCM Correlation Method, Stumpf
%    9. Wideband Maximum Likelihood Estimator
%    Should be a scalar, defaults to Periodogram
%  .window
%    Periodogram only: window to be applied to din before the FFT
%    Should be a function handle if set, defaults to boxcar
%  .sv
%    Steering vectors (should be empty or not exist if uniformly spaced,
%    linear array is used). An example of how to produce the sv is included
%    at the bottom of the file. When sv is set, array processing can be much
%    slower since the FFT cannot be used.
%  .diag_load
%    MVDR only: Diagonal loading since MVDR requires matrix inversion.
%    This is especially relevant when length(bin_rng)*length(rline_rng) < Nc
%    Should be a scalar, defaults to zero
%  .Nsig
%    Method MUSIC, EIG, MLE, WB, and WBMLE: Expected number of signals in
%    the data. For a
%    surface with no layover this is 2, one for left and one for right.
%    For depth sounding with air/ice clutter, this should be 3 for bed
%    sounding only and 4 for 3D tomography. The increase in two is to
%    handle left and right clutter.
%    Should be a scalar, defaults to 3.
%  .Nsv = 64;
%    Choose how many spatial frequencies will be estimated by the beam
%    forming process. Given a normalized frequency axis (-0.5 to 0.5), this
%    is like freq = -0.5:1/Nsv:0.5-1/Nsv. The larger Nsv is the
%    finer the spatial frequency sampling will be.
%    Should be a scalar, defaults to 64
%  .bin_rng:
%    Vector of relative range bin indices to use in multilooking (aka snapshots).
%    Defaults to [-2:2] which causes five range bins to be used in the
%    generation of the array processing output of each pixel.  This means
%    two pixels before (-2) to two pixels after (2) the current pixel.
%  .rline_rng
%    Vector of relative range line indices to use in multilooking (aka snapshots).
%    Defaults to [-5:5] which causes 11 range lines to be used in the
%    generation of the array processing output of each pixel.
%  .debug_level
%    A debug level of 1 or greater produces stdout output and requires
%    tic to have been run before the command.
%    Should be a scalar, defaults to zero (no stdout)
%  .dbin
%    Specify spacing of range bins, setting to one causes every bin to be
%    output
%    Should be a positive integer, defaults to length of bin_rng
%  .dline
%    Specify spacing of range bins, setting to one causes every bin to be
%    output
%    Should be a positive integer, defaults to length of rline_rng
%  .rlines
%    2 element vector which specifies the first and last output range line.
%    This is an optional that overrules the "support" parameter.
%    It is used by combine_wf_chan_task.m to make sure all of the
%    chunks of SAR data run through array_proc can be seamlessly stiched
%    in combine_wf_chan.
%  .freq_rng
%    Vector containing each spatial frequency to be used when computing
%    the max value for each range bin.  If array_param.Nsv = 64, then
%    33 is typically broad side, 1 is looking end fire, 64 is looking
%    end fire in the opposite direction. The default is broadside:
%    array_param.freq_rng = round(0.49*array_param.Nsv) : round(0.51*array_param.Nsv);
%  .chan_equal: Nc by 1 vector. Default is all ones (so no effect on din).
%    data_used_in_array_processing = din / chan_equal
%  .support: string that describes the required data input support for each
%     output. Default is 'full'.
%     'short': drops the last few lines to make Nx a multiple of dline
%              and behaves like 'same' otherwise.
%     'same': output lines are array_param.dline/2+0.5 : array_param.dline : Nx;
%             Edges of data may have less support to create the correlation
%             matrix/power spectrum (e.g. correlation matrix estimate may
%             have reduced rank and less averaging)
%     'full': output lines are numel(array_param.rline_rng)/2+0.5 :
%             array_param.dline : Nx-(numel(array_param.rline_rng)/2-0.5);
%             The idea is that every output will have full support when
%             creating the correlation matrix/power spectrum.
%  .fcs: Flight Coordinate System cell vector used by array_param.sv_fh
%    Only requires the fields describing the sensor positions for each
%    multilook and snapshot. The roll must also be specified and will be
%    read from the first channel.
%    array_param.fcs{1...Ns}{1...Nc}.pos
%      Ns = snapshots
%      Nc = wf-adc pairs (i.e. sensors or channels)
%      pos(2,1...Nx): cross track with positive pointing to the left
%      pos(3,1...Nx): elevation with positive pointing up
%      roll
%    array_param.fcs{1...Ns}{1}.roll
%      Roll angle in degrees
%  .surface: Vector of time delays to the surface for each range line. This
%    is used:
%    1. in the construction of DOA constraints for methods 7-9
%    2. for plotting debug results
%    Default is NaN.
%  .three_dim: Structure describing the 3D processing (default is to not
%    assign)
%    .en: boolean to indicate the 3D outputs will be generated (default is disabled)
%  .sv_fh: Steering vector function handle
%  .fs: sampling frequency in Hz, used for method 8 to determine fast-time
%    impulse response
%  .fc: center frequency in Hz, used during steering vector generation
%  .time: used in DOA methods (7, 8, and 9) for creating the doa
%    constraints and for debug plots
%  .imgs: a list of the images to be used for incoherent array processing
%  .DCM: Structure with bin_rng and rline_rng fields that are used to
%     construct the data covariance matrix (DCM). Use this when you want the
%     DCM to be estimated from a different set of data than the
%     multilooking/snapshots.  The default is to have DCM.bin_rng and
%     DCM.rline_rng equal bin_rng and rline_rng.
%  .bin_restriction: 2 element structure array to limit array processing to
%    the specified range bins (main purpose is to speed up processing by
%    not processing all the range bins)
%    (1).bin: 1 by Nx vector specifying the start range bin for each range
%             line
%    (2).bin: 1 by Nx vector specifying the stop range bin for each range
%             line
%  .doa_constraints: structure array restricting the DOA for each source
%    .method: string indicating the constraint method
%      'fixed': a fixed range of DOAs
%      'surface-left': a fixed range of DOAs around a flat surface on the
%        left
%      'surface-right': a fixed range of DOAs around a flat surface on the
%        right
%      'layer-left': a fixed range of DOAs around a flat layer on the
%        left
%      'layer-right': a fixed range of DOAs around a flat layer on the
%        right
%      'er': dielectric to use for refraction to the layer
%    .init_src_limits: 2 element vector of minimum and maximum allowed
%       angle of arrival for each source (in degrees) during initialization
%    .src_limits: 2 element vector of minimum and maximum allowed
%       angle of arrival for each source (in degrees) during optimization
%  .reg_bins: bins to use for registration in method 8. This specifies the
%    bins that are to be used for the sinc interpolation that is done
%    to determine the pixel power after DOA estimation.
%  .theta_guard: the minimum source separation in radians. Used with DOA
%    methods. Default is 1.5 degrees.
%  .init: string specifying the initialization method:
%    'ap': alternating projection (not global)
%    'grid': sparse grid search (slower, but global)
%  .imp_resp: structure describing the radar's fast time impulse response.
%    Main lobe should be centered on zero time.
%    .vals: N by 1 vector of complex amplitudes
%    .time_vec: N by 1 vector of time
%  .W: widening factor used by wideband methods 8 and 9, analagous to
%    the number of subbands to use
%
% =========================================================================
% OUTPUTS
% array_param: same as input, but with some additional fields added
% dout is structure
%   .val = Nt/length(bin_rng) by Nx/length(rline_rng) real linear power
%     scale matrix
%   .freq = same size as .val and gives the corresponding spatial
%     frequency for each element in .val
%   .valR = Nsv by Nx/length(rline_rng) real linear power scale matrix
%   .bins = same size as .valR and gives the corresponding bin for each
%     element of .valR
%
% See also sim.crosstrack.m for an example of how to run.
    
% =====================================================================
%% Process and check input arguments
% =====================================================================
physical_constants;

% Number of fast-time samples in the din
Nt = size(din{1},1);

% Number of cross-track channels in the din
Nc = size(din{1},5);

% Number of subapertures
Na = size(din{1},3);

% Number of subbands
Nb = size(din{1},4);

if ~isfield(array_param,'method') || isempty(array_param.method)
  array_param.method = 0;
end

if ischar(array_param.method)
  if strcmpi(array_param.method,'standard')
    array_param.method = 0;
  elseif strcmpi(array_param.method,'period')
    array_param.method = 0;
  elseif strcmpi(array_param.method,'mvdr')
    array_param.method = 1;
  elseif strcmpi(array_param.method,'music')
    array_param.method = 2;
  elseif strcmpi(array_param.method,'eig')
    array_param.method = 3;
  elseif strcmpi(array_param.method,'risr')
    array_param.method = 4;
  elseif strcmpi(array_param.method,'geonull')
    array_param.method = 5;
  elseif strcmpi(array_param.method,'robust2')
    array_param.method = 6;
  elseif strcmpi(array_param.method,'mle')
    array_param.method = 7;
  elseif strcmpi(array_param.method,'wideband')
    array_param.method = 8;
  elseif strcmpi(array_param.method,'wbmle')
    array_param.method = 9;
  else
    error('Invalid method %s', array_param.method);
  end
end

if ~isfield(array_param,'window') || isempty(array_param.window)
  array_param.window = @boxcar;
end

if ~isfield(array_param,'diag_load') || isempty(array_param.diag_load)
  array_param.diag_load = 0;
end

if ~isfield(array_param,'Nsig') || isempty(array_param.Nsig)
  array_param.Nsig = 3;
end

if ~isfield(array_param,'bin_rng') || isempty(array_param.bin_rng)
  array_param.bin_rng = -2:1:2;
end

if ~isfield(array_param,'rline_rng') || isempty(array_param.rline_rng)
  array_param.rline_rng = -5:1:5;
end

if ~isfield(array_param,'debug_level') || isempty(array_param.debug_level)
  array_param.debug_level = 0;
end

if ~isfield(array_param,'dbin') || isempty(array_param.dbin)
  array_param.dbin = length(array_param.bin_rng);
end

if ~isfield(array_param,'dline') || isempty(array_param.dline)
  array_param.dline = length(array_param.rline_rng);
end

if ~isfield(array_param,'rlines') || isempty(array_param.rlines)
  array_param.rlines = [];
end

if ~isfield(array_param,'theta_rng') || isempty(array_param.theta_rng)
  array_param.theta_rng = [0 0];
end

if ~isfield(array_param,'chan_equal') || isempty(array_param.chan_equal)
  for ml_idx = 1:length(din)
    array_param.chan_equal{ml_idx} = ones(Nc,1);
  end
end

if ~isfield(array_param,'support') || isempty(array_param.support)
  array_param.support = 'full';
end

if ~isfield(array_param,'DCM') || isempty(array_param.DCM)
  array_param.DCM.bin_rng = array_param.bin_rng;
  array_param.DCM.rline_rng = array_param.rline_rng;
end

if ~isfield(array_param,'bin_restriction') || isempty(array_param.bin_restriction)
  restrict_bins = false;
else
  restrict_bins = true;
end

if ~isfield(array_param,'doa_constraints') || isempty(array_param.doa_constraints)
  for src_idx = 1:array_param.Nsig
    array_param.doa_contraints(src_idx).method = 'fixed';
    array_param.doa_contraints(src_idx).init_src_limits = [-90 90];
    array_param.doa_contraints(src_idx).src_limits = [-90 90];
  end
end

if ~isfield(array_param,'reg_bins') || isempty(array_param.reg_bins)
  array_param.reg_bins = -10:10;
end

if ~isfield(array_param,'theta_guard') || isempty(array_param.theta_guard)
  array_param.theta_guard = 1.5/180*pi;
end

if ~isfield(array_param,'init') || isempty(array_param.init)
  array_param.init = 'grid';
end

if ~isfield(array_param,'three_dim')
  array_param.three_dim.en = [];
end

if ~isfield(array_param.three_dim,'en') || isempty(array_param.three_dim.en)
  array_param.three_dim.en = false;
end

% ====================================================================
% Process the data
% =====================================================================

%% Check to see if data covariance matrix (DCM) generation uses same
% snapshots/pixels as multilooking
if length(array_param.bin_rng) == length(array_param.DCM.bin_rng) && all(array_param.bin_rng == array_param.DCM.bin_rng) ...
    && length(array_param.rline_rng) == length(array_param.DCM.rline_rng) && all(array_param.rline_rng == array_param.DCM.rline_rng)
  DCM_ML_match = true;
else
  DCM_ML_match = false;
end

%% Set up which range bins will be output bins
switch (array_param.method)
  case {0,1,2,3,4,5,6,7}
    % Narrowband methods
    array_param.bins = 1-min(array_param.bin_rng) : array_param.dbin ...
      : Nt-max(array_param.bin_rng);
  case {8}
    array_param.bins = 1 - min(array_param.bin_rng) - min(array_param.reg_bins) - floor((array_param.W-1)/2) : array_param.dbin ...
      : Nt - max(array_param.bin_rng) - max(array_param.reg_bins) + floor((array_param.W-1)/2);
    
  case {9}
    % Wideband MLE
    array_param.bins = 1-min(array_param.bin_rng - floor((array_param.W-1)/2)) : array_param.dbin ...
      : Nt-max(array_param.bin_rng + floor((array_param.W-1)/2));
end

%% Set up which range lines will be output lines

% Number of slow-time/along-track samples in the data
if strcmpi(array_param.support,'short')
  Nx = array_param.dline * floor(size(din{1},2)/array_param.dline)
else
  Nx = size(din{1},2);
end

if ~isempty(array_param.rlines)
  % Start/stop output range lines passed in (typical operation from
  % combine_wf_chan_task)
  array_param.lines = array_param.rlines(1): array_param.dline ...
    : min(array_param.rlines(2),size(din{1},2)-max(array_param.rline_rng));
elseif strcmpi(array_param.support,'full')
  array_param.lines = 1-min(array_param.rline_rng) : array_param.dline ...
    : Nx-max(array_param.rline_rng);
else
  array_param.lines = array_param.dline/2+0.5 : array_param.dline : Nx;
end

%% Preallocate/Initialize the output variables
% dout is output structure
%   .val = Nt/length(bin_rng) by Nx/length(rline_rng) real linear power
%     scale matrix
%   .freq = same size as .val and gives the corresponding spatial
%     frequency for each element in .val
if iscell(array_param.Nsv)
  theta = array_param.Nsv{2};
else
  theta = array_param.sv_fh(array_param.Nsv,array_param.wfs.fc);
end
Nsv = length(theta);
if ~exist('dout','var')
  dout.val ...
    = NaN*zeros(length(array_param.bins), length(array_param.lines),'single');
  dout.freq ...
    = NaN*zeros(length(array_param.bins),length(array_param.lines),'single');
  switch (array_param.method)
    case {7,8,9} % MLE DOA, Wideband DOA
      dout.doa = ...
        NaN*zeros(length(array_param.bins),array_param.Nsig,length(array_param.lines),'single');
      dout.hessian = ...
        NaN*zeros(length(array_param.bins),array_param.Nsig,length(array_param.lines),'single');
      dout.power = ...
        NaN*zeros(length(array_param.bins),array_param.Nsig,length(array_param.lines),'single');
      dout.cost = ...
        NaN*zeros(length(array_param.bins), length(array_param.lines),'single');
      if restrict_bins
        dout.tomo_top = ...
          NaN*zeros(1,length(array_param.lines),'single');
        dout.tomo_bot = ...
          NaN*zeros(1,length(array_param.lines),'single');
      end
      % dout.mle = ...
      %   NaN*zeros(length(array_param.bins),array_param.Nsig,length(array_param.lines),'single');
      
    otherwise % All other methods (beam-forming)
      if array_param.three_dim.en
        dout.img ...
          = NaN*zeros(length(array_param.bins),Nsv,length(array_param.lines),'single');
      end
  end
end

%% Apply channel equalization coefficients. Also apply
% window if periodogram.
Hwindow = array_param.window(Nc);
for ml_idx = 1:length(din)
  for chan = 1:Nc
    if array_param.method == 0
      % Periodogram method
      din{ml_idx}(:,:,:,:,chan) = din{ml_idx}(:,:,:,:,chan) * (Hwindow(chan) / array_param.chan_equal{ml_idx}(chan));
    else
      % All others
      din{ml_idx}(:,:,:,:,chan) = din{ml_idx}(:,:,:,:,chan) / array_param.chan_equal{ml_idx}(chan);
    end
  end
end

%% Setup parameterization for mle
if array_param.method == 7
  doa_param.fs              = array_param.wfs.fs;
  doa_param.fc              = array_param.wfs.fc;
  doa_param.Nsig            = array_param.Nsig;
  doa_param.options         = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);
  doa_param.doa_constraints = array_param.doa_constraints;
  doa_param.theta_guard     = array_param.theta_guard;
  doa_param.search_type     = array_param.init;
end

%% Setup parameterization for wideband doa estimator
if array_param.method == 8
  doa_param.h               = array_param.imp_resp.vals(:);
  doa_param.t0              = array_param.imp_resp.time_vec(1);
  doa_param.dt              = array_param.imp_resp.time_vec(2)-doa_param.t0;
  doa_param.fs              = array_param.wfs.fs;
  doa_param.fc              = array_param.wfs.fc;
  doa_param.Nsig            = array_param.Nsig;
  %doa_param.options         = optimset('Display','off'); %fminsearch
  doa_param.options         = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);
  doa_param.doa_constraints = array_param.doa_constraints;
  doa_param.theta_guard     = array_param.theta_guard;
  doa_param.search_type     = array_param.init;
end

%% Setup parameterization for wideband mle
if array_param.method == 9
  doa_param.fs              = array_param.wfs.fs;
  doa_param.fc              = array_param.wfs.fc;
  doa_param.Nsig            = array_param.Nsig;
  doa_param.nb_filter_banks = array_param.NB;
  doa_param.options         = optimoptions(@fmincon,'Display','off','Algorithm','sqp','TolX',1e-3);
  doa_param.doa_constraints = array_param.doa_constraints;
end

%% Setup Steering Vector variables
% Beam-forming approach: Determine which ky indices will be used in
% selecting the value for each range bin
theta_fftshift = fftshift(theta);
array_param.freq_rng = find(theta_fftshift >= array_param.theta_rng(1) & theta_fftshift <= array_param.theta_rng(2));
if isempty(array_param.freq_rng)
  [tmp array_param.freq_rng] = min(abs(theta_fftshift-mean(array_param.theta_rng)));
end


%% Array process each pixel by looping through each range line and each
%% range bin
tic;
line = 0;
for lineIdx = 1:1:length(array_param.lines)
  if array_param.debug_level >= 3
    % Insert debug code as needed
    if lineIdx == inf
      figure(1); clf;
      imagesc(Smusic);
      figure(2); clf;
      imagesc(20*log10(abs(squeeze(din{1}(:,line,:))).^2));
      keyboard
    end
  end
  
  line = array_param.lines(lineIdx);
  if array_param.debug_level >= 0
    fprintf('    Record %.0f (%.0f of %.0f) (%.1f sec)\n', line, lineIdx, ...
      length(array_param.lines), toc);
  end
  
  %% Check for edge conditions (not enough data to create dataSample
  % properly)
  if line+array_param.rline_rng(1) < 1
    rline_rng = 1-line : array_param.rline_rng(end);
  elseif  line+array_param.rline_rng(end) > Nx
    rline_rng = array_param.rline_rng(1) : Nx-line;
  else
    rline_rng = array_param.rline_rng;
  end
  if ~DCM_ML_match
    if line+array_param.DCM.rline_rng(1) < 1
      DCM_rline_rng = 1-line : array_param.DCM.rline_rng(end);
    elseif  line+array_param.DCM.rline_rng(end) > Nx
      DCM_rline_rng = array_param.DCM.rline_rng(1) : Nx-line;
    else
      DCM_rline_rng = array_param.DCM.rline_rng;
    end
  end
  
  %% Steering Vector setup
  for ml_idx = 1:length(array_param.imgs)
    % Make column vectors of y and z-positions
    for wf_adc_idx = 1:length(array_param.fcs{ml_idx})
      y_pos{ml_idx}(wf_adc_idx,1) = array_param.fcs{ml_idx}{wf_adc_idx}.pos(2,line);
      z_pos{ml_idx}(wf_adc_idx,1) = array_param.fcs{ml_idx}{wf_adc_idx}.pos(3,line);
    end
  end
  % Determine Steering Vector for beam-forming approaches
  for ml_idx = 1:length(array_param.imgs)
    % Make column vectors of y and z-positions
    [~,array_param.sv{ml_idx}] = array_param.sv_fh(array_param.Nsv,array_param.wfs.fc,y_pos{ml_idx},z_pos{ml_idx});
  end
  
  %% Setup line-varying parameters for wideband doa estimator
  if array_param.method == 7
    doa_param.y_pc  = y_pos{1};
    doa_param.z_pc  = z_pos{1};
    doa_param.theta = fftshift(theta);
    doa_param.SV    = fftshift(array_param.sv{1},2);
  end
  
  %% Setup line-varying parameters for wideband doa estimator
  if array_param.method == 8
    doa_param.y_pc   = y_pos{1};
    doa_param.z_pc   = z_pos{1};
    doa_param.theta  = fftshift(theta);
  end
  
  %% Setup line-varying parameters for wideband doa estimator
  if array_param.method == 9
    doa_param.y_pc  = y_pos{1};
    doa_param.z_pc  = z_pos{1};
    doa_param.theta = fftshift(theta);
    doa_param.SV    = fftshift(array_param.sv{1},2);
  end
  
  %% Setup array_proc_bin_idxs (these are the idxs into array_param.bins that will actually be evaluated)
  if ~restrict_bins
    % Do all the bins
    array_proc_bin_idxs = 1:length(array_param.bins);
  else
    % Restrict bins based on array_param.bin_restriction
    if array_param.method == 8 || array_param.method == 9
      tomo_top_bin = array_param.bin_restriction(1).bin(line) - ...
        min(array_param.bin_rng - floor((array_param.W - 1)/2));
      tomo_bot_bin = array_param.bin_restriction(2).bin(line) - ...
        max(array_param.bin_rng + floor((array_param.W-1)/2));
    else
      tomo_top_bin = array_param.bin_restriction(1).bin(line) - ...
        min(array_param.bin_rng);
      tomo_bot_bin = array_param.bin_restriction(2).bin(line) - ...
        max(array_param.bin_rng);
    end
    array_proc_bin_idxs = find((array_param.bins >= tomo_top_bin) & (array_param.bins <= tomo_bot_bin));
    if array_param.method == 7 || array_param.method == 8 || array_param.method == 9
      dout.tomo_top(lineIdx) = array_proc_bin_idxs(1);
      dout.tomo_bot(lineIdx) = array_proc_bin_idxs(end);
    end
  end
  
  %% Preallocate/Initialize Sarray
  switch array_param.method
    case {0,1,2,3,4,5,6}
      Sarray = zeros(Nsv,length(array_param.bins));
    case {7,8,9}
      Sarray = zeros(length(array_param.bins),Nsv);
  end
  
  %% DEBUG CODE FOR TESTING DOA CONSTRAINTS
  if 0 && lineIdx > 1
    warning off
    hist_bins = dout.tomo_top(line)+(150:700).';
    hist_poly = polyfit(hist_bins,dout.doa(hist_bins,lineIdx-1),2);
    warning on;
  end
  
  %% Iterate through each range bin and calculate third dimension (y-dimension) of image.
  for binIdx_idx = 1:1:length(array_proc_bin_idxs)
    binIdx = array_proc_bin_idxs(binIdx_idx);
    bin = array_param.bins(binIdx);
    
    % Handle the case when the data covariance matrix support pixels and
    % pixel neighborhood multilooking do not match. Note that this is only
    % supported for MVDR.
    if ~DCM_ML_match
      if bin+array_param.DCM.bin_rng(1) < 1
        DCM_bin_rng = 1-bin : array_param.DCM.bin_rng(end);
      elseif  bin+array_param.DCM.bin_rng(end) > Nt
        DCM_bin_rng = array_param.DCM.bin_rng(1) : Nt-bin;
      else
        DCM_bin_rng = array_param.DCM.bin_rng;
      end
    end
    
    %% Apply the array processing method to determine the current output pixel
    switch array_param.method
      case 0
        %% Periodogram
        dataSample = din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:);
        dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb, Nc]);
        Sarray(:,binIdx) = mean(abs(array_param.sv{1}(:,:)'*dataSample.').^2,2);
        for ml_idx = 2:length(din)
          dataSample = din{ml_idx}(bin+array_param.bin_rng,line+rline_rng,:,:,:);
          dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
          Sarray(:,binIdx) = Sarray(:,binIdx) ...
            + mean(abs(array_param.sv{ml_idx}(:,:)'*dataSample.').^2,2);
        end
        Sarray(:,binIdx) = Sarray(:,binIdx) / length(din);
        
      case 1
        %% MVDR
        if DCM_ML_match
          % The data covariance matrix creation uses the same set of snapshots/pixels
          % than the multilook does. Implement efficient algorithm.
          
          dataSample = double(din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:));
          dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
          Rxx = 1/size(dataSample,1) * (dataSample' * dataSample);
          %         imagesc(lp(Rxx))
          %         pause;
          diagonal = sqrt(mean(mean(abs(Rxx).^2))) * diag(ones(Nc,1),0);
          Rxx_inv = inv(Rxx + array_param.diag_load*diagonal);
          for freqIdx = 1:size(array_param.sv{1},2)
            Sarray(freqIdx,binIdx) = single(real(array_param.sv{1}(:,freqIdx).' * Rxx_inv * conj(array_param.sv{1}(:,freqIdx))));
          end
          for ml_idx = 2:length(din)
            dataSample = double(din{ml_idx}(bin+array_param.bin_rng,line+rline_rng,:,:,:));
            dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
            Rxx = 1/size(dataSample,1) * (dataSample' * dataSample);
            %         imagesc(lp(Rxx))
            %         pause;
            diagonal = sqrt(mean(mean(abs(Rxx).^2))) * diag(ones(Nc,1),0);
            Rxx_inv = inv(Rxx + array_param.diag_load*diagonal);
            for freqIdx = 1:size(array_param.sv{ml_idx},2)
              Sarray(freqIdx,binIdx) = Sarray(freqIdx,binIdx) ...
                + single(real(array_param.sv{ml_idx}(:,freqIdx).' * Rxx_inv * conj(array_param.sv{ml_idx}(:,freqIdx))));
            end
          end
          Sarray(:,binIdx) = 1 ./ (Sarray(:,binIdx) / length(din));
          
        else
          % The data covariance matrix creation uses a different set of snapshots/pixels
          % than the multilook does.
          dataSample = double(din{1}(bin+DCM_bin_rng,line+DCM_rline_rng,:,:,:));
          dataSample = reshape(dataSample,[length(DCM_bin_rng)*length(DCM_rline_rng)*Na*Nb Nc]).';
          Rxx = 1/size(dataSample,2) * (dataSample * dataSample');
          
          %         imagesc(lp(Rxx))
          %         pause;
          diagonal = sqrt(mean(mean(abs(Rxx).^2))) * diag(ones(Nc,1),0);
          Rxx_inv = inv(Rxx + array_param.diag_load*diagonal);
          for freqIdx = 1:size(array_param.sv{1},2)
            w = array_param.sv{1}(:,freqIdx)' * Rxx_inv;
            dataSample = double(din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:));
            dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]).';
            
            Sarray(freqIdx,binIdx) = mean(abs(w * dataSample).^2);
            Sarray(freqIdx,binIdx) = Sarray(freqIdx,binIdx) / abs(w * array_param.sv{1}(:,freqIdx)).^2;
          end
          for ml_idx = 2:length(din)
            dataSample = double(din{ml_idx}(bin+DCM_bin_rng,line+DCM_rline_rng,:,:,:));
            dataSample = reshape(dataSample,[length(DCM_bin_rng)*length(DCM_rline_rng)*Na*Nb Nc]).';
            Rxx = 1/size(dataSample,2) * (dataSample * dataSample');
            %         imagesc(lp(Rxx))
            %         pause;
            
            diagonal = sqrt(mean(mean(abs(Rxx).^2))) * diag(ones(Nc,1),0);
            Rxx_inv = inv(Rxx + array_param.diag_load*diagonal);
            for freqIdx = 1:size(array_param.sv{ml_idx},2)
              w = array_param.sv{ml_idx}(:,freqIdx)' * Rxx_inv;
              dataSample = double(din{ml_idx}(bin+array_param.bin_rng,line+rline_rng,:,:,:));
              dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]).';
              
              Sarray(freqIdx,binIdx) = Sarray(freqIdx,binIdx) ...
                + mean(abs(w * dataSample).^2) / abs(w * array_param.sv{ml_idx}(:,freqIdx)).^2;
            end
            Sarray(freqIdx,binIdx) = Sarray(freqIdx,binIdx) / size(array_param.sv{1},2);
          end
        end
        
      case 2
        %% MUSIC pseudospectrum
        %  The music algorithm fIdxs the eigenvectors of the correlation
        %  matrix. The inverse of the incoherent average of the magnitude
        %  squared spectrums of the smallest eigenvectors are used.
        dataSample = din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:);
        dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
        
        if isempty(array_param.sv)
          Sarray(:,binIdx) = pmusic(dataSample,array_param.Nsig,Nsv);
        else
          Rxx = 1/size(dataSample,1) * (dataSample' * dataSample);
          [V,D] = eig(Rxx);
          eigenVals = diag(D);
          [eigenVals noiseIdxs] = sort(eigenVals);
          
          % DEBUG CODE TO SLOWLY BUILD UP MUSIC SOLUTION, ONE EIGEN VECTOR
          % AT A TIME
          %           if 0
          %             if binIdx >162
          %               figure(1); clf;
          %               acc = 0;
          %               Nsig
          %               keyboard
          %               for sig_idx = 1:size(V,2)
          %                 acc = acc + abs(array_param.sv(:,:,lineIdx)'*V(:,sig_idx)).^2;
          %                 plot(fftshift(lp(1./acc)),'r')
          %                 plot(fftshift(lp(1./acc)))
          %                 hold on
          %               end
          %             end
          %             SmusicEV(:,binIdx) = eigenVals;
          %           end
          
          noiseIdxs = noiseIdxs(1:end-array_param.Nsig);
          Sarray(:,binIdx) = mean(abs(array_param.sv{1}(:,:).'*V(:,noiseIdxs)).^2,2);
        end
        for ml_idx = 2:length(din)
          dataSample = din{ml_idx}(bin+array_param.bin_rng,line+rline_rng,:,:,:);
          dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
          
          if isempty(array_param.sv)
            Sarray(:,binIdx) = pmusic(dataSample,array_param.Nsig,Nsv);
          else
            Rxx = 1/size(dataSample,1) * (dataSample' * dataSample);
            [V,D] = eig(Rxx);
            eigenVals = diag(D);
            [eigenVals noiseIdxs] = sort(eigenVals);
            
            noiseIdxs = noiseIdxs(1:end-array_param.Nsig);
            Sarray(:,binIdx) = Sarray(:,binIdx) ...
              + mean(abs(array_param.sv{ml_idx}(:,:).'*V(:,noiseIdxs)).^2,2);
          end
        end
        if isempty(array_param.sv)
          Sarray(:,binIdx) = Sarray(:,binIdx) / length(din);
        else
          Sarray(:,binIdx) = 0.5 ./ (Sarray(:,binIdx) / length(din));
        end
        
      case 3
        %% EIGENVECTOR pseudospectrum
        %  Same as MUSIC except the Idxividual noise subspace eigenvectors
        %  are weighted by the inverse of their corresponding eigenvalue
        %  when the incoherent averaging is done.
        error('Eigenvector not supported');
        dataSample = din(bin+eigmeth.bin_rng,line+rline_rng,:,:,:);
        dataSample = reshape(dataSample,[length(eigmeth.bin_rng)*length(rline_rng)*Na*Nb Nchan]);
        if uniformSampled
          Sarray(:,Idx) = peig(dataSample,Nsig,Ntheta);
        else
          Rxx = 1/size(dataSample,1) * (dataSample' * dataSample);
          [V,D] = eig(Rxx);
          eigenVals = diag(D).';
          [eigenVals noiseIdxs] = sort(eigenVals);
          noiseIdxs = noiseIdxs(1+Nsig:end);
          Sarray(:,Idx) = 1./mean(repmat(1./eigenVals,[size(sv,1) 1]).*abs(sv(:,:,lineIdx)'*V(:,noiseIdxs)).^2,2);
        end
        
      case 4
        %% RISR (re-iterative super resolution)
        %   See IEEE Transactions of Aerospace and Electronics Society
        %   2011?
        error('RISR not supported');
        dataSample = din(bin+array_param.bin_rng,line+rline_rng,:,:,:);
        dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
        
        M = size(array_param.sv,2);
        N = size(dataSample,2);
        array_param.num_iter = 15;
        array_param.alpha = 1;
        array_param.sigma_z = 0;
        array_param.R = diag([1 1 1 2 1 1 1]*1.25e-14);
        L = size(dataSample,2);
        
        dataS = dataSample.';
        W = array_param.sv(:,:);
        
        for iter = 1:array_param.num_iter
          %           fprintf('Iteration %d\n', iter);
          
          x_est = W'*dataS;
          SPD = zeros(M);
          for l = 1:L
            SPD = SPD + 1/L * (x_est(:,l) * x_est(:,l)');
          end
          SPD = SPD .* eye(M);
          
          if 0
            sig_est = sqrt(diag(SPD));
            if iter == 1
              figure(2); clf;
            else
              set(h_x_est,'Color','b');
              set(h_x_est,'Marker','none');
            end
            h_x_est = plot(lp(sig_est,2),'r.-');
            %ylim('manual');
            hold on;
            grid on;
            pause;
            %plot(lp(fftshift(fft(din,M)))); % DEBUG
          end
          
          AA = array_param.sv(:,:) * SPD * array_param.sv(:,:)';
          W = (AA + array_param.sigma_z*eye(N)*AA + array_param.alpha*array_param.R)^-1 * array_param.sv(:,:) * SPD;
        end
        sig_est = sqrt(diag(SPD));
        Sarray(:,binIdx) = sig_est;
        
      case 6
        %% Robust MVDR 2
        % Shahram Shahbazpanahi, Alex B. Gershman, Zhi-Quan Luo,
        % Kon Max Wong, ???Robust Adaptive Beamforming for General-Rank
        % Signal Models????, IEEE Transactions on Signal Processing, vol 51,
        % pages 2257-2269, Sept 2003
        %
        % Isaac Tan Implementation
        
        % The data covariance matrix creation uses a different set of snapshots/pixels
        % than the multilook does.
        dataSample = double(din{1}(bin+DCM_bin_rng,line+DCM_rline_rng,:,:,:));
        dataSample = reshape(dataSample,[length(DCM_bin_rng)*length(DCM_rline_rng)*Na*Nb Nc]).';
        Rxx = 1/size(dataSample,2) * (dataSample * dataSample');
        
        Rxx_inv = inv(Rxx);
        for freqIdx = 1:size(array_param.sv{1},2)
          sv_mat = array_param.sv{1}(:,freqIdx) * array_param.sv{1}(:,freqIdx)';
          sv_mat = sv_mat - 0.12*norm(sv_mat,'fro')*eye(size(sv_mat));
          
          % Get the biggest Eigenvector
          [eigen_vectors,eigen_values] = eig(Rxx*sv_mat);
          %           [~,max_idx] = max(real(diag(eigen_values)));
          [~,max_idx] = max(abs(diag(eigen_values)));
          w = eigen_vectors(:,max_idx)';
          dataSample = double(din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:));
          dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]).';
          
          Sarray(freqIdx,binIdx) = mean(abs(w * dataSample).^2);
          Sarray(freqIdx,binIdx) = Sarray(freqIdx,binIdx) / abs(w * array_param.sv{1}(:,freqIdx)).^2;
        end
        for ml_idx = 2:length(din)
          dataSample = double(din{ml_idx}(bin+DCM_bin_rng,line+DCM_rline_rng,:,:,:));
          dataSample = reshape(dataSample,[length(DCM_bin_rng)*length(DCM_rline_rng)*Na*Nb Nc]).';
          Rxx = 1/size(dataSample,2) * (dataSample * dataSample');
          
          Rxx_inv = inv(Rxx);
          for freqIdx = 1:size(array_param.sv{ml_idx},2)
            sv_mat = array_param.sv{ml_idx}(:,freqIdx) * array_param.sv{ml_idx}(:,freqIdx)';
            sv_mat = sv_mat - 0.25*norm(sv_mat,'fro')*eye(size(sv_mat));
            
            % Get the biggest Eigenvector
            [w,D] = eig(Rxx*sv_mat);
            w = w(:,1);
            dataSample = double(din{ml_idx}(bin+array_param.bin_rng,line+rline_rng,:,:,:));
            dataSample = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]).';
            
            Sarray(freqIdx,binIdx) = Sarray(freqIdx,binIdx) ...
              + mean(abs(w * dataSample).^2) / abs(w * array_param.sv{ml_idx}(:,freqIdx)).^2;
          end
          Sarray(freqIdx,binIdx) = Sarray(freqIdx,binIdx) / size(array_param.sv{1},2);
        end
        
      case 7
        %% MLE algorithm
        dataSample  = din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:);
        dataSample  = reshape(dataSample,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
        array_data  = dataSample.';
        Rxx         = (1/size(array_data,2)) * (array_data * array_data');
        doa_param.Rxx = Rxx; % put Rxx in doa_param (to pass to fiminsearch
        
        
        % Setup DOA Constraints
        for src_idx = 1:array_param.Nsig
          % Determine src_limits for each constraint
          doa_res = doa_param.doa_constraints(src_idx);
          switch (doa_res.method)
            case 'surfleft' % Incidence angle to surface clutter on left
              mid_doa(src_idx) = acos(array_param.fcs{1}{1}.surface(line) / array_param.wfs.time(bin));
            case 'surfright'% Incidence angle to surface clutter on right
              mid_doa(src_idx) = -acos(array_param.fcs{1}{1}.surface(line) / array_param.wfs.time(bin));
            case 'layerleft'
              table_doa   = [0:89.75]/180*pi;
              table_delay = array_param.fcs{1}{1}.surface(line) ./ cos(table_doa) ...
                + (doa_res.layer.twtt(line)-array_param.fcs{1}{1}.surface(line)) ./ cos(asin(sin(table_doa)/sqrt(er_ice)));
              doa_res.layer.twtt(line) = max(doa_res.layer.twtt(line),array_param.fcs{1}{1}.surface(line));
              if array_param.wfs.time(bin) <= doa_res.layer.twtt(line)
                mid_doa(src_idx) = 0;
              else
                mid_doa(src_idx) = interp1(table_delay, table_doa, array_param.wfs.time(bin));
              end
            case 'layerright'
              table_doa = [0:89.75]/180*pi;
              table_delay = array_param.fcs{1}{1}.surface(line) ./ cos(table_doa) ...
                + (doa_res.layer.twtt(line)-array_param.fcs{1}{1}.surface(line)) ./ cos(asin(sin(table_doa)/sqrt(er_ice)));
              doa_res.layer.twtt(line) = max(doa_res.layer.twtt(line),array_param.fcs{1}{1}.surface(line));
              if array_param.wfs.time(bin) <= doa_res.layer.twtt(line)
                mid_doa(src_idx) = 0;
              else
                mid_doa(src_idx) = -interp1(table_delay, table_doa, array_param.wfs.time(bin));
              end
            otherwise % 'fixed'
              mid_doa(src_idx) = 0;
          end
        end
        
        % Initialize search
        for src_idx = 1:array_param.Nsig
          doa_param.src_limits{src_idx} = mid_doa(src_idx) ...
            + doa_param.doa_constraints(src_idx).init_src_limits/180*pi;
          %doa_param.src_limits{src_idx} = [-pi/2 pi/2]; % DEBUG
        end
        
        doa_nonlcon_fh = eval(sprintf('@(x) doa_nonlcon(x,%f);', doa_param.theta_guard));
        
        theta0 = mle_initialization(Rxx,doa_param);
        
        %% Minimization of wb_cost_function
        % Set source limits
        LB = zeros(array_param.Nsig,1);
        UB = zeros(array_param.Nsig,1);
        for src_idx = 1:array_param.Nsig
          doa_param.src_limits{src_idx} = mid_doa(src_idx) ...
            + doa_param.doa_constraints(src_idx).src_limits/180*pi;
          %doa_param.src_limits{src_idx} = [-pi/2 pi/2]; % DEBUG
          LB(src_idx) = doa_param.src_limits{src_idx}(1);
          UB(src_idx) = doa_param.src_limits{src_idx}(2);
        end
        
        doa = [];
        
        warning off;
        [doa,Jval,exitflag,OUTPUT,~,~,HESSIAN] = ...
          fmincon(@(theta_hat) mle_cost_function(theta_hat,doa_param), theta0,[],[],[],[],LB,UB,doa_nonlcon_fh,doa_param.options);
        
        [doa,sort_idxs] = sort(doa);
        dout.doa(binIdx,:,lineIdx) = doa;
        dout.cost(binIdx,lineIdx) = Jval;
        dout.hessian(binIdx,:,lineIdx) = HESSIAN(sort_idxs + length(doa)*(sort_idxs-1));
        
        % Apply pseudoinverse and estimate power for each source
        k               = 4*pi*doa_param.fc/c;
        A               = exp(1i*k*(doa_param.y_pc*sin(doa(:)).' - doa_param.z_pc*cos(doa(:)).'));
        Weights         = inv(A'*A)*A';
        S_hat           = Weights*array_data;
        P_hat           = mean(abs(S_hat).^2,2);
        warning on;

        dout.power(binIdx,:,lineIdx)  = P_hat;
        
        if 0
          %% DEBUG code to plot cost function
          Ngrid     = 128;
          dNgrid    = 2/Ngrid;
          uy        = dNgrid*[0 : floor((Ngrid-1)/2), -floor(Ngrid/2) : -1];
          uz        = sqrt(1 - uy.^2);
          grid_vec  = atan2(uy,uz);
          grid_vec  = fftshift(grid_vec);
          switch array_param.Nsig
            case 1 % 1 source
              cf_vals = zeros(Ngrid,1);
              for eval_idx = 1:length(grid_vec);
                eval_theta = grid_vec(eval_idx);
                eval_theta = eval_theta(:);
                cf_vals(eval_idx) = mle_cost_function(eval_theta,doa_param);
              end
              figure(700);clf
              plot(grid_vec.*180/pi,cf_vals)
              grid on
            case 2  % 2 sources
              [grid1,grid2] = meshgrid(grid_vec);
              cf_vals = zeros(Ngrid,Ngrid);
              for row_index = 1:Ngrid
                for col_index = 1:Ngrid
                  eval_theta = [grid1(row_index,col_index) grid2(row_index,col_index)];
                  cf_vals(row_index,col_index) = mle_cost_function(eval_theta,doa_param);
                end
              end
              figure(701);clf
              grid_mask = grid1 <= grid2;
              cf_vals(grid_mask) = NaN;
              figure(101);mesh(grid1.*180/pi,grid2.*180/pi,-1.*cf_vals)
              xlabel('\theta_1')
              ylabel('\theta_2')
            otherwise
              error('Not supported')
          end
        end
        
        if 0 %% DEBUG code for bin restriction
          hist_bins = dout.tomo_top(line)+(150:700).';
          hist_poly = polyfit(hist_bins,dout.doa(hist_bins,lineIdx-1),2);
          plot(hist_bins,dout.doa(hist_bins,lineIdx-1),'.');
          hist_val = polyval(hist_poly,hist_bins);
          hold on;
          plot(hist_bins, hist_val,'r');
          hold off;
          
          hist_bins = dout.tomo_top(line)+(150:1700).';
          hist3([ hist_bins, dout.doa(hist_bins,lineIdx-1)],[round(length(hist_bins)/20) 30])
          set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        end
        
        
      case 8
        %% Parametric, space-time doa estimator for wideband or wide aperture direction of arrival estimation
        
        % Estimate space-time covariance matrix
        % ----------------------------------------------------------------
        dataSample = [];
        for W_offset = -floor((array_param.W-1)/2):floor((array_param.W-1)/2)
          offset_bin      = bin + W_offset;
          dataSample_tmp  = double(din{1}(offset_bin + array_param.bin_rng,line+rline_rng,:,:,:));
          dataSample_tmp  = reshape(dataSample_tmp,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]).';
          dataSample      = cat(1,dataSample,dataSample_tmp);
        end
        
        DCM             = (1/size(dataSample,2))*dataSample*dataSample';
        doa_param.DCM   = DCM;
        
        % Setup DOA Constraints
        for src_idx = 1:array_param.Nsig
          % Determine src_limits for each constraint
          doa_res = doa_param.doa_constraints(src_idx);
          switch (doa_res.method)
            case 'surfleft' % Incidence angle to surface clutter on left
              mid_doa(src_idx) = acos(array_param.fcs{1}{1}.surface(line) / array_param.wfs.time(bin));
            case 'surfright'% Incidence angle to surface clutter on right
              mid_doa(src_idx) = -acos(array_param.fcs{1}{1}.surface(line) / array_param.wfs.time(bin));
            case 'layerleft'
              table_doa   = [0:89.75]/180*pi;
              table_delay = array_param.fcs{1}{1}.surface(line) ./ cos(table_doa) ...
                + (doa_res.layer.twtt(line)-array_param.fcs{1}{1}.surface(line)) ./ cos(asin(sin(table_doa)/sqrt(er_ice)));
              doa_res.layer.twtt(line) = max(doa_res.layer.twtt(line),array_param.fcs{1}{1}.surface(line));
              if array_param.wfs.time(bin) <= doa_res.layer.twtt(line)
                mid_doa(src_idx) = 0;
              else
                mid_doa(src_idx) = interp1(table_delay, table_doa, array_param.wfs.time(bin));
              end
            case 'layerright'
              table_doa = [0:89.75]/180*pi;
              table_delay = array_param.fcs{1}{1}.surface(line) ./ cos(table_doa) ...
                + (doa_res.layer.twtt(line)-array_param.fcs{1}{1}.surface(line)) ./ cos(asin(sin(table_doa)/sqrt(er_ice)));
              doa_res.layer.twtt(line) = max(doa_res.layer.twtt(line),array_param.fcs{1}{1}.surface(line));
              if array_param.wfs.time(bin) <= doa_res.layer.twtt(line)
                mid_doa(src_idx) = 0;
              else
                mid_doa(src_idx) = -interp1(table_delay, table_doa, array_param.wfs.time(bin));
              end
            otherwise % 'fixed'
              mid_doa(src_idx) = 0;
          end
        end
        
        % Initialize search
        for src_idx = 1:array_param.Nsig
          doa_param.src_limits{src_idx} = mid_doa(src_idx) ...
            + doa_param.doa_constraints(src_idx).init_src_limits/180*pi;
          %doa_param.src_limits{src_idx} = [-pi/2 pi/2]; % DEBUG
        end
        theta0 = wb_initialization(DCM,doa_param);
        
        %% Minimization of wb_cost_function
        % Set source limits
        LB = zeros(array_param.Nsig,1);
        UB = zeros(array_param.Nsig,1);
        for src_idx = 1:array_param.Nsig
          doa_param.src_limits{src_idx} = mid_doa(src_idx) ...
            + doa_param.doa_constraints(src_idx).src_limits/180*pi;
          %doa_param.src_limits{src_idx} = [-pi/2 pi/2]; % DEBUG
          LB(src_idx) = doa_param.src_limits{src_idx}(1);
          UB(src_idx) = doa_param.src_limits{src_idx}(2);
        end
        
        % Transform intputs into constrained domain (fminsearch only)
        % for src_idx = 1:array_param.Nsig
        %   theta0(src_idx) = acos((theta0(src_idx) ...
        %     - sum(doa_param.src_limits{src_idx})/2)*2/diff(doa_param.src_limits{src_idx})); % fminsearch
        % end
        
        % Transform input out of constrained domain (fminsearch only)
        % for src_idx = 1:length(doa_param.src_limits)
        %   if ~isempty(doa_param.doa_constraints(src_idx).src_limits)
        %     doa(src_idx) = cos(doa(src_idx))*diff(doa_param.src_limits{src_idx})/2 ...
        %       + sum(doa_param.src_limits{src_idx})/2;
        %   end
        % end
        
        doa_nonlcon_fh = eval(sprintf('@(x) doa_nonlcon(x,%f);', doa_param.theta_guard));
        
        % Perform minimization
        %[doa,Jval,exitflag,OUTPUT] = ...
        %  fminsearch(@(theta_hat) wb_cost_function(theta_hat,doa_param), theta0,doa_param.options);
        [doa,Jval,exitflag,OUTPUT,~,~,HESSIAN] = ...
          fmincon(@(theta_hat) wb_cost_function(theta_hat,doa_param), theta0,[],[],[],[],LB,UB,doa_nonlcon_fh,doa_param.options);
        
        %% Estimate relative power for each source
        % -----------------------------------------------------------------
        % This section does the following:
        % 1) Creates an Nc x Nsnap complex valued matrix (where Nc is the
        % number of receive elements and Nsnap is the number of snapsots,
        % NOTE that Nsnap = length(array_param.bin_rng) +
        % length(array_param.line_rng) + Na + Nb), denoted by array_data
        %
        % 2) Computes the delays needed to registered to doas obtained for
        % a particular pixel and stores them in tau_reg,
        %
        % 3) Loops over sources and
        %       i.  Uses tau_reg to create sinc interpolation filters for
        %           each channel and tapers edges with Hamming weights,
        %       ii. Registers each channel to the particular source,
        %       iii.Applies pseudo-inverse to registered data,
        %       iv. Estimates signal and power,
        
        % Setup matrix of array data used to estimate relative power of
        % each source
        %
        %         array_data  = din{1}(bin+array_param.bin_rng,line+rline_rng,:,:,:);
        %         array_data  = reshape(array_data,[length(array_param.bin_rng)*length(rline_rng)*Na*Nb Nc]);
        %         array_data  = array_data.';
        
        if strcmpi('layerleft',doa_res.method)
          source_indexes = 1;
        elseif strcmpi('layerright',doa_res.method)
          source_indexes = 2;
        else
          source_indexes = 1:array_param.Nsig;
        end
        
        % S_hat: Nsig by bin_snapshots*rline_snapshots
        S_hat       = NaN*zeros(length(doa),length(array_param.bin_rng)*length(rline_rng));
        % tau_reg: longer delay to sensor means more negative
        tau_reg     = (2/c)*(doa_param.y_pc*sin(doa(:)).' - doa_param.z_pc*cos(doa(:)).');
        for src_idx = source_indexes;
          
          % Do each fast time snapshot separately and accumulate the result
          % into array_data
          array_data = [];
          for nt_idx = array_param.bin_rng
            offset_bin = bin + nt_idx;
            
            % For each channel
            registered_data = zeros(Nc,length(rline_rng)*Na*Nb);
            for nc_idx = 1:Nc
              tmp_chan_data = double(din{1}(offset_bin + array_param.reg_bins, line + rline_rng,:,:,nc_idx));
              tmp_chan_data = reshape(tmp_chan_data, [length(array_param.reg_bins) length(rline_rng)*Na*Nb]);
              % Create sinc interpolation coefficients with hamming window
              % tapir
              %   E.g. tau_reg negative means this sensor is delayed relative to
              %   the others and therefore the sinc peak should show up at > 0.
              Hinterp       = sinc(tau_reg(nc_idx,src_idx).*doa_param.fs + array_param.reg_bins.') ...
                .* hamming(length(array_param.reg_bins));

              % Apply the sinc filter to the data
              registered_data(nc_idx,:) ...
                = sum(repmat(Hinterp, [1 size(tmp_chan_data,2)]) .* tmp_chan_data);
            end
            
            % Accumulate the result for each fast time snapshot
            array_data = cat(2,array_data,registered_data);
          end
          A               = (exp(1i*2*pi*doa_param.fc*tau_reg));
          Weights         = inv(A'*A)*A';
          Weights         = Weights(src_idx,:);
          S_hat(src_idx,:)= Weights*array_data;
          
        end
        
        if 0
          H_interp = [];
          for nc_idx = 1:Nc
            H_interp(:,nc_idx) = sinc(tau_reg(nc_idx,src_idx).*doa_param.fs + array_param.reg_bins.') ...
              .*hamming(length(array_param.reg_bins));
          end
          
          figure;imagesc(lp(H_interp));
          figure;imagesc(lp(squeeze(din{1}(bin + array_param.reg_bins,line,1,1,:))))
          keyboard
        end
        
        
        if 0
          unreg_data = [];
          for debug_idx = array_param.bin_rng;
            tmp_unreg_data = double(din{1}(bin + debug_idx,line + rline_rng,:,:,:));
            tmp_unreg_data = reshape(tmp_unreg_data, [length(rline_rng)*Na*Nb Nc]).';
            unreg_data = cat(2,unreg_data,tmp_unreg_data);
          end
          
          
          haxis_list = [];
          h_fig = src_idx*100;
          figure(h_fig);clf;
          haxis_list(end+1) = axes('Parent', h_fig);
          hold(haxis_list(end),'on');
          plot(lp(interpft(unreg_data.',10*size(unreg_data,2))),'parent',haxis_list(end))
          grid on
          title(sprintf('Source %d, Before Registration',src_idx))
          
          h_fig = src_idx*100 + 1;
          figure(h_fig);clf;
          haxis_list(end+1) = axes('Parent', h_fig);
          hold(haxis_list(end),'on');
          plot(lp(interpft(array_data.',10*size(array_data,2))),'parent',haxis_list(end))
          grid on
          title(sprintf('Source %d, After Registration',src_idx))
          linkaxes(haxis_list, 'xy')
          
        end
        
        P_hat = mean(abs(S_hat).^2,2);
        
        % Collect outputs
        dout.doa(binIdx,:,lineIdx)      = doa;
        dout.cost(binIdx,lineIdx)       = Jval;
        dout.hessian(binIdx,:,lineIdx)  = diag(HESSIAN); % Not available with fminsearch
        dout.power(binIdx,:,lineIdx)    = P_hat;
        
        if 0
          %% DEBUG code for bin restriction
          hist_bins = dout.tomo_top(line)+(150:700).';
          hist_poly = polyfit(hist_bins,dout.doa(hist_bins,lineIdx-1),2);
          plot(hist_bins,dout.doa(hist_bins,lineIdx-1),'.');
          hist_val = polyval(hist_poly,hist_bins);
          hold on;
          plot(hist_bins, hist_val,'r');
          hold off;
          
          hist_bins = dout.tomo_top(line)+(150:1700).';
          hist3([ hist_bins, dout.doa(hist_bins,lineIdx-1)],[round(length(hist_bins)/20) 30])
          set(get(gca,'child'),'FaceColor','interp','CDataMode','auto');
        end
        
        if 0
          %% DEBUG code to plot cost function
          Ngrid     = 128;
          dNgrid    = 2/Ngrid;
          uy        = dNgrid*[0 : floor((Ngrid-1)/2), -floor(Ngrid/2) : -1];
          uz        = sqrt(1 - uy.^2);
          grid_vec  = atan2(uy,uz);
          grid_vec  = fftshift(grid_vec);
          switch array_param.Nsig
            case 1 % 1 source
              cf_vals = zeros(Ngrid,1);
              for eval_idx = 1:length(grid_vec);
                eval_theta = grid_vec(eval_idx);
                eval_theta = eval_theta(:);
                cf_vals(eval_idx) = wb_cost_function(eval_theta,doa_param);
              end
              figure(700);clf
              plot(grid_vec.*180/pi,cf_vals)
              grid on
            case 2  % 2 sources
              [grid1,grid2] = meshgrid(grid_vec);
              cf_vals = zeros(Ngrid,Ngrid);
              for row_index = 1:Ngrid
                for col_index = 1:Ngrid
                  eval_theta = [grid1(row_index,col_index) grid2(row_index,col_index)];
                  cf_vals(row_index,col_index) = wb_cost_function(eval_theta,doa_param);
                end
              end
              figure(701);clf
              grid_mask = grid1 <= grid2;
              cf_vals(grid_mask) = NaN;
              figure(101);mesh(grid1.*180/pi,grid2.*180/pi,-1*cf_vals)
              xlabel('\theta_1')
              ylabel('\theta_2')
            otherwise
              error('Not supported')
          end
        end
        
        
      case 9
        %% WBMLE Wideband Maximum Likelihood Estimator algorithm
        
        % Create data covariance matrix (DCM)
        Nsnap_td = length(array_param.bin_rng);
        Nsnap_other = length(rline_rng)*Na*Nb;
        NB = array_param.NB;
        DCM_fd = complex(zeros(Nc*NB,Nc));
        % Perform DFT for each set of data samples
        for idx = 1:NB:Nsnap_td-NB+1
          x_nb = fft(din{1}(offset_bin + array_param.bin_rng(idx + (0:NB-1)), ...
            line+rline_rng,:,:,:));
          for nb = 1:NB
            x_nb_snaps = reshape(x_nb(nb,:,:,:,:),[Nsnap_other Nc]);
            DCM_fd((nb-1)*Nc+(1:Nc),:) = DCM_fd((nb-1)*Nc+(1:Nc),:) + x_nb_snaps.'*conj(x_nb_snaps);
          end
        end
        DCM_fd = 1/(Nsnap_td*Nsnap_other) * DCM_fd;
        doa_param.DCM  = DCM_fd;
        
        % DOA Constraints
        for src_idx = 1:array_param.Nsig
          % Determine src_limits for each constraint
          doa_res = doa_param.doa_constraints(src_idx);
          switch (doa_res.method)
            case 'surfleft'
              mid_doa(src_idx) = acos(array_param.fcs{1}{1}.surface(line) / array_param.wfs.time(bin));
            case 'surfright'
              mid_doa(src_idx) = -acos(array_param.fcs{1}{1}.surface(line) / array_param.wfs.time(bin));
            case 'layerleft'
              table_doa = [0:89.75]/180*pi;
              table_delay = array_param.fcs{1}{1}.surface(line) ./ cos(table_doa) ...
                + (doa_res.layer.twtt(line)-array_param.fcs{1}{1}.surface(line)) ./ cos(asin(sin(table_doa)/sqrt(er_ice)));
              doa_res.layer.twtt(line) = max(doa_res.layer.twtt(line),array_param.fcs{1}{1}.surface(line));
              if array_param.wfs.time(bin) <= doa_res.layer.twtt(line)
                mid_doa(src_idx) = 0;
              else
                mid_doa(src_idx) = interp1(table_delay, table_doa, array_param.wfs.time(bin));
              end
            case 'layerright'
              table_doa = [0:89.75]/180*pi;
              table_delay = array_param.fcs{1}{1}.surface(line) ./ cos(table_doa) ...
                + (doa_res.layer.twtt(line)-array_param.fcs{1}{1}.surface(line)) ./ cos(asin(sin(table_doa)/sqrt(er_ice)));
              doa_res.layer.twtt(line) = max(doa_res.layer.twtt(line),array_param.fcs{1}{1}.surface(line));
              if array_param.wfs.time(bin) <= doa_res.layer.twtt(line)
                mid_doa(src_idx) = 0;
              else
                mid_doa(src_idx) = -interp1(table_delay, table_doa, array_param.wfs.time(bin));
              end
            otherwise % 'fixed'
              mid_doa(src_idx) = 0;
          end
        end
        
        % Initialize search
        for src_idx = 1:array_param.Nsig
          doa_param.src_limits{src_idx} = mid_doa(src_idx) ...
            + doa_param.doa_constraints(src_idx).init_src_limits/180*pi;
        end
        theta0 = wbmle_initialization(DCM,doa_param);
        
        %% Minimization of wb_cost_function
        % Set source limits
        LB = zeros(array_param.Nsig,1);
        UB = zeros(array_param.Nsig,1);
        for src_idx = 1:array_param.Nsig
          doa_param.src_limits{src_idx} = mid_doa(src_idx) ...
            + doa_param.doa_constraints(src_idx).src_limits/180*pi;
          LB(src_idx) = doa_param.src_limits{src_idx}(1);
          UB(src_idx) = doa_param.src_limits{src_idx}(2);
        end
        
        % Perform minimization
        [doa,Jval,exitflag,OUTPUT,~,~,HESSIAN] = ...
          fmincon(@(theta_hat) wbmle_cost_function(theta_hat,doa_param), theta0,[],[],[],[],LB,UB,[],doa_param.options);
        
        % Collect outputs
        dout.doa(binIdx,:,lineIdx)  = doa;
        dout.cost(binIdx,lineIdx)   = Jval;
        dout.func_counts(binIdx,lineIdx)  = OUTPUT.funcCount;
        dout.hessian(binIdx,:,lineIdx)  = diag(HESSIAN);
        
    end
  end
  
  if 0 %SEARCH
    scatterPlot = search_range(xaxis)*180/pi;
    figure(2);clf;
    %     plot(scatterPlot,'*')
    plot(scatterPlot,1:size(scatterPlot,1),'b*')
    %     plot(scatterPlot(650:850,:),650:850,'b*')
    %     for row = 1:size(scatterPlot,1)
    %       plot(scatterPlot(row,:),row,'*')
    %       hold on
    %     end
    set(gca,'ydir','reverse')
    %     axis([-90 90 0 size(scatterPlot,1)])
    %     axis([-90 90 950 1150])
  end
  
  %% Reformat output for this range line into a single slice of a 3D echogram
  switch array_param.method
    case {0,1,2,3,4,5,6}
      Sarray = fftshift(Sarray.',2);
    case {7,8,9}
      %% Not supported yet.
      
  end
  
  %% Reformat output to store full 3-D image (if enabled)
  if array_param.three_dim.en
    switch array_param.method
      case {0,1,2,3,4,5,6}
        dout.img(:,:,lineIdx) = Sarray;
        
      case {7,8,9}
        % Nothing to be done at this point
    end
  end
  
  %% Find maximum frequency (DOA) bin and the value
  % This pair provides the standard echogram (.val) and an echoram of the
  % DOA (.freq).
  % Search for maximum value in each range bin and keep track of
  % which cross-track frequency the max falls in
  [dout.val(:,lineIdx) dout.freq(:,lineIdx)] = max(Sarray(:,array_param.freq_rng),[],2);
  dout.freq(:,lineIdx) = array_param.freq_rng(dout.freq(:,lineIdx));
  
  if 0 && (~mod(lineIdx,size(dout.val,2)) || lineIdx == 1)
    % DEBUG CODE: change 0&& to 1&& on line above to run it
    switch array_param.method
      case {0,1,2,3,4,5,6}
        figure(1); clf;
        imagesc(10*log10(Sarray));
        hold on
        plot(dout.freq(:,lineIdx),1:size(dout.freq,1),'k')
        hold off
        
        figure(2); clf;
        imagesc(10*log10(dout.val));
        
      case {7,8,9}
        figure(1); clf;
        plot(dout.doa(:,:,lineIdx)*180/pi,'.')
        hold on;
        surf = interp1(array_param.wfs.time,1:length(array_param.wfs.time), ...
          array_param.fcs{1}{1}.surface);
        surf = interp1(array_param.bins, 1:length(array_param.bins), surf);
        plot(surf(line)*ones(1,2),[-90 90],'k')
        ylim([-90 90])
        surf_curve = acosd(array_param.fcs{1}{1}.surface(line) ./ array_param.wfs.time(array_param.bins));
        bad_mask = array_param.wfs.time(array_param.bins) < array_param.fcs{1}{1}.surface(line);
        surf_curve(bad_mask) = NaN;
        plot(surf_curve,'r')
        hold on
        plot(-1.*surf_curve,'r')
        
        if isfield(array_param.doa_constraints,'layer')
          table_doa = [0:89.75]/180*pi;
          table_delay = array_param.fcs{1}{1}.surface(line) ./ cos(table_doa) ...
            + (array_param.doa_constraints(2).layer.twtt(line)-array_param.fcs{1}{1}.surface(line)) ./ cos(asin(sin(table_doa)/sqrt(er_ice)));
          plot(interp1(array_param.wfs.time(array_param.bins), 1:length(array_param.wfs.time(array_param.bins)), ...
            table_delay), table_doa*180/pi, 'k');
          plot(interp1(array_param.wfs.time(array_param.bins), 1:length(array_param.wfs.time(array_param.bins)), ...
            table_delay), table_doa*180/pi+doa_param.doa_constraints(1).src_limits(1), 'k');
          plot(interp1(array_param.wfs.time(array_param.bins), 1:length(array_param.wfs.time(array_param.bins)), ...
            table_delay), table_doa*180/pi+doa_param.doa_constraints(1).src_limits(2), 'k');
        end
    end
    
    toc
    keyboard
  end
  
end

% Save theta variable
array_param.theta = theta_fftshift;

return;
