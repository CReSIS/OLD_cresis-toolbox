function lidar = read_lidar_lvis(fns,varargin)
% Read LVIS Lidar Data. (LGE and TXT format)
% 
% lidar = read_lidar_lvis(fns,samp_type,samp_rate);
% 
%   fns = cell array of string filenames
%
%   samp_type = Optional type of sampling.
%
%     Options:
%       (1) 'factor' : Keeps every nth point
%       (2) 'track'  : Keeps a point every n meters
%
%   samp_rate = Optional rate of sampling.
%
%    (1) If samp_type = 'factor'. Specify "n" to keep every nth point.
%    (2) If samp_type = 'track'. Specify "n" to keep a point ever n meters.
%
% lidar = Output structure containing all data from LVIS file.
% 
% Disclaimers:
%   - This tool is written to support reading of LVIS data served from
%   NSIDC at: ftp://n4ftl01u.ecs.nasa.gov/SAN2/ICEBRIDGE_FTP/ILVIS2_LVISxyz_v01/
%   - This tool relies on two custom MATLAB functions written by CReSIS
%   staff at the University of Kansas.
%       (1) geodetic_to_along_track.m
%       (2) get_equal_alongtrack_spacing_idxs.m
%   -This tool is given with no implied guarentee of accuracy or
%   completeness. 
%
% Additional Information:
%   (1) The tool can only read one format (LGE or TXT) at a time.
%   (2) Downsampling LVIS data creates a data product NOT supported by LVIS.
%
% Authors:
%   (1) Kyle W. Purdon, CReSIS UGRA
%   (2) John King, CReSIS UGRA
%
% See also geodetic_to_along_track.m get_equal_alongtrack_spacing_idxs.m

% Check Function Input
if nargin > 3
  error('Too many inputs to function read_lidar_lvis.');
end

% Check for varargin or set defaults
if nargin > 1
  samp = true;
  samp_type = varargin{1};
  samp_rate = varargin{2};
else
  samp = false;
end

% Check samp type input
if samp
  if ~strcmp(samp_type,{'factor','track'})
    error('Invalid input: Samling Type can be "factor" or "track".');
  end
end

% Get all filenames and filetypes
name = cell(length(fns),1);
ext = cell(length(fns),1);
for f_idx = 1:length(fns)
  [tmp name{f_idx} ext{f_idx}] = fileparts(fns{f_idx});
end

% If any filetypes except LGE/TXT exist error out.
if sum((strcmpi(ext,'.lge'))+(strcmpi(ext,'.TXT'))) ~= length(ext)
  error('Only ".lge" OR ".TXT" are supported.');
end

% If both LGE and TXT exist error out.
if sum(strcmpi(ext,'.TXT')) ~= length(ext) && sum(strcmpi(ext,'.lge')) ~= length(ext)
  error('Can only process one format (TXT or LGE) at a time.');
end

% Initiate holder for output structure.
if strcmpi(ext{f_idx},'.TXT')
  lidar = struct('FID',[],'shotnumber',[],'time',[],'lon',[],'lat',[],'elev',[],'lonLow',[],...
    'latLow',[],'elevLow',[],'lonHigh',[],'latHigh',[],'elevHigh',[]);
elseif strcmpi(ext{f_idx},'.lge')
  lidar = struct('lfid',[],'shotnumber',[],'time',[],'glon',[],'glat',[],'zg',[],'rh25',[],...
    'rh50',[],'rh75',[],'rh100',[]);
end

% Read LVIS data (TXT or LGE)
for f_idx = 1:length(fns)
  % Read TXT Data
  if strcmpi(ext{f_idx},'.TXT')
    tic;
    fprintf('Opening %s ...\n',name{f_idx});
    fid = fopen(fns{f_idx},'r');
    fprintf('Reading contents (May take awhile)\n');
    data = textscan(fid,'%u %f %f %f %f %f %f %f %f %f %f %f','headerlines',2);
    fclose(fid);
    % Store the data into an array (Append for each file)
    lidar.FID = cat(1,lidar.FID,data{1});
    lidar.shotnumber = cat(1,lidar.shotnumber,data{2});
    lidar.time = cat(1,lidar.time,data{3});
    lidar.lon = cat(1,lidar.lon,data{4});
    lidar.lat = cat(1,lidar.lat,data{5});
    lidar.elev = cat(1,lidar.elev,data{6});
    lidar.lonLow = cat(1,lidar.lonLow,data{7});
    lidar.latLow = cat(1,lidar.latLow,data{8});
    lidar.elevLow = cat(1,lidar.elevLow,data{9});
    lidar.lonHigh = cat(1,lidar.lonHigh,data{10});
    lidar.latHigh = cat(1,lidar.latHigh,data{11});
    lidar.elevHigh = cat(1,lidar.elevHigh,data{12});
    fprintf('  Done (%.1f sec)\n', toc);
  elseif strcmpi(ext,'.lge')
    recordType = {'ulong','ulong','double','double','double','float','float','float','float','float'};
    recordLen = [4 4 8 8 8 4 4 4 4 4];
    data = cell(1,numel(recordType));
    tic;
    fprintf('Opening %s ...\n',name{f_idx});
    fid = fopen(fns{f_idx},'rb');
    fprintf('  Done (%.1f sec)\n', toc);
    fprintf('Reading contents (May take awhile)\n');
    for b_idx = 1:numel(recordType);
      fseek(fid, sum(recordLen(1:b_idx-1)),'bof');
      data{b_idx} = fread(fid, Inf, [recordType{b_idx}], sum(recordLen)-recordLen(b_idx),'ieee-be');
    end
    fclose(fid);
    % Store the data into an array (Append for each file)
    lidar.lfid = cat(1,lidar.lfid,data{1});
    lidar.shotnumber = cat(1,lidar.shotnumber,data{2});
    lidar.time = cat(1,lidar.time,data{3});
    lidar.glon = cat(1,lidar.glon,data{4});
    lidar.glat = cat(1,lidar.glat,data{5});
    lidar.zg = cat(1,lidar.zg,data{6});
    lidar.rh25 = cat(1,lidar.rh25,data{7});
    lidar.rh50 = cat(1,lidar.rh50,data{8});
    lidar.rh75 = cat(1,lidar.rh75,data{9});
    lidar.rh100 = cat(1,lidar.rh100,data{10});
    fprintf('  Done (%.1f sec)\n', toc);
  end
end

% Optionally sample data (factor or track method)
if samp
  if strcmpi(ext{f_idx},'.TXT')
    if strcmp(samp_type,'factor')
      tic;
      fprintf('Downsampling Contents (Factor Method)\n');
      lidar.FID = downsample(lidar.FID,samp_rate);
      lidar.shotnumber = downsample(lidar.shotnumber,samp_rate);
      lidar.time = downsample(lidar.time,samp_rate);
      lidar.lon = downsample(lidar.lon,samp_rate);
      lidar.lat = downsample(lidar.lat,samp_rate);
      lidar.elev = downsample(lidar.elev,samp_rate);
      lidar.lonLow = downsample(lidar.lonLow,samp_rate);
      lidar.latLow = downsample(lidar.latLow,samp_rate);
      lidar.elevLow = downsample(lidar.elevLow,samp_rate);
      lidar.lonHigh = downsample(lidar.lonHigh,samp_rate);
      lidar.latHigh = downsample(lidar.latHigh,samp_rate);
      lidar.elevHigh = downsample(lidar.elevHigh,samp_rate);
      fprintf('  Done (%.1f sec)\n', toc);
    elseif strcmp(samp_type,'track')
      tic;
      fprintf('Downsampling Contents (Along-Track Method)\n');
      along_track = geodetic_to_along_track(lidar.lat,lidar.lon,lidar.elev);
      decim_idxs  = get_equal_alongtrack_spacing_idxs(along_track,samp_rate);
      lidar.FID = lidar.FID(decim_idxs);
      lidar.shotnumber = lidar.shotnumber(decim_idxs);
      lidar.time = lidar.time(decim_idxs);
      lidar.lon = lidar.lon(decim_idxs);
      lidar.lat = lidar.lat(decim_idxs);
      lidar.elev = lidar.elev(decim_idxs);
      lidar.lonLow = lidar.lonLow(decim_idxs);
      lidar.latLow = lidar.latLow(decim_idxs);
      lidar.elevLow = lidar.elevLow(decim_idxs);
      lidar.lonHigh = lidar.lonHigh(decim_idxs);
      lidar.latHigh = lidar.latHigh(decim_idxs);
      lidar.elevHigh = lidar.elevHigh(decim_idxs);
      fprintf('  Done (%.1f sec)\n', toc);
    end
  elseif strcmpi(ext,'.lge')
    if strcmp(samp_type,'factor')
      tic;
      fprintf('Downsampling Contents (Factor Method)\n');
      lidar.lfid = downsample(lidar.lfid,samp_rate);
      lidar.shotnumber = downsample(lidar.shotnumber,samp_rate);
      lidar.time = downsample(lidar.time,samp_rate);
      lidar.glon = downsample(lidar.glon,samp_rate);
      lidar.glat = downsample(lidar.glat,samp_rate);
      lidar.zg = downsample(lidar.zg,samp_rate);
      lidar.rh25 = downsample(lidar.rh25,samp_rate);
      lidar.rh50 = downsample(lidar.rh50,samp_rate);
      lidar.rh75 = downsample(lidar.rh75,samp_rate);
      lidar.rh100 = downsample(lidar.rh100,samp_rate);
      fprintf('  Done (%.1f sec)\n', toc);
    elseif strcmp(samp_type,'track')
      tic;
      fprintf('Downsampling Contents (Along-Track Method)\n');
      along_track = geodetic_to_along_track(lidar.glat,lidar.glon,lidar.zg);
      decim_idxs  = get_equal_alongtrack_spacing_idxs(along_track,samp_rate);
      lidar.lfid = lidar.lfid(decim_idxs);
      lidar.shotnumber = lidar.shotnumber(decim_idxs);
      lidar.time = lidar.time(decim_idxs);
      lidar.glon = lidar.glon(decim_idxs);
      lidar.glat = lidar.glat(decim_idxs);
      lidar.zg = lidar.zg(decim_idxs);
      lidar.rh25 = lidar.rh25(decim_idxs);
      lidar.rh50 = lidar.rh50(decim_idxs);
      lidar.rh75 = lidar.rh75(decim_idxs);
      lidar.rh100 = lidar.rh100(decim_idxs);
      fprintf('  Done (%.1f sec)\n', toc);
    end
  end
end
return;