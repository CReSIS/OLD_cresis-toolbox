function [data] = runOpsGetCrossoversWithinPolygon(varargin)
% script runOpsGetCrossoversWithinPolygon.m
%
% Example script for running opsGetCrossoversWithinPolygon.m.
%
% Authors: Reece Mathews

% SET THE SYSTEM
% sys = 'accum';
% 
% clear param;
% param.properties.location = 'arctic';
% param.properties.bound = 'POLYGON((-47.2967767172642 67.9471696002916,-59.09608185472462 67.00540191709096,-52.71050061854119 58.62701549679592,-35.59641918466397 58.04436438623961,-21.754770528971047 66.95314921567221,-22.58875553994411 67.09033936337771,-23.688058116710163 67.26160670251588,-47.2967767172642 67.9471696002916))';
% param.properties.seasons = {'2014_Greenland_P3', '2013_Greenland_P3'};

global gRadar;

switch nargin
  case 0
    sys = 'rds';
    
    clear param;
    param.properties.location = 'arctic';
    param.properties.bound = 'POLYGON((-47.6998416363594 69.29327309180026,-48.03263869201905 69.2388893651732,-47.99316494165699 69.13627210479136,-47.420511242020055 69.13617626790551,-47.6998416363594 69.29327309180026))';
    % POLYGON((-41.50736205994954 63.094240466556094,-41.53020306022716 62.8157382852153,-41.00806137448067 62.80233791769747,-40.98097889949436 63.07202824018716,-41.50736205994954 63.094240466556094))
    param.properties.seasons = {'2009_Greenland_TO', '2019_Greenland_P3', '2014_Greenland_P3'};
  case 2
    sys = varargin{1};
    param = varargin{2};
end

%% AUTOMATED SECTION BELOW

start = tic;
fprintf('Running opsGetCrossoversWithinPolygon ...\n');
[status,data] = opsGetCrossoversWithinPolygon(sys,param);

if status ~= 1
  error(data);
  return;
end

% url = ['http://ops.cresis.ku.edu/' data.properties];
url = [gRadar.ops.url data.properties];


pause(.5);

% Check to see if file is ready to download
opts = delimitedTextImportOptions('Delimiter', ',');

disp 'Waiting for file to be available to download. This may take a few minutes...';
while true
  try
    data = readmatrix(url, opts);
    break;
  catch ME
    pause(1);
  end
end

if size(data, 1) == 1
  warning("No crossovers found for selected seasons in selected polygon");
  return;
end

if startsWith(data{1}, 'Error')
  opts = weboptions("ContentType", "text");
  err = webread(url, opts);
  error(err);
end

stop = toc(start);
fprintf('   Completed in %f seconds\n',stop);

disp 'Data loaded into the data variable';