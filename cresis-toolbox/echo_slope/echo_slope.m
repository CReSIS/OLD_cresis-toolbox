function [slope, slope_corr] = echo_slope(param, param_override)

%keyboard





%% General Setup
% =====================================================================

param = merge_structs(param, param_override);

% fprintf('=====================================================================\n');
% fprintf('%s: %s (%s)\n', mfilename, param.day_seg, datestr(now));
% fprintf('=====================================================================\n');




rows = param.echo_slope.rows;
cols = param.echo_slope.cols;
max_slope = param.echo_slope.max_slope;
min_slope = param.echo_slope.min_slope;
n = param.echo_slope.n;
sigma_factor = param.echo_slope.sigma_factor;

%sigma X and Y based on the rows and columns
sigmaX = 1/(15 * sigma_factor);
sigmaY = 1.5/(sigma_factor);


%
x = -cols/2:cols/2;
y = -rows/2:rows/2;

%array of theta angles 
theta = linspace(min_slope, max_slope, n);

%cell array of tiles 
tiles = {n};

%populate each tile with slope data 
for i = 1:n
  [X,Y] = meshgrid(x,y);
  a = 12 * (pi) * sqrt(sigmaX^2 * sigmaY^2);
  B = (((X.*cosd(-theta(i)) - Y*sind(-theta(i))).^2) ./ 2*sigmaX^2);
  C = (((X.*sind(-theta(i)) - Y*cosd(-theta(i))).^2) ./ 2*sigmaY^2);
  F = a*exp(-(B+C));
  F = F / sum(F(:));
  
  tile.array = F;
  tile.slope = theta(i);
  
  tiles{i} = tile;

%   figure(i);
%   imagesc(tiles{i}.array);
%   image_label = sprintf('Tile # %d, theta = %.3f', i, theta(i));
%   title(image_label);
  
  
end
param.echo_slope.tiles = tiles;

mdata = load('/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_post/CSARP_standard/20140508_01/Data_20140508_01_057.mat');

figure(101);
imagesc(lp(mdata.Data))

colormap(1-gray(256))

echo_slope_task(mdata, param);


% xcorr

% %% Input Checks: cmd
% % =====================================================================
% 
% %Remove frames that do not exist from param.cmd.frms list
% frames = frames_load(params);
% params.cmd.frms = frames_param_cmd_frms(params,frames);
% 
%  % Load the current frame
%     frm_str{frm_idx} = sprintf('%s_%03d',param.day_seg,frm);
%     data_fn = fullfile(in_fn_dir, sprintf('Data_%s.mat',frm_str{frm_idx}));
%     %if frm_idx == 1
%       mdata = load_L1B(data_fn);
%       frm_start(frm_idx) = 1;
%       frm_stop(frm_idx) = length(mdata.GPS_time);
%       dt = mdata.Time(2) - mdata.Time(1);





