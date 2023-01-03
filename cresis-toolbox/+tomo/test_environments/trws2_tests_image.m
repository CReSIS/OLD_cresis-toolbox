MAX_LOOPS = 1;

% 'None' displays no indices
% 'ALL' displays index numbers for every cell
% 'FT' displays traversal order     for a fast-time-based search
% 'CT' displays traversal order for a cross-track-based search
INDEX_LABEL = 'None';
% Display traversal order for even loops of TRWS
INDEX_EVEN_LOOP = true;

% Display on this figure
FIGURE_NUM = 1;

% Downsampling interval
CROP = 50;

% Find surface of Greenland data
param.radar_name = 'rds';
param.radar.lever_arm_fh = '@lever_arm';
surfdata_source = 'surfData_paden';

param.season_name = '2014_Greenland_P3';
param.day_seg = '20140502_01';
out_type = 'multipass';
frm = 41;
echogram_fn = fullfile(ct_filename_out(param,out_type,''), 'summit_2012_2014_allwf_2012_music.mat');

if ~exist('params', 'var')
  params = read_param_xls(ct_filename_param('rds_param_2014_Greenland_P3.xls'), '20140502_01');
  params = ct_set_params(params,'cmd.generic',1);
end

layer_params.name = 'bottom';
layer_params.source = 'layerdata';
layer_params.layerdata_source = 'layer';
layer_params.existence_check = false;

if ~exist('mdata', 'var')
  mdata = load(echogram_fn);
end
trws_data = mdata.Tomo.img(1:CROP, 1:CROP, 1:CROP);
trws_data = trws_data - min(trws_data(:)) +.1;
trws_data = 10*log10(trws_data);

[Nt, Nsv, Nx] = size(trws_data);
Z = 1:Nt;
Y = 1:Nsv;
X = 1:Nx;

bounds = ones(2, Nx);
bounds(2, :) = Nt;

min_bounds = ones(Nt, Nx);
max_bounds = ones(Nt, Nx)*Nsv;
    

at_slope  = zeros(1, Nx);
at_weight = 1;
ct_slope  = zeros(Nsv, Nx);
ct_weight = ones(1, Nsv);

correct_surface = tomo.trws2(single(trws_data),single(at_slope), ...
  single(at_weight),single(ct_slope),single(ct_weight), ...
  uint32(MAX_LOOPS), uint32(bounds - 1), uint32(1), uint32(min_bounds - 1), uint32(max_bounds - 1));

figure(FIGURE_NUM);
clf;
hold on;


surf(X, Y, correct_surface, 'FaceAlpha', .2);
% surf(X, Y, repmat(bounds(1, :), Nsv, 1), 'FaceColor', [86, 135, 214]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);
% surf(X, Y, repmat(bounds(2, :), Nsv, 1), 'FaceColor', [232, 145, 90]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);

xlim([1 Nx]);
%     xticks(1:Nx);
xlabel('X : Along-Track (Nx, DIM 2)');
set(gca, 'XColor', 'b');
set(gca, 'xdir', 'reverse');

ylim([1 Nsv]);
%     yticks(1:Nsv);
ylabel('Y : Cross-Track (Nsv. DIM 1)');
set(gca, 'YColor', 'g');
set(gca, 'ydir', 'reverse');

zlim([0 Nt]);
%     zticks(0:Nt);
zlabel('Z : Fast-Time (Nt, DIM 0)');
set(gca, 'ZColor', 'r');
set(gca, 'zdir', 'reverse');

view([-45 90 90]);
cameratoolbar('SetCoordSys', 'none');
cameratoolbar('SetMode', 'nomode');
rotate3d on;

camva(10);
colormap(bone);

Nsvs_center = floor(Nsv/2);
Nsvs_array = [Nsvs_center:-1:1 (Nsvs_center+1):Nsv];

