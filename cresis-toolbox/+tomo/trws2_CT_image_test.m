FIGURE_NUM = 6;
SAVE_PATH = 'C:/Users/mathe/Documents/MATLAB/TRWS_CT results/';
RNG_SEED = 1;

% Start index
St  = NaN;
Ssv = NaN;
Sx  = NaN;

% End index
Et  = NaN;
Esv = NaN;
Ex  = NaN;

% Num sampled points (resolutions)
Nt  = 100;
Nsv = 64;
Nx  = 100;

MAX_LOOPS = 10; % Number of iterations of TRWS to perform
INDEX_EVEN_LOOP = true; % Plot index traversal order for an even loop
PLOT_POINTS = false; % Plot data points as sphere's
PLOT_INDICES = false; % Plot indices on which dt is performed
PLOT_THRESHOLD = 7; % The minimum intensity to plot if PLOT_POINTS is true
RELOAD_DATA = false; % When false, only load trws_data if not already loaded
MAX_TICKS = 15; % Max number of axis tick labels
LOAD_DATA = true; % Use image data to plot intensities

%% Load Data
param.radar_name = 'rds';
surfdata_source = 'surfData_paden';
out_type = 'music3D_paden';
param.season_name = '2019_Antarctica_Ground';
param.day_seg = '20200107_01';
frm = 1;

if 1
    echogram_fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    sb_param.surfdata_fn = fullfile(ct_filename_out(param,surfdata_source,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    surfaces = tomo.surfdata.update_file(sb_param.surfdata_fn,sb_param.surfdata_fn,echogram_fn);
    surface_names = {surfaces.surf.name};
    top_surface = surfaces.surf(strcmp(surface_names, 'top'));
    bottom_surface = surfaces.surf(strcmp(surface_names, 'bottom'));
    
    FIND_SURF = true;
    LOAD_DATA = true;
elseif 0
    correct_surface = load('C:/Users/mathe/Documents/MATLAB/TRWS_CT results/perm Antarctica/30 iterations/entire_matrix_surf.mat');
    correct_surface = correct_surface.correct_surface;
    echogram_fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    FIND_SURF = false;
    LOAD_DATA = true;
%     Nt  = NaN;
%     Nsv = NaN;
%     Nx  = NaN;
    Esv = 64;
elseif 0
    param.season_name = '2014_Greenland_P3';
    param.day_seg = '';
    out_type = 'multipass';
    frm = 4;
    echogram_fn = fullfile(ct_filename_out(param,out_type,''), 'summit_2012_2014_allwf_2012_music.mat');
    FIND_SURF = true;
    LOAD_DATA = true;
end

if LOAD_DATA
    if ~exist('mdata', 'var') || RELOAD_DATA
      mdata = load(echogram_fn);
    end
    trws_data = 10*log10(mdata.Tomo.img);
end


%% Create length vars

% Set default start index
if ~exist('St', 'var') || isnan(St) || isempty(St)
    St  = 1;
end
if ~exist('Ssv', 'var') || isnan(Ssv) || isempty(Ssv)
    Ssv = 1;
end
if ~exist('Sx', 'var') || isnan(Sx) || isempty(Sx)
    Sx  = 1;
end

% Set default end indices
if ~exist('Et', 'var') || isnan(Et) || isempty(Et)
    if FIND_SURF
        Et = size(trws_data, 1);
    else
        Et = size(correct_surface, 1); 
    end
end
if ~exist('Esv', 'var') || isnan(Esv) || isempty(Esv)
    if ~FIND_SURF
        error('Cannot automatically determine Esv without loading echo_gram.');
    end
    Esv = size(trws_data, 2);
end
if ~exist('Ex', 'var') || isnan(Ex) || isempty(Ex)
    if FIND_SURF
        Ex = size(trws_data, 3);
    else
        Ex = size(correct_surface, 2);
    end
end

% Set default resolutions
if ~exist('Nt', 'var') || isnan(Nt) || isempty(Nt)
    Nt  = Et-St+1;
end
if ~exist('Nsv', 'var') || isnan(Nsv) || isempty(Nsv)
    Nsv = Esv-Ssv+1;
end
if ~exist('Nx', 'var') || isnan(Nx) || isempty(Nx)
    Nx  = Ex-Sx+1;
end

% Calculate step size
Tt  = round((Et-St)/Nt);
Tsv = round((Esv-Ssv)/Nsv);
Tx  = round((Ex-Sx)/Nx);

Tt  = max(Tt, 1);
Tsv = max(Tsv, 1);
Tx  = max(Tx, 1);

% Resample and window data to given start and end indices and resolutions
if LOAD_DATA
    trws_data = trws_data(St:Tt:Et, Ssv:Tsv:Esv, Sx:Tx:Ex);
    % Normalize data
    trws_data = trws_data - min(trws_data(:));
    trws_data = trws_data/max(trws_data(:)) * 10;
    
    min_bounds = interp1(mdata.Time, 1:length(mdata.Time), top_surface.y, 'nearest', 'extrap');
    max_bounds = interp1(mdata.Time, 1:length(mdata.Time), bottom_surface.y, 'nearest', 'extrap');
    min_bounds = min_bounds(Ssv:Tsv:Esv, Sx:Tx:Ex) - St + 1;
    max_bounds = max_bounds(Ssv:Tsv:Esv, Sx:Tx:Ex) - St + 1;
    min_bounds = interp1(St:Tt:Et, 1:length(St:Tt:Et), min_bounds, 'nearest');
    max_bounds = interp1(St:Tt:Et, 1:length(St:Tt:Et), max_bounds, 'nearest');
end

% Correct num sampled points to match step size constraint
Nt  = length(St:Tt:Et);
Nsv = length(Ssv:Tsv:Esv);
Nx  = length(Sx:Tx:Ex);

% Create default slopes, weights, and bounds
at_slope  = zeros(1, Nx);
at_weight = 1;
CT_bounds = ones(2, Nx);
CT_bounds(2, :) = Nsv;

if ~LOAD_DATA
    min_bounds = ones(Nsv, Nx)*1;
    max_bounds = ones(Nsv, Nx)*Nt;
end

% Perform TRWS and save surface
if FIND_SURF
  correct_surface = tomo.trws2_CT_perm(single(trws_data),single(at_slope), ...
    single(at_weight), uint32(MAX_LOOPS), uint32(CT_bounds-1), min_bounds-1, max_bounds-1);
  save([SAVE_PATH 'entire_matrix_surf.mat'], 'correct_surface');
else
 correct_surface = correct_surface(St:Tt:Et, Sx:Tx:Ex);
end

figure(FIGURE_NUM);
clf;
hold on;

Z = 1:Nt;
Y = 1:Nsv;
X = 1:Nx;

xlims = [1 Nx];
zlims = [0 Nsv];
ylims = [1 Nt];

clear surf;
surf(correct_surface, 'FaceAlpha', .6, 'FaceColor', 'interp', 'LineStyle', 'none');
blue   = [86 , 135, 214]./255;
orange = [232, 145, 90 ]./255;
surf(X, Z, repmat(CT_bounds(1, :), Nt, 1), 'FaceColor', blue, 'LineStyle', 'none', 'FaceAlpha', 0.1);
surf(X, Z, repmat(CT_bounds(2, :), Nt, 1), 'FaceColor', orange, 'LineStyle', 'none', 'FaceAlpha', 0.1);

xlim(xlims);
x_points = round(linspace(1, Nx, min(Nx, MAX_TICKS)));
xticks(x_points);
x_labels = (x_points-1)*Tx+Sx;
xticklabels(x_labels);
xlabel('X : Along-Track (Nx, DIM 2)');
set(gca, 'XColor', 'b');
set(gca, 'xdir', 'reverse');

zlim(zlims);
z_points = round(linspace(0, Nsv, min(Nsv, MAX_TICKS)));
zticks(z_points);
z_labels = (z_points-1)*Tsv+Ssv;
zticklabels(z_labels);
zlabel('Z : Cross-Track (Nsv. DIM 1)');
set(gca, 'ZColor', 'g');

ylim(ylims);
y_points = round(linspace(1, Nt, min(Nt, MAX_TICKS)));
yticks(y_points);
y_labels = (y_points-1)*Tt+St;
yticklabels(y_labels);
ylabel('Y : Fast-Time (Nt, DIM 0)');
set(gca, 'YColor', 'r');
set(gca, 'ydir', 'reverse');

view([90 45 90]);
camorbit(90, 0, 'data', [0 0 1]);
camorbit(90, 0, 'data', [1 0 0]);
cameratoolbar('SetCoordSys', 'y');
cameratoolbar('SetMode', 'orbit');

camva(10);
% colormap(bone);


Nsvs_center = floor(Nsv/2);
Nsvs_array = [Nsvs_center:-1:1 (Nsvs_center+1):Nsv];

ylim(zlims);
zlim(ylims);
minh = surf(min_bounds, 'FaceColor', '#304163', 'LineStyle', '-', 'FaceAlpha', 0.1);
rotate(minh, [1 0 0], 90, [0 Nt 0]);
set(minh, 'YData', Nt - get(minh, 'YData'));
set(minh, 'ZData', get(minh, 'ZData') - min(get(minh, 'ZData')) + 1);
maxh = surf(max_bounds, 'FaceColor', '#87763d', 'LineStyle', '-', 'FaceAlpha', 0.1);
rotate(maxh, [1 0 0], 90, [0 Nt 0]);
set(maxh, 'YData', Nt - get(maxh, 'YData'));
set(maxh, 'ZData', get(maxh, 'ZData') - min(get(maxh, 'ZData')) + 1);
ylim(ylims);
zlim(zlims);

for w_idx = 1:Nx
  if ~INDEX_EVEN_LOOP
      w = Nx-w_idx;
  else
      w = w_idx-1;
  end
  for h_idx = 1:Nsv
      h = Nsvs_array(h_idx)-1;

      for d_idx = 1:Nt
          d = d_idx-1;

          if PLOT_POINTS
            intensity = trws_data(d_idx, h_idx, w_idx);
            color = 'b';
            if intensity > 1
                color = 'r';
            end
            if intensity >= PLOT_THRESHOLD
                plot3(w_idx, d_idx, h_idx, sprintf('%s.', color), 'MarkerSize', intensity);
            end
          end

          if h_idx == 1 && PLOT_INDICES
              text(w_idx, d_idx, h_idx, sprintf('%d', d+w*Nt + 1));
          end
      end
  end
end
