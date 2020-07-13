FIGURE_NUM = 1;
SAVE_PATH = 'C:/Users/mathe/Documents/MATLAB/TRWS_CT results/';

% Start index
St  = 1;
Ssv = 1;
Sx  = 1;

% End index
Et  = NaN;
Esv = NaN;
Ex  = NaN;

% Num sampled points (resolutions)
Nt  = 10;
Nsv = 10;
Nx  = 10;

MAX_LOOPS = 2; % Number of iterations of TRWS to perform
INDEX_EVEN_LOOP = true; % Plot index traversal order for an even loop
PLOT_POINTS = false; % Plot data points as sphere's
PLOT_INDICES = false; % Plot indices on which dt is performed
PLOT_THRESHOLD = 7; % The minimum intensity to plot if PLOT_POINTS is true
RELOAD_DATA = false; % When false, only load trws_data if not already loaded


%% Load Data
param.radar_name = 'rds';
surfdata_source = 'surfData_paden';

if 1
    param.season_name = '2019_Antarctica_Ground';
    param.day_seg = '20200107_01';
    out_type = 'music3D_paden';
    frm = 1;
    echogram_fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
    FIND_SURF = true;
elseif 0
    correct_surface = load('C:/Users/mathe/Documents/MATLAB/TRWS_CT results/perm Antarctica/entire_matrix_surf.mat');
    correct_surface = correct_surface.correct_surface;
    FIND_SURF = false;
    Nt  = NaN;
    Nsv = NaN;
    Nx  = NaN;
    Esv = 64;
elseif 0
    param.season_name = '2014_Greenland_P3';
    param.day_seg = '';
    out_type = 'multipass';
    frm = 4;
    echogram_fn = fullfile(ct_filename_out(param,out_type,''), 'summit_2012_2014_allwf_2012_music.mat');
    FIND_SURF = true;
end

if FIND_SURF
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
    Nt  = Et;
end
if ~exist('Nsv', 'var') || isnan(Nsv) || isempty(Nsv)
    Nsv = Esv;
end
if ~exist('Nx', 'var') || isnan(Nx) || isempty(Nx)
    Nx  = Ex;
end

% Calculate step size
Tt  = round((Et-St)/Nt);
Tsv = round((Esv-Ssv)/Nsv);
Tx  = round((Ex-Sx)/Nx);

Tt  = max(Tt, 1);
Tsv = max(Tsv, 1);
Tx  = max(Tx, 1);

% Resample and window data to given start and end indices and resolutions
if FIND_SURF
    trws_data = trws_data(St:Tt:Et, Ssv:Tsv:Esv, Sx:Tx:Ex);
    % Normalize data
    trws_data = trws_data - min(trws_data(:));
    trws_data = trws_data/max(trws_data(:)) * 10;
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

top_bounds = ones(Nsv, Nx) * 5;
bottom_bounds = ones(Nsv, Nx) * Nt - 5;

% Perform TRWS and save surface
if FIND_SURF
  correct_surface = tomo.trws2_CT_perm(single(trws_data),single(at_slope), ...
    single(at_weight), uint32(MAX_LOOPS), uint32(CT_bounds-1), top_bounds, bottom_bounds);
  save([SAVE_PATH 'entire_matrix_surf.mat'], 'correct_surface');
end

figure(FIGURE_NUM);
clf;
hold on;

Z = 1:Nt;
Y = 1:Nsv;
X = 1:Nx;

clear surf;
surf(correct_surface, 'FaceAlpha', .6, 'FaceColor', 'interp', 'LineStyle', 'none');
blue   = [86 , 135, 214]./255;
orange = [232, 145, 90 ]./255;
surf(X, Z, repmat(CT_bounds(1, :), Nt, 1), 'FaceColor', blue, 'LineStyle', 'none', 'FaceAlpha', 0.1);
surf(X, Z, repmat(CT_bounds(2, :), Nt, 1), 'FaceColor', orange, 'LineStyle', 'none', 'FaceAlpha', 0.1);
xlim([1 Nx]);
%     xticks(1:Nx);
xlabel('X : Along-Track (Nx, DIM 2)');
set(gca, 'XColor', 'b');
set(gca, 'xdir', 'reverse');

zlim([0 Nsv]);
%     zticks(0:Nsv);
zlabel('Z : Cross-Track (Nsv. DIM 1)');
set(gca, 'ZColor', 'g');

ylim([1 Nt]);
%     yticks(1:Nt);
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



if PLOT_POINTS || PLOT_INDICES
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
                if intensity > PLOT_THRESHOLD
                    plot3(w_idx, d_idx, h_idx, sprintf('%s.', color), 'MarkerSize', intensity);
                end
              end

              if h_idx == 1 && PLOT_INDICES
                  text(w_idx, d_idx, h_idx, sprintf('%d', d+w*Nt + 1));
              end
          end
      end
  end
end
