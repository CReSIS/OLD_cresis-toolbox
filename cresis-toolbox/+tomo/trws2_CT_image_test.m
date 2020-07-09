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

% Num sampled points
Nt  = 100;
Nsv = 100;
Nx  = 100;

MAX_LOOPS = 2;
INDEX_EVEN_LOOP = true;
PLOT_POINTS = false;
PLOT_INDICES = false;
PLOT_THRESHOLD = 7;

param.radar_name = 'rds';
surfdata_source = 'surfData_paden';

if 0
    param.season_name = '2019_Antarctica_Ground';
    param.day_seg = '20200107_01';
    out_type = 'music3D_paden';
    frm = 1;
    echogram_fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
elseif 1
    param.season_name = '2014_Greenland_P3';
    param.day_seg = '';
    out_type = 'multipass';
    frm = 4;
    echogram_fn = fullfile(ct_filename_out(param,out_type,''), 'summit_2012_2014_allwf_2012_music.mat');
end

if ~exist('mdata', 'var')
  mdata = load(echogram_fn);
end
trws_data = 10*log10(mdata.Tomo.img);

if ~exist('St', 'var') || isnan(St) || isempty(St)
    St  = 1;
end
if ~exist('Ssv', 'var') || isnan(Ssv) || isempty(Ssv)
    Ssv = 1;
end
if ~exist('Sx', 'var') || isnan(Sx) || isempty(Sx)
    Sx  = 1;
end

if ~exist('Et', 'var') || isnan(Et) || isempty(Et)
    Et  = size(trws_data, 1);
end
if ~exist('Esv', 'var') || isnan(Esv) || isempty(Esv)
    Esv = size(trws_data, 2);
end
if ~exist('Ex', 'var') || isnan(Ex) || isempty(Ex)
    Ex  = size(trws_data, 3);
end

if ~exist('Nt', 'var') || isnan(Nt) || isempty(Nt)
    Nt  = Et;
end
if ~exist('Nsv', 'var') || isnan(Nsv) || isempty(Nsv)
    Nsv = Esv;
end
if ~exist('Nx', 'var') || isnan(Nx) || isempty(Nx)
    Nx  = Ex;
end

% Step size
Tt  = round((Et-St)/Nt);
Tsv = round((Esv-Ssv)/Nsv);
Tx  = round((Ex-Sx)/Nx);

Tt  = max(Tt, 1);
Tsv = max(Tsv, 1);
Tx  = max(Tx, 1);

trws_data = trws_data(St:Tt:Et, Ssv:Tsv:Esv, Sx:Tx:Ex);
Nt  = size(trws_data, 1);
Nsv = size(trws_data, 2);
Nx  = size(trws_data, 3);

trws_data = trws_data - min(trws_data(:));
trws_data = trws_data/max(trws_data(:)) * 10;

at_slope  = zeros(1, Nx);
at_weight = 1;
bounds = ones(2, Nx);
bounds(2, :) = Nsv;
correct_surface = tomo.trws2_CT_perm(single(trws_data),single(at_slope), ...
  single(at_weight), uint32(MAX_LOOPS), uint32(bounds-1));

save([SAVE_PATH 'entire_matrix_surf.mat'], 'correct_surface');

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
surf(X, Z, repmat(bounds(1, :), Nt, 1), 'FaceColor', blue, 'LineStyle', 'none', 'FaceAlpha', 0.1);
surf(X, Z, repmat(bounds(2, :), Nt, 1), 'FaceColor', orange, 'LineStyle', 'none', 'FaceAlpha', 0.1);
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

              msg_idx = NaN;

              if h_idx == 1
                  msg_idx = d+w*Nt + 1;
              end

              if ~isnan(msg_idx) && PLOT_INDICES
                  text(w_idx, d_idx, h_idx, sprintf('%d', msg_idx));
              end
          end
      end
  end
end
