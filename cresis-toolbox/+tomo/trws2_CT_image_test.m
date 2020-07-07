FIGURE_NUM = 1;

% Start index
St = 1;
Ssv = 1;
Sx = 1;

% End index
Et  = 3503;
Esv = 64;
Ex  = 5073;

% Num sampled points
Nt  = 100;
Nsv = 64;
Nx  = 100;

% Step size
Tt = round((Et-St)/Nt);
Tsv = round((Esv-Ssv)/Nsv);
Tx = round((Ex-Sx)/Nx);

MAX_LOOPS = 2;
INDEX_EVEN_LOOP = true;
PLOT_POINTS = false;
PLOT_INDICES = false;
PLOT_THRESHOLD = 7;


param.radar_name = 'rds';
param.season_name = '2019_Antarctica_Ground';
out_type = 'music3D_paden';
surfdata_source = 'surfData_paden';
param.day_seg = '20200107_01';
frm = 1;
geotiff_fn = ct_filename_gis(param,fullfile('antarctica','Landsat-7','Antarctica_LIMA_480m.tif'));
ice_mask_fn = '';
bounds_relative = [0 0 0 0];

echogram_fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
if ~exist('mdata')
  mdata = load(echogram_fn);
end
trws_data = 10*log10(mdata.Tomo.img);
trws_data = trws_data(St:Tt:Et, Ssv:Tsv:Esv, Sx:Tx:Ex);
trws_data = trws_data - min(trws_data(:));
trws_data = trws_data/max(trws_data(:)) * 10;
save('/users/reece/Desktop/entire_matrix.mat', 'trws_data');

at_slope  = zeros(1, Nx);
at_weight = 1;
bounds = ones(2, Nx);
bounds(2, :) = Nsv;
correct_surface = tomo.trws2_CT_perm(single(trws_data),single(at_slope), ...
  single(at_weight), uint32(MAX_LOOPS), uint32(bounds-1));

save('/users/reece/Desktop/entire_matrix_surf.mat', 'correct_surface');

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
