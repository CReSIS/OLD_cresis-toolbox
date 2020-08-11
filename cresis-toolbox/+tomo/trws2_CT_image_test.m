%% CONSTANTS
FIGURE_NUM = NaN;
SAVE_PATH = 'C:/Users/mathe/Documents/MATLAB/TRWS_CT results/perm Antarctica/output images/';
SAVE_IMAGE = false; % Save the normalized (output) image
SAVE_SURF = true; % Save the output surface
SAVE_FIG = true; % Save the created figure
SAVE_NAME = 'forward send 2 perms 1 pass full vert'; % name of save files and directory
RNG_SEED = 1;
TOP_LAYER_GUARD = 0;
BOTTOM_LAYER_GUARD = 0;
NOISE_FLOOR = -40;

% Start index
St  = NaN;
Ssv = NaN;
Sx  = NaN;

% End index
Et  = NaN;
Esv = NaN;
Ex  = NaN;

% Num sampled points (resolutions)
Nt  = NaN;
Nsv = NaN;
Nx  = NaN;

MAX_LOOPS = 1; % Number of iterations of TRWS to perform

PLOT_INDICES = false; % Plot indices on which dt is performed
INDEX_LOOP_NUM = 0; % Plot index traversal order for this loop num
MAX_TICKS = 10; % Max number of axis tick labels

PLOT_POINTS = true; % Plot data points as sphere's if PLOT_MAX_POINTS is false
PLOT_THRESHOLD = 'auto'; % The minimum intensity to plot if PLOT_POINTS is true. 'auto' selects the greatest AUTO_THRESHOLD_MAX_POINTS points;
AUTO_THRESHOLD_MAX_POINTS = 2000; % Number of points to plot when PLOT_THRESHOLD is set to 'auto'

PLOT_MAX_POINTS = false; % Plot maximum data point in the CT dimension for each (AT, FT) row. Overrides PLOT_POINTS
MAX_POINTS_NUM = 3; % Number of max points to plot for every (AT, FT) row.

PLOT_MAX_SURF = false; % Plot the maximum points as a surface. (LOAD_DATA must be true. Does not plot found surface when this is true)

LOAD_DATA = true; % Use image data to plot intensities
RELOAD_DATA = false; % When false, only load trws_data if not already loaded. Ignored if LOAD_DATA is false
FIND_SURF = true; % Calls TRWS when true. Overridden by presets below. correct_surface must be set manually when false.
COLOR_SURF_DATA = 'trws_data'; % Color the surface with intensity values from given variable. generally 'trws_data_norm' for output image or 'trws_data' for input image. 'none' for solid coloring;

PLOT_CT_BOUNDS = false; % Plot the CT-slope surface boundaries
PLOT_FT_BOUNDS = true; % Plot the FT-slope surface boundaries

NORMALIZE_DATA = true; % Normalize data for display identically to the normalization used by trws2_CT_perm.m. LOAD_DATA must be true.
USE_DEBUG_MATRIX = true; % Get the final image values from the TRWS2 algorithm and use these to plot points, etc. Replaces normalized data with debug data.
PLOT_HISTOGRAM = false; % Plot a histogram of the data

%% DISPLAY VARS
DISPLAY_VARS = {'RNG_SEED', 'TOP_LAYER_GUARD', 'BOTTOM_LAYER_GUARD', 'NOISE_FLOOR', 'St', 'Ssv', 'Sx', 'Et', 'Esv', 'Ex', 'Nt', 'Nsv', 'Nx', 'MAX_LOOPS', 'NORMALIZE_DATA', 'USE_DEBUG_MATRIX', 'FIND_SURF', 'COLOR_SURF_DATA', 'SAVE_IMAGE', 'SAVE_SURF', 'SAVE_FIG', 'SAVE_NAME'};
if PLOT_POINTS
  DISPLAY_VARS{end + 1} = 'PLOT_POINTS';
  DISPLAY_VARS{end + 1} = 'PLOT_THRESHOLD';
  if strcmp(PLOT_THRESHOLD, 'auto')
    DISPLAY_VARS{end + 1} = 'AUTO_THRESHOLD_MAX_POINTS';
  end
end
if PLOT_MAX_POINTS
  DISPLAY_VARS{end + 1} = 'PLOT_MAX_POINTS';
  DISPLAY_VARS{end + 1} = 'MAX_POINTS_NUM';
end

%% DATA VARS
param.radar_name = 'rds';
surfdata_source = 'surfData_paden';
out_type = 'music3D_paden';
param.season_name = '2019_Antarctica_Ground';
param.day_seg = '20200107_01';
frm = 1;

%% PRESETS
if 1
  % Find surface of antarctica data
  echogram_fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  
  FIND_SURF = true;
  LOAD_DATA = true;
elseif 0
  % Load surface from file
  correct_surface = load('C:/Users/mathe/Documents/MATLAB/TRWS_CT results/perm Antarctica/bounded/cropped surf traverse entire/no guard/1 iter/entire_matrix_surf.mat');
  correct_surface = correct_surface.correct_surface;
  echogram_fn = fullfile(ct_filename_out(param,out_type,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
  FIND_SURF = false;
  LOAD_DATA = true;
  %     Nt  = NaN;
  %     Nsv = NaN;
  %     Nx  = NaN;
  SAVE_FIG = false;
  SAVE_SURF = false;
  SAVE_IMAGE = false;
  Esv = 64;
elseif 0
  % Find surface of Greenland data
  param.season_name = '2014_Greenland_P3';
  param.day_seg = '';
  out_type = 'multipass';
  frm = 4;
  echogram_fn = fullfile(ct_filename_out(param,out_type,''), 'summit_2012_2014_allwf_2012_music.mat');
  FIND_SURF = true;
  LOAD_DATA = true;
end

%% LOAD SURFACES
surfdata_fn = fullfile(ct_filename_out(param,surfdata_source,''),sprintf('Data_%s_%03d.mat',param.day_seg,frm));
surfaces = tomo.surfdata.update_file(surfdata_fn,surfdata_fn,echogram_fn);
surface_names = {surfaces.surf.name};
top_surface = surfaces.surf(strcmp(surface_names, 'top'));
bottom_surface = surfaces.surf(strcmp(surface_names, 'bottom'));

%% LOAD IMAGE DATA
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

%% RESAMPLE DATA
% Resample and window data to given start and end indices and resolutions
if LOAD_DATA
  trws_data = trws_data(St:Tt:Et, Ssv:Tsv:Esv, Sx:Tx:Ex);
  trws_data_norm = trws_data;
  
  % Position FT-slope surface bounds
  min_bounds = interp1(mdata.Time, 1:length(mdata.Time), top_surface.y, 'nearest', 'extrap')+TOP_LAYER_GUARD;
  max_bounds = interp1(mdata.Time, 1:length(mdata.Time), bottom_surface.y, 'nearest', 'extrap')-BOTTOM_LAYER_GUARD;
  min_bounds = min_bounds(Ssv:Tsv:Esv, Sx:Tx:Ex) - St + 1;
  max_bounds = max_bounds(Ssv:Tsv:Esv, Sx:Tx:Ex) - St + 1;
  min_bounds = interp1(St:Tt:Et, 1:length(St:Tt:Et), min_bounds, 'nearest');
  max_bounds = interp1(St:Tt:Et, 1:length(St:Tt:Et), max_bounds, 'nearest');
  min_bound = min(min_bounds(:));
  max_bound = max(max_bounds(:));
end

% Correct num sampled points to match step size constraint
Nt  = length(St:Tt:Et);
Nsv = length(Ssv:Tsv:Esv);
Nx  = length(Sx:Tx:Ex);

if ~LOAD_DATA
  min_bounds = ones(Nsv, Nx)*1;
  max_bounds = ones(Nsv, Nx)*Nt;
end

%% PERFORM TRWS
% Create default slopes, weights, and bounds
at_slope  = zeros(1, Nx);
at_weight = 1;
CT_bounds = ones(2, Nx);
CT_bounds(2, :) = Nsv;

% Perform TRWS and save surface
debug = [];
if FIND_SURF
  if USE_DEBUG_MATRIX
    [correct_surface, debug] = tomo.trws2_CT_perm(single(trws_data),single(at_slope), ...
      single(at_weight), uint32(MAX_LOOPS), uint32(CT_bounds-1), min_bounds-1, max_bounds-1);
    trws_data_norm = max(debug(:)) - debug + .01;
  else
    correct_surface = tomo.trws2_CT_perm(single(trws_data),single(at_slope), ...
      single(at_weight), uint32(MAX_LOOPS), uint32(CT_bounds-1), min_bounds-1, max_bounds-1);
  end
else
  correct_surface = correct_surface(St:Tt:Et, Sx:Tx:Ex);
end


%% NORMALIZE IMAGE DATA
if LOAD_DATA
  % Normalize data
  if NORMALIZE_DATA && ~USE_DEBUG_MATRIX
    trws_data_norm = echo_norm(trws_data_norm, struct('scale', [NOISE_FLOOR 90]));
    for rline = 1:Nx
      for doa_bin = 1:Nsv
        trws_data_norm(1:Nt < min_bounds(doa_bin, rline), doa_bin, rline) = NOISE_FLOOR;
        trws_data_norm(1:Nt > max_bounds(doa_bin, rline), doa_bin, rline) = NOISE_FLOOR;
      end
    end 
    trws_data_norm = [nan(min_bound - 1, Nsv, Nx); trws_data_norm(min_bound:max_bound, :, :); nan(size(trws_data_norm, 1) - max_bound, Nsv, Nx)];
  end
  
  min_intensity = min(trws_data_norm(:));
  max_intensity = max(trws_data_norm(:)) - min_intensity;
  trws_data_norm = (trws_data_norm - min_intensity + .01)/max_intensity*10;
  
  %% HISTOGRAM
  if PLOT_HISTOGRAM
    figure;
    histogram(trws_data_norm(:));
    title(sprintf('loops: %d', MAX_LOOPS));
  end

  %% MAX SURFACE
  if PLOT_MAX_SURF
    [max_surf_vals, max_surf] = max(trws_data_norm, [], 2);
    max_surf = squeeze(max_surf);
    max_surf_vals = squeeze(max_surf_vals);
    
    % Remove max surface outside extreme boundaries. Those points do not
    %    affect the result as they are set to nan before TRWS2 is called.
    max_surf = [nan(min_bound - 1, Nx); max_surf(min_bound:max_bound, :); nan(size(trws_data_norm, 1) - max_bound, Nx)];
%     for rline = 1:Nx
%       for rbin = 1:Nt
%         doa_bin = max_surf(rbin, rline);
%         if rbin < min_bounds(doa_bin, rline) + 1 || rbin > max_bounds(doa_bin, rline) + 1
%           max_surf_vals(rbin, rline) = NaN;
%           max_surf(rbin, rline) = NaN;
%         end
%       end
%     end
  end
end

%% CREATE FIGURE
if ~isnan(FIGURE_NUM)
  fig = figure(FIGURE_NUM);
else
  fig = figure;
end

clf;
hold on;

Z = 1:Nt;
Y = 1:Nsv;
X = 1:Nx;

xlims = [1 Nx];
zlims = [0 Nsv];
ylims = [1 Nt];

clear surf;

%% PLOT MAX SURF
if LOAD_DATA && PLOT_MAX_SURF
  h_max_surf = surf(max_surf, 'FaceAlpha', .6, 'FaceColor', 'interp', 'LineStyle', 'none');
  set(h_max_surf, 'CData', max_surf_vals);
  caxis([min(max_surf_vals(:)) max(max_surf_vals(:))]);
  c_bar = colorbar;
  ylabel(c_bar, 'Normalized Intensity');
end

%% PLOT SURFACE
h_correct_surf = surf(correct_surface, 'FaceAlpha', .6, 'FaceColor', 'black', 'LineStyle', 'none');
%   set(h_correct_surf, 'CData', get(h_correct_surf, 'ZData'));

%% PLOT CT BOUNDS
if PLOT_CT_BOUNDS
  blue   = [86 , 135, 214]./255;
  orange = [232, 145, 90 ]./255;
  surf(X, Z, repmat(CT_bounds(1, :), Nt, 1), 'FaceColor', blue, 'LineStyle', 'none', 'FaceAlpha', 0.1);
  surf(X, Z, repmat(CT_bounds(2, :), Nt, 1), 'FaceColor', orange, 'LineStyle', 'none', 'FaceAlpha', 0.1);
end

%% Set up camera view and axes
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

view([40 90 15]);
camorbit(90, 0, 'data', [0 0 1]);
camorbit(90, 0, 'data', [1 0 0]);
cameratoolbar('SetCoordSys', 'y');
cameratoolbar('SetMode', 'orbit');

camva(10);
% colormap(bone);

%% COLOR SURFACE
if ~strcmp(COLOR_SURF_DATA, 'none')
  color_data = eval(COLOR_SURF_DATA);
  intensities = nan(Nt, Nx);
  for d = 1:Nt
    for w = 1:Nx
      h = correct_surface(d, w);
      if isnan(h)
        continue
      end
      intensities(d, w) = color_data(d, h, w);
    end
  end
  set(h_correct_surf, 'CData', intensities);
  set(h_correct_surf, 'FaceColor', 'interp');
  caxis([min(color_data(:)) max(color_data(:))]);
  c_bar = colorbar;
  ylabel(c_bar, 'Input Intensity');
end

%% PLOT FT BOUNDS
if PLOT_FT_BOUNDS
  % Remove data entirely outside ft_bounds
  ylim(zlims);
  zlim(ylims);
  minh = surf(min_bounds, 'FaceColor', 'none', 'LineStyle', ':', 'FaceAlpha', 0.1);
  rotate(minh, [1 0 0], 90, [0 Nt 0]);
  set(minh, 'YData', Nt - get(minh, 'YData'));
  % set(minh, 'CData', get(minh, 'YData'));
  set(minh, 'ZData', get(minh, 'ZData') - min(get(minh, 'ZData')) + 1);
  maxh = surf(max_bounds, 'FaceColor', 'none', 'LineStyle', ':', 'FaceAlpha', 0.1);
  rotate(maxh, [1 0 0], 90, [0 Nt 0]);
  set(maxh, 'YData', Nt - get(maxh, 'YData'));
  set(maxh, 'ZData', get(maxh, 'ZData') - min(get(maxh, 'ZData')) + 1);
  % set(maxh, 'CData', get(maxh, 'YData'));
  ylim(ylims);
  zlim(zlims);
end

%% Plot max points if PLOT_MAX_POINTS is true
if PLOT_MAX_POINTS && LOAD_DATA
  for w_idx = 1:Nx
    for d_idx = 1:Nt
      points = trws_data_norm(d_idx, :, w_idx);
      points = points(~isnan(points));
      num_points = MAX_POINTS_NUM;
      if numel(points) < MAX_POINTS_NUM && numel(points) > 0
        num_points = numel(points); % Only a few points that are not NaN
      elseif numel(points) == 0
        continue; % All points in this row are NaN
      end
        
      [sorted_points, indices] = sort(points);
      greatest_points = indices((numel(indices)-num_points + 1):end);
      intensities = points(greatest_points);
      color = 'm';
      for point = 1:numel(intensities)
        if point == numel(intensities)
          color = 'r';
        end
        plot3(w_idx, d_idx, greatest_points(point), sprintf('%s.', color), 'MarkerSize', intensities(point));
      end
    end
  end
end

%% Set threshold for plot_points if set to auto
if PLOT_POINTS && strcmp(PLOT_THRESHOLD, 'auto') && LOAD_DATA
  sorted_data = trws_data_norm(~isnan(trws_data_norm));
  [values, indices] = sort(sorted_data(:));
  PLOT_THRESHOLD = values(numel(values)-AUTO_THRESHOLD_MAX_POINTS+1);  

%   [counts, edges] = histcounts(trws_data_norm(:), 'BinMethod', 'sqrt');
%   [~, biggest] = max(counts);
%   if biggest < numel(counts)
%     biggest = biggest + 1;
%   end
%   PLOT_THRESHOLD = edges(biggest)
end

%% Plot indices and points

Nsvs_center = floor(Nsv/2);
Nsvs_array = [Nsvs_center:-1:1 (Nsvs_center+1):Nsv];

for w_idx = 1:Nx
  if mod(INDEX_LOOP_NUM, 2) == 1
    w = Nx-w_idx;
  else
    w = w_idx-1;
  end
  for h_idx = 1:Nsv
    if mod(floor(INDEX_LOOP_NUM/2), 2) == 1
      h = Nsv-h_idx;
    else
      h = h_idx-1;
    end
    
    for d_idx = 1:Nt
      d = d_idx-1;
      
      % Plot points
      if PLOT_POINTS && ~PLOT_MAX_POINTS && LOAD_DATA
        intensity = trws_data_norm(d_idx, h_idx, w_idx);
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

%% Add variable values to plot
current_axes = gca;
current_axes.Position(3) = .5;
current_axes.Position(1) = .3;


str = '';
for display_cell = DISPLAY_VARS
  display_var = display_cell{1};
  value = eval(display_var);
  
  if islogical(value)
    if value
      str = sprintf('%s%s\n', str, display_var);
    end
  else
    if isnumeric(value)
      template = '%s = %.2f';
    elseif isstring(value) || ischar(value)
      template = '%s = %s';
    else
      warning('No template for %s (class %s)', display_var, class(value));
      continue;
    end
    str = sprintf('%s%s\n', str, sprintf(template, display_var, value));
  end
end
annotation('textbox', [0, 0, 0.3, 1], 'String', str, 'Interpreter', 'none', 'FontSize', 6, 'LineStyle', 'none');

%% SAVE OUTPUTS

if exist([SAVE_PATH SAVE_NAME], 'dir') ~= 7 % 7 is a folder
  mkdir([SAVE_PATH SAVE_NAME]);
end

output_path = [SAVE_PATH SAVE_NAME '/'];

if SAVE_IMAGE
  Tomo.img = trws_data_norm;
  Time = mdata.Time(St:Tt:Et);
  save([output_path 'image.mat'], 'Tomo', 'Time');
end
if SAVE_SURF
  save([output_path 'surf.mat'], 'correct_surface');
end
if SAVE_FIG
  savefig(fig, [output_path 'fig.fig']);
end


