Nt  = 5;
Nsv = 6;
Nx  = 7;

% 'ALL' displays index numbers for every cell
% 'FT' displays traversal order for a fast-time-based search
% 'CT' displays traversal order for a cross-track-based search
INDEX_LABEL = 'CT';
% 'NONE' does not plot a surface
% 'FT' plots a surface normal to the fast-time axis
% 'CT' plots a surface normal to the cross-track axis
SURFACE = 'CT';
INDEX_EVEN_LOOP = true;

% Display on this figure
FIGURE_NUM = 2;


trws_data = zeros(Nt, Nsv, Nx)/5;



Z = 1:Nt;
Y = 1:Nsv;
X = 1:Nx;

if strcmp(SURFACE, 'FT')
    find_surf = @tomo.trws2;
    
    trws_data(3, :, :) = 10;
    trws_data(2, 3:4, 2:6) = 10;
    trws_data(3, 3:4, 2:6) = 0;
    % trws_data(1, :, :) = 9;
    trws_data(1, 4, 4) = 30;
    
    bounds = zeros(2, Nx);
    bounds(2, :) = Nt - 1;
    bounds(1, 3:5) = 2;
elseif strcmp(SURFACE, 'CT')
    find_surf = @tomo.trws2_CT;
    
    trws_data(3, 3, 3) = 0;
    trws_data(3:4, 2, 2:6) = 10;
    trws_data(3:4, 3, 2:6) = 0;
    % trws_data(1, :, :) = 9;
%     trws_data(1, 1, 6) = 30;
%     trws_data(3, 4, 3) = 40;
    
    bounds = zeros(2, Nx);
    bounds(2, :) = Nsv - 1;
elseif ~strcmp(SURFACE, 'NONE')
    error('Invalid surface selection.');
end
    
% trws_data = reshape(1:Nt*Nsv*Nx, Nt, Nsv, Nx);


if ~strcmp(SURFACE, 'NONE')
    at_slope  = zeros(1, Nx);
    at_weight = 1;
    ct_slope  = zeros(Nsv, Nx);
    ct_weight = ones(1, Nsv);

    max_loops = 3;

    [correct_surface, debug] = find_surf(single(trws_data),single(at_slope),single(at_weight),single(ct_slope),single(ct_weight), uint32(max_loops), uint32(bounds));
end

figure(FIGURE_NUM);
clf;
hold on;

if ~strcmp(SURFACE, 'CT')
    if ~strcmp(SURFACE, 'NONE')
        surf(X, Y, correct_surface);
        surf(X, Y, repmat(bounds(1, :), Nx-1, 1), 'FaceColor', [86, 135, 214]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);
        surf(X, Y, repmat(bounds(2, :), Nx-1, 1), 'FaceColor', [232, 145, 90]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);
    end
    xlim([1 Nx]);
    xticks(1:Nx);
    xlabel('X : Along-Track (Nx, DIM 2)');
    set(gca, 'XColor', 'b');
    set(gca, 'xdir', 'reverse');

    ylim([1 Nsv]);
    yticks(1:Nsv);
    ylabel('Y : Cross-Track (Nsv. DIM 1)');
    set(gca, 'YColor', 'g');
    set(gca, 'ydir', 'reverse');

    zlim([0 Nt]);
    zticks(0:Nt);
    zlabel('Z : Fast-Time (Nt, DIM 0)');
    set(gca, 'ZColor', 'r');
    set(gca, 'zdir', 'reverse');

    view([-45 90 90]);
    cameratoolbar('SetCoordSys', 'none');
    cameratoolbar('SetMode', 'nomode');
    rotate3d on;
else
    if ~strcmp(SURFACE, 'NONE')
        surf(X, Z, correct_surface);
        surf(X, Z, repmat(bounds(1, :), Nsv-1, 1), 'FaceColor', [86, 135, 214]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);
        surf(X, Z, repmat(bounds(2, :), Nsv-1, 1), 'FaceColor', [232, 145, 90]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);
    end
    xlim([1 Nx]);
    xticks(1:Nx);
    xlabel('X : Along-Track (Nx, DIM 2)');
    set(gca, 'XColor', 'b');
    set(gca, 'xdir', 'reverse');

    zlim([0 Nsv]);
    zticks(0:Nsv);
    zlabel('Z : Cross-Track (Nsv. DIM 1)');
    set(gca, 'ZColor', 'g');
    
    ylim([1 Nt]);
    yticks(1:Nt);
    ylabel('Y : Fast-Time (Nt, DIM 0)');
    set(gca, 'YColor', 'r');
    set(gca, 'ydir', 'reverse');
    
    view([90 45 90]);
    camorbit(90, 0, 'data', [0 0 1]);
    camorbit(90, 0, 'data', [1 0 0]);
    cameratoolbar('SetCoordSys', 'y');
    cameratoolbar('SetMode', 'orbit');
end
camva(10);
shading interp;
colormap(bone);


Nsvs_center = floor(Nsv/2);
Nsvs_array = [Nsvs_center:-1:1 (Nsvs_center+1):Nsv];

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
            
            intensity = trws_data(d_idx, h_idx, w_idx);
            color = 'b';
            if intensity > 1
                color = 'r';
            end
            if intensity > 0
                if ~strcmp(INDEX_LABEL, 'CT')
                    plot3(w_idx, h_idx, d_idx, sprintf('%s.', color), 'MarkerSize', intensity);
                else
                    plot3(w_idx, d_idx, h_idx, sprintf('%s.', color), 'MarkerSize', intensity);
                end
            end
            
            msg_idx = NaN;
           
            if strcmp(INDEX_LABEL, 'ALL')
                msg_idx = d+(h_idx-1)*Nt+w*Nsv*Nt + 1;
            elseif strcmp(INDEX_LABEL, 'FT') && d_idx == 1
%                 idx = (h + w*Nsv)*Nt;
                msg_idx = h + w*Nsv + 1;
            elseif strcmp(INDEX_LABEL, 'CT') && h_idx == 1
%                 idx = d+w*Nsv*Nt;
                msg_idx = d+w*Nt + 1;
            end
            
%             if ~isnan(msg_idx)
                if ~strcmp(SURFACE, 'CT')
                    text(w_idx, h_idx, d_idx, sprintf("%f", intensity));
                else
                    value = debug(d_idx, h_idx, w_idx);
                    text_colormap = parula;
                    debug_min = min(debug(:));
                    debug_max = max(debug(:));
                    text_color = round((value-debug_min)/(debug_max-debug_min)*length(text_colormap));
                    text_color = ind2rgb(text_color+1,text_colormap);
                    text(w_idx, d_idx, h_idx, sprintf("%.0f", value), 'color', text_color);
                end
%             end
        end
    end
end
