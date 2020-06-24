Nt  = 5;
Nsv = 6;
Nx  = 7;

% 'ALL' displays index numbers for every cell
% 'FT' displays traversal order for a fast-time-based search
% 'CT' displays traversal order for a cross-track-based search
INDEX_LABEL = 'ALL';
INDEX_EVEN_LOOP = true;


trws_data = rand(Nt, Nsv, Nx);

trws_data(3, :, :) = 10;
trws_data(2, 3:4, 2:6) = 10;
trws_data(3, 3:4, 2:6) = 0;
% trws_data(1, :, :) = 9;
trws_data(1, 4, 4) = 30;


at_slope  = zeros(1, Nx);
at_weight = 1;
ct_slope  = zeros(Nsv, Nx);
ct_weight = ones(1, Nsv);

bounds = zeros(2, Nx);
bounds(2, :) = Nt - 1;
bounds(1, 3:5) = 2;

max_loops = 2;

correct_surface = tomo.trws2(single(trws_data),single(at_slope),single(at_weight),single(ct_slope),single(ct_weight), uint32(max_loops), uint32(bounds));

[X, Y] = meshgrid(1:Nx, 1:Nsv);

figure(1);
clf;
hold on;

surf(X, Y, correct_surface);
shading interp;
colormap(parula);

xlim([1 Nx]);
% xticks(1:Nx);
xlabel('X : Along-Track (Nx, DIM 2)');
set(gca, 'XColor', 'b');

ylim([1 Nsv]);
% yticks(1:Nsv);
ylabel('Y : Cross-Track (Nsv. DIM 1)');
set(gca, 'YColor', 'g');

zlim([0 Nt]);
% zticks(0:Nt);
zlabel('Z : Fast-Time (Nt, DIM 0)');
set(gca, 'ZColor', 'r');
set(gca, 'zdir', 'reverse');


surf(X, Y, repmat(bounds(1, :), Nx-1, 1), 'FaceColor', [86, 135, 214]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);
surf(X, Y, repmat(bounds(2, :), Nx-1, 1), 'FaceColor', [232, 145, 90]./255, 'LineStyle', 'none', 'FaceAlpha', 0.2);

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
                plot3(w_idx, h_idx, d_idx, sprintf('%s.', color), 'MarkerSize', intensity);
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
            
            if ~isnan(msg_idx)
                text(w_idx, h_idx, d_idx, sprintf("%d", msg_idx));
            end
        end
    end
end


