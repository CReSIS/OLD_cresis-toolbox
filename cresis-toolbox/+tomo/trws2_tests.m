Nt  = 5;
Nsv = 6;
Nx  = 7;


trws_data = rand(Nt, Nsv, Nx);

trws_data(4, :, :) = 10;
trws_data(3, 3:4, 2:6) = 10;
trws_data(4, 3:4, 2:6) = 0;
% trws_data(1, :, :) = 9;
trws_data(1, 4, 4) = 30;



at_slope  = zeros(1, Nx);
at_weight = 1;
ct_slope  = zeros(Nsv, Nx);
ct_weight = ones(1, Nsv);

bounds = zeros(2, Nx);
bounds(2, :) = Nt - 1;

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
xlabel('X : Along-Track (Nx)');

ylim([1 Nsv]);
% yticks(1:Nsv);
ylabel('Y : Cross-Track (Nsv)');

zlim([0 Nt]);
% zticks(0:Nt);
zlabel('Z : Fast-Time (Nt)');

for x = 1:Nx
    for y = 1:Nsv
        for z = 1:Nt
            intensity = trws_data(z, y, x);
            color = 'b';
            if intensity > 1
                color = 'r';
            end
            if intensity > 0
                plot3(x, y, z, sprintf('%s.', color), 'MarkerSize', intensity);
            end
        end
    end
end


