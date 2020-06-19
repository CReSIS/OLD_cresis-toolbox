Nt  = 5;
Nsv = 6;
Nx  = 7;


trws_data = zeros(Nt, Nsv, Nx);

trws_data(4, :, :) = 1;
trws_data(2, 3:4, 2:6) = 1;
trws_data(4, 3:4, 2:6) = 0;

at_slope  = zeros(1, Nx);
at_weight = 0;
ct_slope  = zeros(Nsv, Nx);
ct_weight = zeros(1, Nsv);

bounds = zeros(2, Nx);
bounds(2, :) = Nt - 1;

max_loops = 2;

correct_surface = tomo.trws2(single(trws_data),single(at_slope),single(at_weight),single(ct_slope),single(ct_weight), uint32(max_loops), uint32(bounds));

[X, Y] = meshgrid(1:Nx, 1:Nsv);

figure(1);
clf;
hold on;

surf(X, Y, correct_surface);
shading interp
colormap(parula);

xlim([1 Nx]);
xticks(1:Nx);
xlabel('X : Along-Track (Nx)');

ylim([1 Nsv]);
yticks(1:Nsv);
ylabel('Y : Cross-Track (Nsv)');

zlim([0 Nt]);
zticks(0:Nt);
zlabel('Z : Fast-Time (Nt)');

for x = 1:Nx
    for y = 1:Nsv
        for z = 1:Nt
            intensity = trws_data(z, y, x);
            if intensity > 0
                plot3(x, y, z, 'r.', 'MarkerSize', intensity*30);
            end
        end
    end
end


