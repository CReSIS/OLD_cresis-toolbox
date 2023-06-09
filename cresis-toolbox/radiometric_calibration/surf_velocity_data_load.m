hm;

fn = fullfile('X:\ct_data\ct_tmp\radcal_itslive\surf_velocity\2017\GRE_G0240_2017.nc');

% raw_v = importdata();
% 
% vars = fieldnames(raw_v);
% for i = 1:length(vars)
%     assignin('base', vars{i}, raw_v.(vars{i}));
% end

%% laod
info = ncinfo(fn);
var = struct();
for idx = 1:length(info.Variables)
    var.(info.Variables(idx).Name) = ncread(fn, info.Variables(idx).Name);
end

%% 
projection = projcrs(3413); % EPSG
[X, Y] = meshgrid(var.x,var.y);
[lat, lon] = projinv(projection, X, Y);
save(gcf, 'temp');
save(gcf, 'temp', '-dpng');