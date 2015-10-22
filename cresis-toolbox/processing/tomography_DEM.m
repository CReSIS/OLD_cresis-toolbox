% Run after tomography_post.m

date

if 0
filenames = {};

base_path = '/N/dc/projects/cresis/output/mcords/2011_Antarctica_DC8/CSARP_tomo/';
for frm = [1:26,28:30]
% for frm = [28:30]
  filenames{end+1} = sprintf('20111014_07/Data_img_01_20111014_07_%03i',frm);
end

% base_path = '/N/dcwan/projects/cresis/output/mcords2/2011_Greenland_P3/CSARP_music/';
% First section (north/east to south/west)
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_010';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_029';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_011';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_026';
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_007';
% Second section
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_011';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_030';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_012';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_025';
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_006';
% % Third
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_012';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_031';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_013';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_024';
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_005';
% % Fourth
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_013';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_032';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_014';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_023';
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_004';
% % Fifth
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_014';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_033';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_015';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_022';
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_003';
% % Sixth
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_015';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_034';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_016';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_021';
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_002';
% % Seventh
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_016';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_035';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_017';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_020';
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_001';
% % Eighth
% filenames{end+1} = '20110506_02/Data_img_02_20110506_02_017';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_036';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_018';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_019';
% filenames{end+1} = '20110506_01/Data_img_02_20110506_01_037';

zlimits = [-1500 500];
% zlimits = [-500 1200];
% zlimits = [-300 175];

out_fn = '/N/u/jpaden/Quarry/PIG_20111014_07.mat';
% out_fn = '/cresis/scratch2/mdce/mcords2/2011_Greenland_P3/CSARP_tomography/russell.mat';
% out_fn = '/N/u/jpaden/Quarry/NEEM_NGRIP_3d_1to2.mat';

minEast = inf;
maxEast = -inf;
minNorth= inf;
maxNorth = -inf;
meas = {};
for measInd = 1:length(filenames)
  meas{measInd} = load(fullfile(base_path,filenames{measInd}));
  if min(meas{measInd}.east) < minEast
    minEast = min(meas{measInd}.east);
  end
  if max(meas{measInd}.east) > maxEast
    maxEast = max(meas{measInd}.east);
  end
  if min(meas{measInd}.north) < minNorth
    minNorth = min(meas{measInd}.north);
  end
  if max(meas{measInd}.north) > maxNorth
    maxNorth = max(meas{measInd}.north);
  end
end
minEast = minEast - 1000;
maxEast = maxEast + 1000;
minNorth = minNorth - 1000;
maxNorth = maxNorth + 1000;

tic;
debugLevel = 1;
physical_constants;

grid_spacing = 25;
eastAxis = minEast:grid_spacing:maxEast;
northAxis = (minNorth:grid_spacing:maxNorth).';
[eastMesh,northMesh] = meshgrid(eastAxis,northAxis);

for measInd = 1:length(meas)
  fprintf('Creating DEM for measurement %d (%.1f sec)\n', measInd, toc);
  xAngle = -(meas{measInd}.angleMean + meas{measInd}.yAngle);
  pc = [meas{measInd}.east; meas{measInd}.north; meas{measInd}.elev - meas{measInd}.surface*c/2];
  yUnit = [sin(xAngle); cos(xAngle)];
  pnts = zeros(3,size(meas{measInd}.body,1),size(meas{measInd}.body,2));
  cpnts = zeros(1,size(meas{measInd}.body,1),size(meas{measInd}.body,2));
  for lineInd = 1:size(meas{measInd}.body,2)
    pnts(1,:,lineInd) = pc(1,lineInd) + yUnit(1,lineInd)*meas{measInd}.yAxis;
    pnts(2,:,lineInd) = pc(2,lineInd) + yUnit(2,lineInd)*meas{measInd}.yAxis;
    pnts(3,:,lineInd) = pc(3,lineInd) + meas{measInd}.body(:,lineInd);
    cpnts(1,:,lineInd) = meas{measInd}.yAxis;
  end

  % Create a constrained delaunay triangulization that forces edges
  % along the boundary (concave_hull) of our swath
  concave_hull = [ [(1:size(pnts,2)-1).' (2:size(pnts,2)).'];
    size(pnts,2)*[(1:size(pnts,3)-1).' (2:size(pnts,3)).'];
    size(pnts,2)*(size(pnts,3)-1) + [(size(pnts,2):-1:2).' (size(pnts,2)-1:-1:1).'];
    size(pnts,2)*[(size(pnts,3):-1:2).' (size(pnts,3)-1:-1:1).']-size(pnts,2)+1 ];
  dt = DelaunayTri(pnts(1,:).',pnts(2,:).',concave_hull);

  % Remove too large triangles: DOES NOT WORK
%   bad_mask = zeros(size(dt.Triangulation,1),1);
%   for tri_idx = 1:size(dt.Triangulation,1)
%     if mod(tri_idx,1000)
%       tri_idx
%     end
%     if sqrt(sum(abs(dt.X(dt.Triangulation(tri_idx,1),:)-dt.X(dt.Triangulation(tri_idx,2),:)).^2)) > 2*grid_spacing ...
%         || sqrt(sum(abs(dt.X(dt.Triangulation(tri_idx,1),:)-dt.X(dt.Triangulation(tri_idx,3),:)).^2)) > 2*grid_spacing ...
%         || sqrt(sum(abs(dt.X(dt.Triangulation(tri_idx,3),:)-dt.X(dt.Triangulation(tri_idx,2),:)).^2)) > 2*grid_spacing
%       bad_mask(tri_idx) = 1
%     end
%   end
%   dt.Triangulation = dt.Triangulation(~bad_mask,:); % Not allowed
  
  F = TriScatteredInterp(dt,pnts(3,:).');
  %F = TriScatteredInterp(dt,pnts(3,:).','natural');
  meas{measInd}.DEM = F(eastMesh,northMesh);
  
  fprintf('  Creating boundary and removing outside points (%.1f sec)\n', toc);
  if 0
    % Slow method using inpolygon
    boundary = [pnts(1:2,:,1) squeeze(pnts(1:2,end,:)) squeeze(pnts(1:2,end:-1:1,end)) squeeze(pnts(1:2,1,end:-1:1))];
    good_mask = ~isnan(meas{measInd}.DEM);
    good_mask_idx = find(good_mask);
    in = inpolygon(eastMesh(good_mask),northMesh(good_mask),boundary(1,:),boundary(2,:));
    good_mask(good_mask_idx(~in)) = 0;
    meas{measInd}.DEM(~good_mask) = NaN;
  else
    % Faster method using inOutStatus/pointLocation with edge constraints
    % Finds the triangles outside the concave hull (the bad ones)
    bad_tri_mask = ~inOutStatus(dt);
    % Selects just the points that were in the convex hull
    good_mask = ~isnan(meas{measInd}.DEM);
    good_mask_idx = find(good_mask);
    % For each point in the convex hull, find the triangle enclosing it
    tri_list = pointLocation(dt,eastMesh(good_mask),northMesh(good_mask));
    % Use the bad triangle list to find the bad points
    bad_mask = bad_tri_mask(tri_list);
    % Set all the points in bad triangles to NaN
    good_mask(good_mask_idx(bad_mask)) = 0;
    meas{measInd}.DEM(~good_mask) = NaN;
  end
  fprintf('  Done creating boundary and removing outside points (%.1f sec)\n', toc);
  
  F = TriScatteredInterp(dt,cpnts(1,:).');
  meas{measInd}.COST = F(eastMesh,northMesh);
  meas{measInd}.COST(~good_mask) = NaN;
  if debugLevel == 1
    figure(measInd); clf;
    set(measInd,'Position',[200 200 560 350]);
    imagesc((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,meas{measInd}.DEM,zlimits);
    set(gca,'YDir','normal');
    xlabel('Easting (km)');
    ylabel('Northing (km)');
    grid on;
    hold on;
    plot((meas{measInd}.east-eastAxis(1))/1e3,(meas{measInd}.north-northAxis(1))/1e3,'k-');
    plot((meas{measInd}.east(1)-eastAxis(1))/1e3,(meas{measInd}.north(1)-northAxis(1))/1e3,'ko');
%     h = plot((neem.east-eastAxis(1))/1e3,(neem.north-northAxis(1))/1e3,'r*');
%     set(h, 'LineWidth', 3);
%     set(h, 'MarkerSize', 7);
    hold off;
  end
  save(fullfile(base_path,filenames{measInd}),'-append','COST','DEM');
end
fprintf('Done (%.1f sec)\n', toc);

if debugLevel == 2
  figure(20); clf;
  imagesc((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,abs(meas{8}.DEM-meas{3}.DEM),[10 50]);
  set(gca,'YDir','normal');
  xlabel('Easting (km)');
  ylabel('Northing (km)');
  grid on;
  figure(21); clf;
  imagesc((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,abs(meas{5}.DEM-meas{3}.DEM),[10 50]);
  xlabel('Easting (km)');
  ylabel('Northing (km)');
  grid on;
end

end

DEM = NaN*zeros(size(eastMesh),'single');
NUM_SAMPLES = zeros(size(eastMesh),'single');
COST = zeros(size(eastMesh),'single');
val = zeros(11,1);
cost = zeros(11,1);
%for ind = sub2ind(size(eastMesh),105,96):numel(DEM)
for ind = 1:numel(DEM)
  if ~mod(ind-1,10000)
    fprintf('Index %d of %d (%.1f sec)\n', ind, numel(DEM), toc);
    %figure(10);clf;
    %imagesc(DEM);
    %keyboard;
  end
  for measInd = 1:length(meas)
    if abs(meas{measInd}.COST(ind)) < 1050
      NUM_SAMPLES(ind) = NUM_SAMPLES(ind) + 1;
      val(NUM_SAMPLES(ind)) = meas{measInd}.DEM(ind);
      cost(NUM_SAMPLES(ind)) = 1101-abs(meas{measInd}.COST(ind));
    end
  end
  if NUM_SAMPLES(ind) == 1
    DEM(ind) = val(1);
    COST(ind) = cost(1);
  elseif NUM_SAMPLES(ind) >= 2
    sumCost = sum(cost(1:NUM_SAMPLES(ind)));
    COST(ind) = max(cost(1:NUM_SAMPLES(ind)));
    DEM(ind) = sum(val(1:NUM_SAMPLES(ind)).*cost(1:NUM_SAMPLES(ind)))/sumCost;
  end
end
DEM_orig = DEM;

% Remove and/or filter areas where DEM cost function suggests poor
% performance
DEM(NUM_SAMPLES<2 & COST<200) = NaN;
DEM_filt1 = filter2(ones(3)/9, DEM);
DEM_filt2 = filter2(ones(7,3)/(7*3), DEM);
DEM(COST<450) = DEM_filt1(COST<450);
DEM(COST<300) = DEM_filt2(COST<300);

% Find GRIP and GISP2 bed elevations
% [minEDist minEInd] = min(abs(eastAxis-neem.east));
% [minNDist minNInd] = min(abs(northAxis-neem.north));
% neem.DEM = DEM(minNInd,minEInd);

out_fn_dir = fileparts(out_fn);
if ~exist(out_fn_dir,'dir')
  mkdir(out_fn_dir);
end
%save(out_fn,'-v6','meas','DEM_orig','DEM','eastAxis','northAxis');
save(out_fn,'-v6','DEM','eastAxis','northAxis');
zlims = [-110 -30];

% Fill holes in DEM
mask = ~isnan(DEM);
mask = shrink(grow(mask));
DEM_interp = DEM;
DEM_interp(isnan(DEM_interp)) = 1;
DEM_interp = medfilt2(DEM_interp,[3 3]);
DEM(isnan(DEM) & mask) = DEM_interp(isnan(DEM) & mask);

figure(20); clf;
hA2 = axes;
hC = surf((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,double(medfilt2(DEM,[3 11])));
set(hC,'EdgeAlpha',0.2); grid off;
set(hA2,'View',[190 75]);
hold on;
% h = plot3((neem.east-eastAxis(1))/1e3,(neem.north-northAxis(1))/1e3,neem.DEM,'ko');
% set(h, 'LineWidth', 3);
% set(h, 'MarkerSize', 7);
hold off;
hx = xlabel('Easting (km)');
set(hx,'Position',[11 16 -1500]);
set(hx,'Rotation',-3);
hy = ylabel('Northing (km)');
set(hy,'Position',[-2 10 -1500]);
set(hy,'Rotation',67);
zlabel('WGS-84 bed elevation (m)');
grid on;
hC = colorbar;
set(get(hC,'YLabel'),'String','WGS-84 bed elevation (m)');
set(hA2,'Position',[0.12 0.11 0.72 0.815])
set(hC,'Position',[0.9 0.11 0.022 0.815])

figure(21); clf;
set(21,'Position',[200 200 560 350]);
imagesc((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,double(DEM));
hold on;
% h = plot((neem.east-eastAxis(1))/1e3,(neem.north-northAxis(1))/1e3,'ko');
% set(h, 'LineWidth', 3);
% set(h, 'MarkerSize', 7);
hold off;
xlabel('Easting (km)');
ylabel('Northing (km)');
grid on;
set(gca,'YDir','normal');
h = colorbar;
set(get(h,'YLabel'),'String','WGS-84 bed elevation (m)');

return;

saveas(20, '/users/paden/tmp/DEM_surface.fig');
saveas(21, '/users/paden/tmp/DEM_image.fig');































% DEM surface
% phong value for the FaceLighting and EdgeLighting 
% =========================================================================
% =========================================================================
figure(1); clf;
set(1,'Position',[50 50 700 400]);
set(1,'Color',[1 1 1]);

hA2 = axes;
%hC = contourf((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,double(DEM),12);
hC = surf((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,double(DEM)*0+25,double(DEM));
set(hC(1),'EdgeAlpha',0); grid off;
%axis([2 38 0 10 zlims]);
set(hA2,'Box','off');
set(hA2,'View',[7.5   74]);
set(hA2,'Position',[0.025 0.1 0.875 0.53]);
set(hA2,'YTick',[0 5 10]);
set(hA2,'ZColor',[1 1 1]);
colormap(jet(256))
hc = colorbar;
caxis(zlims);
set(hc,'Position',[0.95 0.1 0.015 0.8]);
set(get(hc,'YLabel'),'String','Bed height (m)');

hA = axes;
hS = surf((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3-1,double(DEM),double(1*DEM));
hA = gca; grid off;
%axis([2 38 0 10 zlims]);
set(hA,'Box','off');
% View ["" "elevation angle of observer"]
set(hA,'View',[7.5 74]);
set(hA,'Position',[0.025 0.45 0.875 0.66]);
set(hA,'XTick',[]);
set(hA,'YTick',[]);
set(hA,'ZTick',[]);
set(hA,'XColor',[1 1 1]);
set(hA,'YColor',[1 1 1]);
set(hA,'ZColor',[1 1 1]);
set(hS(1),'EdgeAlpha',0.2);

if 0 % surf with color
  colormap(jet(256))
  hc = colorbar;
  caxis(zlims);
  set(hc,'Position',[0.95 0.1 0.015 0.8]);
  set(get(hc,'YLabel'),'String','Bed height (m)');
end

axes(hA2);
hx = xlabel('Easting (km)');
%set(hx,'Position',[19 -1.8 0]);
hy = ylabel('Northing (km)');
%set(hy,'Position',[40 3.5 0]);
%set(hy,'Rotation',54);

axes(hA2);
hold on;
h = plot3((neem.east-eastAxis(1))/1e3,(neem.north-northAxis(1))/1e3,30,'ko');
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', 5);
hold off;

axes(hA);
hold on;
h = plot3((neem.east-eastAxis(1))/1e3,(neem.north-northAxis(1))/1e3-1,neem.DEM+12,'ro');
set(h, 'LineWidth', 2);
set(h, 'MarkerSize', 5);
hold off;

% =========================================================================
% =========================================================================
figure(3); clf;
set(3,'Position',[50 50 500 500]);
set(3,'Color',[1 1 1]);

hS = surf((eastAxis-eastAxis(1))/1e3-2,(northAxis-northAxis(1))/1e3-1,double(DEM),double(DEM));
hA = gca; grid off;
%axis([0 36 0 10 zlims]);
set(hA,'Box','off');
% View ["" "elevation angle of observer"]
set(hA,'YTick',[0 5 10]);
set(hA,'ZColor',[1 1 1]);
set(hA,'View',[7.5 74]);
set(hA,'Position',[0.025 0.135 0.875 1.1]);
set(hA,'ZColor',[1 1 1]);
set(hS(1),'EdgeAlpha',0.2);

colormap(jet(256))
hc = colorbar;
set(hc,'Position',[0.95 0.1 0.015 0.8]);
set(get(hc,'YLabel'),'String','Bed height (m)');
caxis(zlims);

hx = xlabel('Easting (km)');
set(hx,'Position',[19 -1.2 0]);
set(hx,'Rotation',-2);
hy = ylabel('Northing (km)');
set(hy,'Position',[38 3.5 0]);
set(hy,'Rotation',68);

hold on;
h1 = plot3((neem.east-eastAxis(1))/1e3-2,(neem.north-northAxis(1))/1e3-1,neem.DEM*0+diff(zlims),'kx');
set(h1, 'LineWidth', 3);
set(h1, 'MarkerSize', 7);
hold off;


% =========================================================================
% =========================================================================
DEM_contour = DEM;
DEM_contour(isnan(DEM_contour)) = 10;
figure(2); clf
contourf((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,DEM_contour,linspace(25,125,10));
axis([22 32 1.5 9]);
hx = xlabel('Easting (km)');
hy = ylabel('Northing (km)');
set(2,'Position',[50 50 650 475]);
colormap(jet(256));
hc = colorbar;
set(get(hc,'YLabel'),'String','Bed height (m)');

hold on;
h = plot((neem.east-eastAxis(1))/1e3,(neem.north-northAxis(1))/1e3,'ko');
set(h, 'LineWidth', 2);
hold off;

return;

figure(9); clf;
imagesc((eastAxis-eastAxis(1))/1e3,(northAxis-northAxis(1))/1e3,meas{9}.DEM);
hold on;
plot((meas{9}.east-eastAxis(1))/1e3,(meas{9}.north-northAxis(1))/1e3,'k-');
plot((meas{9}.east(1)-eastAxis(1))/1e3,(meas{9}.north(1)-northAxis(1))/1e3,'ko');
hold off;
grid on;










