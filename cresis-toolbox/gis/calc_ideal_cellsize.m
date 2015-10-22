%% Script Details/Information

% =============================================================
% Title: calc_ideal_cellsize.m
% About: Calculates the ideal interpolation cell size based on
%        the "percentage rule" eg. Atleast 50% of the resulting
%        cells must contain atleast one point from the point data.
% Author: Kyle Purdon
% Version: V2.0 11/3/2011
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "calc ideal cellsize tool"
% Toolboxes Used: Standard, Mapping, Signal, Stats
% Known Bugs: None
% Planned Updates: Efficiency Modifications
% Additional Information:
% =============================================================

%% User Input

% Path to the "StudyArea" shapefile (Path must have .shp,.dbf,.shx)

boundpath = 'C:\Users\SomeFolder\SomeShapefileBoundary.shp';

% Path to the "Flightline" shapefile (Path must have .shp,.dbf,.shx)

flightpath = 'C:\Users\SomeFolder\SomeShapefileFlightlines.shp';

% Is the flightline shapefile you provided decimated?
% If "No" the script will downsample it for efficiency.
% 1 = yes and 2 = no

fltshp_dec = 2;

%Factor to decimate by. Keep every Nth point.
%Default = 10
%General Range ~10-30 (Higher=Faster Lower=Precise)

fltshp_dec_factor = 10;

% Maximum Cell Size and Increment to Test
% The closer this is to the resulting cell size the faster the script will
% run. (Make a realistic guess. ~2000m)
% For faster processing do a first run with a large increment. (50-100m)
% then decrease the increment and update the maxCellSize for a precise run.

maxCellSize = 2000;
CellSizeIncr = 25;  %Default ~ 25-50m (Higher=Faster Lower=Precise)

% Percentage Value to Test
% N percent of cells must contain at least 1 point from the flightlines
% CReSIS Default = 50% or 0.50

prctToTest = 0.50;

% Show plots after completion?
% If "Yes" the script will plot the bounding "study Area" the
% Flightlines and the center points of the resulting cells.
% 1 = yes and 2 = no    Default = 1;

showPlots = 1;


%% Automated Section

% Main Program
clear xf yf x y;
fprintf('================================\n');
fprintf('      Cell Size Calculator      \n');
fprintf('================================\n');

%Read the input shapefiles
fprintf('Opening the shapefiles... ');
boundshp = shaperead(boundpath);
flightshp = shaperead(flightpath,'Attributes',{});
fprintf('%s\n',datestr(now,'HH:MM:SS'));

%Create flightline points and downsample to every 10th point for
%efficiency.
xf = zeros(1,length(flightshp));
yf = zeros(1,length(flightshp));
for fdat_idx = 1:length(flightshp)
  xf(1,fdat_idx) = flightshp(fdat_idx).X;
  yf(1,fdat_idx) = flightshp(fdat_idx).Y;
end

% Decimate the flightlines if they are not already.
if fltshp_dec == 2
  xf = downsample(xf,fltshp_dec_factor);
  yf = downsample(yf,fltshp_dec_factor);
end

%Put the bounding box from the "Study Area" into variables
minx = boundshp.BoundingBox(1,1);
miny = boundshp.BoundingBox(1,2);
maxx = boundshp.BoundingBox(2,1);
maxy = boundshp.BoundingBox(2,2);

fprintf('Cell Size Calculations Starting...\n');
fprintf('--------------------------------\n');
fprintf('Cell Size ------ Percentage\n');
goodpct = inf;
while goodpct > prctToTest
  for cellsize = maxCellSize:-CellSizeIncr:0
    clear x y
    xdiff = maxx-minx;
    xlen = (xdiff/cellsize);
    ydiff = maxy-miny;
    ylen = (ydiff/cellsize);
    
    %Derive X and Y vecotrs for grid creation
    for idx_x = 1:xlen-1
      x(1,idx_x+1) = minx+(cellsize*idx_x)+(0.5*cellsize);
    end
    x(1) = x(2)-cellsize;
    for idx_y = 1:ylen-1                                                    %
      y(idx_y+1,1) = miny+(cellsize*idx_y)+(0.5*cellsize);
    end
    y(1) = y(2)-cellsize;
    
    %Create Grids
    yg = repmat(y,1,length(x));
    xg = repmat(x,length(y),1);
    [xg1 xg2] = size(xg);
    
    %Make any points outside the polygon boundary NaN
    IN = inpolygon(xg,yg,boundshp.X,boundshp.Y);
    for ridx = 1:xg1
      for cidx = 1:xg2
        if IN(ridx,cidx) == 0
          xg(ridx,cidx) = NaN;
          yg(ridx,cidx) = NaN;
        end
      end
    end
    
    fprintf('  %dm ---------- ',cellsize);
    %Preallocate matrix for distances
    ptdist = zeros(xg1*xg2,length(xf));
    minDist = zeros(xg1*xg2,1);
    cellsWpts = zeros(xg1*xg2,1);
    
    %Get distances
    for idxg = 1:xg1*xg2
      for idxp = 1:length(xf)
        ptdist(idxg,idxp) = sqrt((xf(idxp)-xg(idxg))^2+(yf(idxp)-yg(idxg))^2);
      end
      minDist(idxg) = min(ptdist(idxg,:));
      if ~isnan(minDist(idxg)) && minDist(idxg) < cellsize
        cellsWpts(idxg) = 1;
      elseif isnan(minDist(idxg))
        cellsWpts(idxg) = NaN;
      end
    end
    goodpct = nanmean(cellsWpts);
    fprintf('%2.2f%%\n',goodpct*100);
    if goodpct <= prctToTest
      break
    end
  end
end

idealCellSize = cellsize+CellSizeIncr;

fprintf('--------------------------------\n');
fprintf('Calculation Complete ... %s \n',datestr(now,'HH:MM:SS'));
fprintf('--------------------------------\n');
fprintf('Ideal Cell Size: %dm\n',idealCellSize);
fprintf('--------------------------------\n');

%Verification plots
if showPlots == 1
  clf;
  fprintf('Creating plots ... ');
  plot(xg,yg,'r+');
  hold on;
  plot(boundshp.X,boundshp.Y,'k-');
  plot(xf,yf,'b.');
  axis([minx maxx miny maxy])
  fprintf('%s\n',datestr(now,'HH:MM:SS'));
end

