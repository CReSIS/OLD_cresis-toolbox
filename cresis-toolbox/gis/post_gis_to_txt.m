%% Script Information/Details

% =============================================
% Writes TXT files for GIS L3 Projects
% Converts 3 ASCII to XYZ TXT and 1 DBF to TXT
% Author: Kyle Purdon
% Version: 05/16/2012
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "POSTGIS TO TXT"
% Additional Information:
% See also dbf_to_txt.m dbfread.m
% =============================================

%% User Input

% Paths to ESRI ASCII Grids(.TXT)
ascii_conv = true;
surface_fn = 'C:\SomeUsers\SomeProjectFolder\post\SomeSurface.txt';
thickness_fn = 'C:\SomeUsers\SomeProjectFolder\post\SomeThickness.txt';
bottom_fn = 'C:\SomeUsers\SomeProjectFolder\post\SomeBottom.txt';

% Path to Flightline DBF File
fl_conv = true;
flightline_dbf_fn = 'C:\SomeUsers\SomeProjectFolder\post\SomeFlightlines.dbf';

% Output Directory and Filenames
project_dir = 'C:\SomeUsers\SomeProjectFolder\post\';
output_xyz_fn = 'GlacierName_StartYear_EndYear_Composite_XYZGrid';
output_flight_fn = 'GlacierName_StartYear_EndYear_Composite_Flightlines';

%% Print Title Block
fprintf('======================================\n');
fprintf('          POST GIS TO TXT             \n');
fprintf('======================================\n');

%% Create Output Filenames
fprintf('Creating Output Filenames ... ');
if ascii_conv
  xyz_dir_fn = fullfile(project_dir,'grids',strcat(output_xyz_fn,'.txt'));
end
if fl_conv
  fl_dir_fn = fullfile(project_dir,'flightlines',strcat(output_flight_fn,'.txt'));
end
fprintf('%s', datestr(now,'HH:MM:SS\n'));

%% Write XYZGrid File
if ascii_conv
  % Load in ESRI ASCII Grids
  fprintf('Loading in the ESRI ASCII Grids ... ');
  [SURF SR] = arcgridread(surface_fn);
  [BED BR] = arcgridread(bottom_fn);
  [THICK TR] = arcgridread(thickness_fn);
  fprintf('%s', datestr(now,'HH:MM:SS\n'));
  
  % Get the X/Y (LAT/LON) Values from the grids
  fprintf('Getting LAT/LON Values ... ');
  [ROWS COLS] = size(SURF);
  LAT = zeros(ROWS,COLS);
  LON = zeros(ROWS,COLS);
  for idxR = 1:ROWS
    for idxC = 1:COLS
      [LAT(idxR,idxC) LON(idxR,idxC)] = pix2latlon(SR,idxR,idxC);
    end
  end
  fprintf('%s', datestr(now,'HH:MM:SS\n'));
  
  % Write the values to an XYZ TXT File
  fprintf('Writing XYZGrid output file ... ');
  %Convert to column vectors
  LAT=LAT(:);LON=LON(:);THICK=THICK(:);SURF=SURF(:);BED=BED(:);
  %Change NaN to -9999
  SURF(isnan(SURF))=-9999;THICK(isnan(THICK))=-9999;BED(isnan(BED))=-9999;
  
  % Save file
  fid_xyz = fopen(xyz_dir_fn,'w');
  fprintf(fid_xyz,'%s,%s,%s,%s,%s\n','LAT','LON','THICK','SURFACE','BOTTOM');
  for idxD = 1:length(LAT)
    fprintf(fid_xyz,'%2.6f,%2.6f,%6.2f,%6.2f,%6.2f\n',LAT(idxD),LON(idxD),THICK(idxD),SURF(idxD),BED(idxD));
  end
  fprintf('%s', datestr(now,'HH:MM:SS\n'));
  fcid = fclose(fid_xyz);
  if fcid == 0
  else
    fprintf('TXT File Release Error.\n');
  end
end
%% Write Flightline TXT File
if fl_conv
  fprintf('Writing flightline TXT from DBF ... ');
  dbf_to_txt(flightline_dbf_fn,fl_dir_fn);
  fprintf('%s', datestr(now,'HH:MM:SS\n'));
end
%% Complete Script
fprintf('Post GIS to TXT Complete ... %s\n',datestr(now,'HH:MM:SS\n'));
