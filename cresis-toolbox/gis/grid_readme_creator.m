%% Script Information/Details

% =============================================================
% About: Creates a formatted TXT README file from custom user input.
% Author: Kyle Purdon, Weibo Liu
% Version Date: 10/04/2014
% Contact Information: cresis_data@cresis.ku.edu
% Contact Subject: "Grid Readme Creator"
% =============================================================

%% User Input
% Output folder
param.outputFolder = 'C:\Users\SomeUser\SomeProjectFolder\post\readme\';

% Base Name of Project Files
% Format: GlacierName_StartYear_EndYear_Composite
param.basename = 'GlacierName_StartYear_EndYear_Composite';

% Cell Size used in interpolation (meters)
% For L3 Data this should be 500
param.finalCellSize = 500;

% Location of Project
% (1)=Greenland (2)=Antarctica (3)=Canada
param.Location = 2;

% Earliest Year In Dataset
param.EARLYYEAR = 2011;

% Optional Data Citations 
% Turn citation prints on with param.citations;
param.citations = true;
% Select individual citations.
param.ATMLIDAR = true; % ** ONLY USE IF YOUR DATA CONTAINS ATM LIDAR (Check "DataType") **
param.ICEFREE = true; % Shows both ICESat and IceMask Citations

%% Automated Section

%Print Opening Header

fprintf('================================\n');
fprintf('         README Creator         \n');
fprintf('================================\n');

%Open new file for writing
param.outputFolder = fullfile(param.outputFolder,filesep);
fileID = strcat(param.outputFolder,param.basename,'_README.txt');
fid = fopen(fileID,'w+');
if fid == -1
  error('Error Opening File. Check path and filename');
else
  fprintf('File Opened. Writing In Progress...');
end

%Create values needed
switch param.Location
  case 1
    param.locationvalue = 'Greenland';
    param.proj = 'WGS_84_NSIDC_Sea_Ice_Polar_Sterographic_North';
  case 2
    param.locationvalue = 'Antarctica';
    param.proj = 'WGS_84_Antarctic_Polar_Stereographic';
  case 3
    param.locationvalue = 'Canada';
    param.proj = 'WGS_84_NSIDC_Sea_Ice_Polar_Stereographic_North';
end

%Write to file
fprintf(fid,'CReSIS Gridded Data Product README \n');
fprintf(fid,'-----------------------------------------------------------------\n\n');
fprintf(fid,'Version Date: %s \n',datestr(now,1));
fprintf(fid,'Contact: cresis_data@cresis.ku.edu \n\n');
fprintf(fid,'Location: %s \n',param.locationvalue);
fprintf(fid,'GlacierName: %s \n',param.basename);
fprintf(fid,'Projection: %s \n\n',param.proj);
fprintf(fid,'Datum: WGS_84 \n\n');
fprintf(fid,'SoftwareUsed: ESRI ArcGIS 10.1 (ArcInfo) \n');
fprintf(fid,'SoftwareUsed: MATLAB R2014a \n\n');
fprintf(fid,'ProductsIncluded: Flightlines \n');
fprintf(fid,'ProductsIncluded: Errors \n');
fprintf(fid,'ProductsIncluded: Boundaries \n');
fprintf(fid,'ProductsIncluded: Grids \n');
fprintf(fid,'ProductsIncluded: PreviewImages \n\n');

fprintf(fid,'Flightline Details: \n');
fprintf(fid,'-----------------------------------\n');
fprintf(fid,'- ESRI Shapefile and TXT of all data points used \n');
fprintf(fid,'- A_SURF and A_BED are the variables used to interpolate. \n');
fprintf(fid,'- NASA ATM data is used, if it exists. See variable "Data_Type" \n');
fprintf(fid,'- ICESat data is used, if it exists for IceFree Areas. See variable "Data_Type" \n');
fprintf(fid,'- Flightlines are clipped to a 10km buffer of the boundary (Study Area) \n\n');

fprintf(fid,'Error Details: \n');
fprintf(fid,'-----------------------------------\n');
fprintf(fid,'- Results of Crossover Analysis \n');
fprintf(fid,'- Note: Crossovers Are corrected season by season, not across season.\n');
fprintf(fid,'        Errors present are across, not within season.\n');
fprintf(fid,'- CSV File of Crossover Analysis Results \n');
fprintf(fid,'- MAT File of Crossover Analysis Results \n');
fprintf(fid,'- TXT File of Statistics from Crossover Analysis Results \n');
fprintf(fid,'- PNG Preview Image of Crossover Analysis Results \n\n');

fprintf(fid,'Boundary Details: \n');
fprintf(fid,'-----------------------------------\n');
fprintf(fid,'- ESRI Shapefile of "Study Area" extent \n\n');

fprintf(fid,'Grid Details: \n');
fprintf(fid,'-----------------------------------\n');
fprintf(fid,'- Surface interpolated using IDW interpolation in ArcGIS \n');
fprintf(fid,'- Bottom interpolated using TopoToRaster interpolation in ArcGIS \n');
fprintf(fid,'   Note: Any bottom past the grounding line represents ice bottom not ocean bottom. \n');
fprintf(fid,'- Thickness is calculated using Surface minus Bottom (above) in ArcGIS \n');
fprintf(fid,'- ASCII rasters for Surface,Thickness,and Bed Elevation are provided \n');
fprintf(fid,'- An XYZ TXT File containing data from all grids is also provided \n\n');
fprintf(fid,'GridNoDataValue: -9999 \n');
fprintf(fid,'GridCellSize: %d x %dm \n\n',param.finalCellSize,param.finalCellSize);

fprintf(fid,'PreviewImage Details: \n');
fprintf(fid,'-----------------------------------\n');
fprintf(fid,'- PNG maps of Flightlines,Surface,Thickness,and Bed Elevation \n\n');

fprintf(fid,'Other Information: \n');
fprintf(fid,'-----------------------------------\n');
fprintf(fid,'- If NASA ATM data exists for a season, ATM Surface is used. (Data_Type) \n');
fprintf(fid,'- To read ASCII grids use Conversion > ASCIItoRaster in ArcGIS. \n');
fprintf(fid,'- To read ASCII grids use arcgridread.m in MATLAB(Mapping). \n');
fprintf(fid,'- To read SHP files use shaperead.m in MATLAB(Mapping) \n\n');

fprintf(fid,'Citing and Acknowledgements : \n');
fprintf(fid,'-----------------------------------\n');
if param.citations
  fprintf(fid,'External Data Used:\n\n');
  % NASA ATM LIDAR CITATION
  if param.ATMLIDAR
    if param.EARLYYEAR <= 2002;
      fprintf(fid,'NASA ATM LIDAR\n');
      fprintf(fid,'IceBridge ATM L2 Icessn Elevation, Slope, and Roughness\n');
      fprintf(fid,'Pre-IceBridge ATM L2 Icessn Elevation, Slope, and Roughness\n');
      fprintf(fid,'http://nsidc.org/data/docs/daac/icebridge/ilatm2/index.html\n\n');
    else
      fprintf(fid,'NASA ATM LIDAR\n');
      fprintf(fid,'IceBridge ATM L2 Icessn Elevation, Slope, and Roughness\n');
      fprintf(fid,'http://nsidc.org/data/docs/daac/icebridge/ilatm2/index.html\n\n');
    end
  end
  % ICE FREE DATA CITATION
  if param.ICEFREE;
    if param.Location == 1 %Greenland
      fprintf(fid,'ICESAT DEM\n');
      fprintf(fid,'DiMarzio, J., A. Brenner, R. Schutz, C. A. Shuman, and H. J. Zwally. 2007. GLAS/ICESat 1 km laser altimetry digital elevation model of Greenland. Boulder, Colorado USA: National Snow and Ice Data Center. Digital media.\n');
      fprintf(fid,'IceFree Mask\n');
      fprintf(fid,'Howat I.M and A. Negrete, in prep, A high-resolution ice mask for the Greenland Ice Sheet and peripheral glaciers and icecaps.\n\n');
    elseif param.Location == 2 %Antarctica
      fprintf(fid,'ICESAT DEM\n');
      fprintf(fid,'DiMarzio, J., A. Brenner, R. Schutz, C. A. Shuman, and H. J. Zwally. 2007. GLAS/ICESat 500 m laser altimetry digital elevation model of Antarctica. Boulder, Colorado USA: National Snow and Ice Data Center. Digital media.\n');
      fprintf(fid,'IceFree Mask\n');
      fprintf(fid,'The Antarctic Digital Database, Scientific Committtee on Antarctic Research 1993 - 2006. E00 Rock Outcrop (Lines and Polygons)\n\n');
    end
  end
end
fprintf(fid,'Whenever the data are used, please include the following acknowledgement:\n\n');
fprintf(fid,'We acknowledge the use of data and/or data products from CReSIS generated with\n');
fprintf(fid,'support from NSF grant ANT-0424589 and NASA grant NNX10AT68G\n\n');
fprintf(fid,'Please cite data according to NSIDC standard\n\n');

%Close file
fprintf('%s\n',datestr(now,13));
cid = fclose(fid);
if cid == -1
  fprintf('Error Closing File.\n');
else
  fprintf('File Close Succesfull\n');
end

fprintf('README Creation Complete ... %s \n',datestr(now,13));
