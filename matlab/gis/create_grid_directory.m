function create_grid_directory(output_dir,folder_name)
% Creates an empty Level-III Project Directory
% Used in Level-III Data Product Creation
%
% Format: create_grid_directory(output_dir,folder_name);
%
% ----- Example -----
% fout = 'C:\Users\SomeUser\Projects\'; 
% fname = 'NWCoast_2010_2011_Composite';
% create_grid_directory(fout,fname);
%
% Above example creates the following in "fout"
%
% NWCoast_2010_2011_Composite/
%   NWCoast_2010_2011_Composite/
%   errrors/
%   boundaries/
%   flightlines/
%   grids/
%    >>   ascii/
%   POST/
%    >>   errors/
%    >>   boundaries/
%    >>   flightlines/
%    >>   grids/
%    >>   preview_images/
%    >>   readme/
%
% Author: Kyle W. Purdon, UGRA

% Print Start Block
fprintf('Creating Directory Structure ... %s\n',datestr(now,13));

% Save the current path
cpath = pwd;

% Make a valid path
out_dir = fullfile(output_dir,filesep);

% Check if the output path exists
if ~exist(out_dir,'dir')
  mkdir(out_dir);
  cd(out_dir);
end

% Make the root folder
root_dir = fullfile(out_dir,folder_name,filesep);

% Check for a pre-existing folder
if ~exist(root_dir,'dir')
  mkdir(root_dir);
  cd(root_dir);
else
  fprintf('The directory already exists. (Delete or Stop?)\n');
  while true
    runstop = input('Delete (1) or Stop (2): ');
    if runstop == 2
      clear param;
      error('Program Terminated By User !!!');
    elseif runstop == 1
      rmdir(root_dir,'s');
      mkdir(root_dir);
      cd(root_dir);
      break;
    else
      fprintf('Invalid Choice. Choose value 1 or 2!\n');
    end
  end
end

% Make first level of sub-dirs
bound_dir = fullfile(root_dir,'boundaries',filesep);
mkdir(bound_dir);
error_dir = fullfile(root_dir,'errors',filesep);
mkdir(error_dir);
fl_dir = fullfile(root_dir,'flightlines',filesep); 
mkdir(fl_dir);
grid_dir = fullfile(root_dir,'grids',filesep);
mkdir(grid_dir);
post_dir = fullfile(root_dir,'post',filesep);
mkdir(post_dir);
final_dir = fullfile(root_dir,folder_name,filesep);
mkdir(final_dir);

% Make the second level of sub_dirs
  cd(grid_dir);
  mkdir(pwd,'ascii');
  cd(post_dir);
  mkdir(pwd,'boundaries');
  mkdir(pwd,'flightlines');
  mkdir(pwd,'preview_images');
  mkdir(pwd,'readme');
  mkdir(pwd,'grids');
  mkdir(pwd,'errors');

% Write the Standards File
cd(root_dir);
fn = 'FILE_NAMING_STANDARDS.txt';
fprintf('Writing Standards File ... %s\n',datestr(now,13));
fid = fopen(fn,'wt');
if fid == -1
  error('Error writing standards file.');
end
fprintf(fid,'CReSIS Grid Project File Structure and Naming Standards\n');
fprintf(fid,'____________________________________________________________________________\n');
fprintf(fid,'\n');

fprintf(fid,'Naming Structure:\n');
fprintf(fid,'>> GlacierName_StartYear_EndYear_Composite_DataType\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Flightlines)\n');
fprintf(fid,'\n');
fprintf(fid,'____________________________________________________________________________\n');
fprintf(fid,'POST directory file structure:\n');
fprintf(fid,'\n');
fprintf(fid,'-- grids --\n');
fprintf(fid,'\n');
fprintf(fid,'ascii\n');
fprintf(fid,'>> 3 ASCII Grids of Surface,Bed,Thickness and 1 TXT XYZ file\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Surface.txt)\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Bottom.txt)\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Thickness.txt)\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_XYZGrid.txt)\n');
fprintf(fid,'\n');
fprintf(fid,'-- errors --\n');
fprintf(fid,'\n');
fprintf(fid,'>> 1 TXT file, 1 MAT, 1 TXT, and 1 PNG of the crossovers \n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Crossovers.csv\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Crossovers.mat\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Crossovers.png\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Crossovers_Errors.txt\n');
fprintf(fid,'\n');
fprintf(fid,'-- flightlines --\n');
fprintf(fid,'\n');
fprintf(fid,'>> 1 TXT file and 1 ESRI shapefile of the flightlines \n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Flightlines.shp\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Flightlines.txt\n');
fprintf(fid,'\n');
fprintf(fid,'-- boundaries --\n');
fprintf(fid,'\n');
fprintf(fid,'>>1 shapefile of the studyarea.\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_StudyArea.shp\n');
fprintf(fid,'\n');
fprintf(fid,'-- preview_images --\n');
fprintf(fid,'\n');
fprintf(fid,'>> 4 PNG maps of the Surface,Bed,Thickness, and Flightlines\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Surface_Preview.png)\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Bottom_Preview.png)\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Thickness_Preview.png)\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_Flightlines_Preview.png)\n');
fprintf(fid,'\n');
fprintf(fid,'-- readme --\n');
fprintf(fid,'\n');
fprintf(fid,'>> 1 .txt README file created using grid_readme_creator.m\n');
fprintf(fid,'>> Ex.(PineIsland_2008_2010_Composite_README.txt)\n');
fprintf(fid,'\n');
fprintf(fid,'____________________________________________________________________________\n');
fprintf(fid,'\n');
fprintf(fid,'OTHER Considerations:\n');
fprintf(fid,'\n');
fprintf(fid,'>> All Shapefiles should include 6+ files. This can be accomplished\n');
fprintf(fid,'   using ArcMap to export the file to the post directory. The file types\n');
fprintf(fid,'   needed are (SHP,SHX,SBN,SBX,DBF,PRJ)\n');
fprintf(fid,'>> All ASCII grids should include 2+ files. This can be accomplished\n');
fprintf(fid,'   using ArcMap "Raster to ASCII" to export the files to the POST directory.\n');
fprintf(fid,'   The files needed are (TXT,PRJ)\n');
cid = fclose(fid);
if cid == -1
  fprintf('Error Closing File.\n');
end

% Print completeion and return to original dir.
cd(cpath);
fprintf('  Done ... %s\n',datestr(now,13));
end