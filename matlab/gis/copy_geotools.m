function copy_geotools(userpath,toolboxpath)
% Copy's all the CReSIS GeoTools from a specified CReSIS-Toolbox
%
% Format: copy_geotools(userpath,toolboxpath)
%
% Supports absolute or relative paths (Uses copyfile.m)
%
% ----- Examples -----
%
% Windows Command Copy
% copy_geotools('H:\scripts\matlab\','H:\scripts\cresis-toolbox\');
%
% UNIX CommandCopy
% copy_geotools('/users/kpurdon/scripts/matlab','/users/kpurdon/scripts/cresis-toolbox/');
%
% Author: Kyle W. Purdon, UGRA

% Print Start Block
fprintf('Copying GeoTools from Toolbox\n');

% Check to make sure the user path exists, if not create it.
if ~exist(userpath,'dir')
  mkdir(userpath);
end

% Set locations of GeoTools
file_locations = {};
file_locations{end+1} = fullfile(toolboxpath,'posting','run_cross_over_analysis.m');
file_locations{end+1} = fullfile(toolboxpath,'gis','post_gis_to_txt.m');
file_locations{end+1} = fullfile(toolboxpath,'gis','grid_readme_creator.m');
file_locations{end+1} = fullfile(toolboxpath,'gis','calc_ideal_cellsize.m');
file_locations{end+1} = fullfile(toolboxpath,'gis','run_gis_dataprep.m');
  
% Copy all files over
for file_idx = 1:length(file_locations);
  copyfile(file_locations{file_idx},userpath,'f');
end

% Print Complete
fprintf('  File Copy Complete.\n');

end