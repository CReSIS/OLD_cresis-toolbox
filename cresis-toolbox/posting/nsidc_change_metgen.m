function nsidc_change_metgen(in_fn, version_id, premet_spatial_dir, data_dir, radar_type, data_type)
% nsidc_change_metgen(in_fn, version_id, premet_spatial_dir, data_dir, radar_type, data_type)
%
% Creates the MetGen configuration file for MetGen.Start automatically.
%
% in_fn: The path to MetGen tool
%   e.g. /cresis/snfs1/scratch1/paden/nsidc/2013_Greenland_P3/rds
% version_id: MCF version ID (currently should be string set to '002')
% premet_spatial_dir: the premet_spatial filename
%   e.g. /cresis/snfs1/scratch1/paden/nsidc/2013_Greenland_P3/rds/IRMCR2_Files/Spatial_Premet_Files/
% data_dir: the output data directory
%   e.g. /cresis/snfs1/scratch1/paden/nsidc/2013_Greenland_P3/rds/IRMCR2_Files/Data_2013_GR
% radar_type: String containing NSIDC radar type which should be one of:
%   IRACC, IRKUB, IRMCR, or IRSNO
% data_type: String containing NSIDC data type which should be '1B' or '2'
%
% See also: type "nsidc_help.m"
%
% Author: Yi Zhu, John Paden

if strcmpi(radar_type,'IRKUB') || strcmpi(radar_type,'IRSNO')
  sidelength = 1000;
else
  sidelength = 200;
end

if strcmpi(data_type,'2')

  %% Changes made into IRMCR2_PrimaryConfig.cfg

  pri_cfg_fn = fullfile(in_fn, '/SIPSMetGen/Examples/IRMCR2_SimpleExample/SIPSMetGenTool/IRMCR2_Configs/IRMCR2_PrimaryConfig.cfg');

  [fid,msg] = fopen(pri_cfg_fn, 'w+');

  if fid < 1
    fprintf('Could not open file %s\n', pri_cfg_fn);
    error(msg);
  end

  file_Ext = '.csv';

  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '#          Required Files          \n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '#***** File Directories\n');
  fprintf(fid, '# PreMet_Dir - directory where SIPSMetGen tool will find the .premet files for each data file\n');
  fprintf(fid, 'PreMet_Dir=%s/\n', premet_spatial_dir);
  fprintf(fid, '# Spatial_Dir - directory where SIPSMetGen tool will find the .spatial files for each data file\n');
  fprintf(fid, 'Spatial_Dir=%s/\n', premet_spatial_dir');
  fprintf(fid, '\n');
  fprintf(fid, '#***** Provided Files\n');
  fprintf(fid, '# Path to the MCF file that is delivered for each data product\n');
  fprintf(fid, 'MCF_File=%s\n', fullfile(in_fn, '/SIPSMetGen/Examples/IRMCR2_SimpleExample/SIPSMetGenTool/IRMCR2_Configs/', ...
    sprintf('%s%s#%s.MCF', radar_type, data_type, version_id)));
  fprintf(fid, '# Path to valids file which contains list of complete list of valid meta data values\n');
  fprintf(fid, 'Valids_File=%s\n', fullfile(in_fn, '/SIPSMetGen/Examples/IRMCR2_SimpleExample/SIPSMetGenTool/IRMCR2_Configs/valids.cfg'));
  fprintf(fid, '\n');
  fprintf(fid, '#***** Push Files\n');
  fprintf(fid, '# Where to find the data push config if using tool to send data to NSIDC with -P or -PO options\n');
  fprintf(fid, 'Push_File=%s\n', fullfile(in_fn, '/SIPSMetGen/Examples/IRMCR2_SimpleExample/SIPSMetGenTool/IRMCR2_Configs/push.cfg'));
  fprintf(fid, '\n');
  fprintf(fid, '#***** Data Files\n');
  fprintf(fid, '# Path to data files to be processed\n');
  fprintf(fid, 'Data_Dir=%s/\n', data_dir);
  fprintf(fid, '# Extension to look for to identify files as a data file to process.\n');
  fprintf(fid, 'Data_Ext=%s\n', file_Ext);
  fprintf(fid, '\n');
  fprintf(fid, '#***** Output Directories\n');
  fprintf(fid, '# Path to stage output PDRs\n'); 
  fprintf(fid, 'PDR_Dir=%s\n', fullfile(in_fn, sprintf('%s%s_Files', ...
    radar_type, data_type), '/output/pdrs/'));
  fprintf(fid, '# Path to stage output Metadata files\n');
  fprintf(fid, 'MET_Dir=%s\n', fullfile(in_fn, sprintf('%s%s_Files', ...
    radar_type, data_type), '/output/mets/'));
  fprintf(fid, '\n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '#          Configuration Parameters          \n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '# Provider version id for tracking changes to data - integer values only\n');
  fprintf(fid, 'VersionID_local=1\n');
  fprintf(fid, '# Cell size parameter for starting values for spatial algorithm\n');
  fprintf(fid,  sprintf('SideLength=%d\n',sidelength));
  fprintf(fid, '# Determines the max number of cells for the spatial algorithm for creating geo-polygon cell sizes will scale up until cell count falls below this max value\n');
  fprintf(fid, 'MaxCellCount=300000\n');
  fclose(fid);
   
elseif strcmpi(data_type,'1B')

  %% Changes made into ILATM1B_MultiBrowseExamplePrimary.cfg

  pri_cfg_fn = fullfile(in_fn, '/SIPSMetGen/Examples/ILATM1B_MultiDataBrowseExample/SIPSMetGenTool/ILATM1B_Configs/ILATM1B_MultiBrowseExamplePrimary.cfg');

  fid= fopen(pri_cfg_fn, 'w+');

  if fid < 1
    fprintf('Could not open file %s\n', pri_cfg_fn);
    error(msg);
  end

  file_Ext = '.nc _Map.jpg _Echogram.jpg _Echogram_Picks.jpg';

  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '#          Required Files          \n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '#***** File Directories\n');
  fprintf(fid, 'PreMet_Dir=%s/\n', premet_spatial_dir);
  fprintf(fid, 'Spatial_Dir=%s/\n', premet_spatial_dir');
  fprintf(fid, '\n');
  fprintf(fid, '#***** Provided Files\n');
  fprintf(fid, '# Path to the MCF file that is delivered for each data product\n');
  fprintf(fid, 'MCF_File=%s\n', fullfile(in_fn, '/SIPSMetGen/Examples/ILATM1B_MultiDataBrowseExample/SIPSMetGenTool/ILATM1B_Configs/', ...
    sprintf('%s%s#%s.MCF',radar_type, data_type, version_id)));
  fprintf(fid, '# Path to valids file which contains list of complete list of valid meta data values\n');
  fprintf(fid, 'Valids_File=%s\n', fullfile(in_fn, '/SIPSMetGen/Examples/ILATM1B_MultiDataBrowseExample/SIPSMetGenTool/ILATM1B_Configs/valids.cfg'));
  fprintf(fid, '\n');
  fprintf(fid, '#***** Push Files\n');
  fprintf(fid, 'Push_File=%s\n', fullfile(in_fn, '/SIPSMetGen/Examples/ILATM1B_MultiDataBrowseExample/SIPSMetGenTool/ILATM1B_Configs/push.cfg'));
  fprintf(fid, '\n');
  fprintf(fid, '#***** Data Files\n');
  fprintf(fid, '# Path to data files to be processed\n');
  fprintf(fid, 'Data_Dir=%s/\n', data_dir);
  fprintf(fid, '# Extension to look for to identify files as a data file to process.\n');
  fprintf(fid, 'Data_Ext=%s\n', file_Ext);
  fprintf(fid, '\n');
  fprintf(fid, '#***** Output Directories\n');
  fprintf(fid, '# Path to stage output PDRs\n');
  fprintf(fid, 'PDR_Dir=%s\n', fullfile(in_fn, sprintf('%s%s_Files', ...
    radar_type,data_type), '/output/pdrs/'));
  fprintf(fid, '# Path to stage output Metadata files\n');
  fprintf(fid, 'MET_Dir=%s\n', fullfile(in_fn, sprintf('%s%s_Files', ...
    radar_type,data_type), '/output/mets/'));
  fprintf(fid, '\n');
  fprintf(fid, '\n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '#          Configuration Parameters          \n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '# Provider version id for tracking changes to data - integer values only\n');
  fprintf(fid, 'VersionID_local=1\n');
  fprintf(fid, '# Cell size parameter for starting values for spatial algorithm\n');
  fprintf(fid,  sprintf('SideLength=%d\n',sidelength));
  fprintf(fid, '# Determines the max number of cells for the spatial algorithm for creating geo-polygon cell sizes will scale up until cell count falls below this max value\n');
  fprintf(fid, 'MaxCellCount=300000\n');
  fprintf(fid, '\n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '#          Optional Parameters          \n');
  fprintf(fid, '#++++++++++++++++++++++++++++++++++++\n');
  fprintf(fid, '# The file extensions are just examples\n');
  fprintf(fid, '# to show the format of multiple file\n');
  fprintf(fid, '# extensions. The file extensions could\n');
  fprintf(fid, '# be anything, as long as they are\n');
  fprintf(fid, '# delimited by a space and begin with\n');
  fprintf(fid, '# a ".".\n');
  fprintf(fid, '#***** Browse Files\n');
  fprintf(fid, '#Browse_Req=N\n');
  fprintf(fid, '#Browse_Dir=%s\n', data_dir);
  fprintf(fid, '#Browse_Ext=.brws\n');
  fprintf(fid, '##### NOTE Browse_Req=N will add browse for each pair of data files if there is a .brws file with matching name\n');
  fprintf(fid, '###### Browse_Req=Y will fail to create metadata and PDR files for data if a matching .brws file does not exist.\n');
  fprintf(fid, '#***** QA Files\n');
  fprintf(fid, '#QA_Req=N\n');
  fprintf(fid, '#QA_Dir=/home/SIPSMetGen/qa/\n');
  fprintf(fid, '#QA_Ext=.QA0 .QA1 .QA2\n');
  fprintf(fid, '#***** PH Files\n');
  fprintf(fid, '#PH_Req=N\n');
  fprintf(fid, '#PH_Dir=/home/SIPSMetGen/ph/\n');
  fprintf(fid, '#PH_Ext=.PH0 .PH1 .PH2\n');
  fprintf(fid, '#\n');
  fprintf(fid, '##################################');
  
  fclose(fid);
  
else 
  error('Unsupported data_type %s (must be 1B or 2)\n', data_type);
end

end

