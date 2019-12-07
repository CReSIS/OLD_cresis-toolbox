%% User Settings
% =========================================================================

param_override = [];
param = [];
passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});

% master_pass_idx = struct('day_seg','20140505_01','frms',17:18,'param_fn','rds_param_2014_Greenland_P3.xls','in_path','CSARP_post/standard');
% master_pass_idx = struct('day_seg','20170331_01','frms',16:17,'param_fn','rds_param_2017_Greenland_P3.xls','in_path','CSARP_post/standard');

%% Petermann Line 1 2002
% pass_name = sprintf('Petermann_line1_2002');
%dist)min = 500;
%master_pass_idx = 1

% param_fn = 'rds_param_2002_Greenland_P3.xls';
%day_seg = '20020528_06';
%frms = 5:6;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
% start = struct('lat',80.492292,'lon',-59.975190); 
% stop = struct('lat',81.116866,'lon',-62.041787); 

%% Petermann Line 1 2011
% pass_name = sprintf('Petermann_line1_2011');
% dist_min = 500;
% master_pass_idx = 1;

% param_fn = 'rds_param_2011_Greenland_P3.xls';
% day_seg = '20110507_02';
% frms = 9:11;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/csarp-combined');
% start = struct('lat',80.499370,'lon',-60.013440); 
% stop = struct('lat',80.952663,'lon',-61.586491); 

%% Petermann Line 1 2014
% pass_name = sprintf('Petermann_line1_2014');
% dist_min = 500;
% master_pass_idx = 1;
% start = struct('lat',80.503546,'lon',-60.034880);
% stop = struct('lat',80.929117,'lon',-61.524422); 
% 
% param_fn = 'rds_param_2014_Greenland_P3.xls';
% day_seg = '20140505_01';
% frms = 17:18;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

%% 2017
% param_fn = 'rds_param_2017_Greenland_P3.xls';
% day_seg = '20170331_01';
% frms = 16:17;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
% start = struct('lat',80.504179,'lon',-60.037068);
% stop = struct('lat',80.962162,'lon',-61.612482);
% pass_name = sprintf('2017_Greenland_P3');

%% Petermann Line 1 2018
% pass_name = sprintf('Petermann_line1_2018');
% dist_min = 500;
% master_pass_idx = 1;
% start = struct('lat',80.521230,'lon',-60.127320);
% stop = struct('lat',80.964739,'lon',-61.618307);
% 
% param_fn = 'rds_param_2018_Greenland_P3.xls';
% day_seg = '20180405_01';
% frms = 15:16;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

%% Petermann Line 1 2011, 2014, 2018
% pass_name = sprintf('Petermann_line1_2011_2014_2018');
% dist_min = 500;
% master_pass_idx = 2;
% start = struct('lat',80.499370,'lon',-60.013440);
% stop = struct('lat',80.964739,'lon',-61.618307);
% 
% param_fn = 'rds_param_2011_Greenland_P3.xls';
% day_seg = '20110507_02';
% frms = 9:11;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
% 
% param_fn = 'rds_param_2014_Greenland_P3.xls';
% day_seg = '20140505_01';
% frms = 17:18;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
% 
% param_fn = 'rds_param_2018_Greenland_P3.xls';
% day_seg = '20180405_01';
% frms = 15:16;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

%% Petermann Line 2 2013, 2014
% pass_name = sprintf('Petermann_line2_2013_2014');
% dist_min = 500;
% master_pass_idx = 2;
% start = struct('lat',80.517191,'lon',-59.844969);
% stop = struct('lat',80.903411,'lon',-61.332159);
% 
% param_fn = 'rds_param_2013_Greenland_P3.xls';
% day_seg = '20130420_02';
% frms = 4:5;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
% 
% param_fn = 'rds_param_2014_Greenland_P3.xls';
% day_seg = '20140512_01';
% frms = 16:17;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

%% Petermann Line 4 2011, 2014, 2018
% pass_name = sprintf('Petermann_line4_2010_2011_2013_2014');
% dist_min = 500;
% master_pass_idx = 4;
% start = struct('lat',80.5295278,'lon',-59.567146);
% stop = struct('lat',80.995874,'lon',-61.357055);
% 
% param_fn = 'rds_param_2010_Greenland_DC8.xls';
% day_seg = '20100324_01';
% frms = 25:26;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
% 
% param_fn = 'rds_param_2011_Greenland_P3.xls';
% day_seg = '20110507_02';
% frms = 13:15;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
% 
% param_fn = 'rds_param_2013_Greenland_P3.xls';
% day_seg = '20130420_02';
% frms = 7:8;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
% 
% param_fn = 'rds_param_2014_Greenland_P3.xls';
% day_seg = '20140512_01';
% frms = 13:14;
% passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

%% 79N Line 1 2010, 2014, 2016, 2018
pass_name = sprintf('79N_line1_2010_2014_2016_2018');
dist_min = 16000;
master_pass_idx = 2;
start = struct('lat',79.346109,'lon',-22.575643);
stop = struct('lat',79.559002,'lon',-19.329911);

param_fn = 'rds_param_2010_Greenland_P3.xls';
day_seg = '20100525_04';
frms = 11:13;
passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

param_fn = 'rds_param_2014_Greenland_P3.xls';
day_seg = '20140429_01';
frms = 43:44;
passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

param_fn = 'rds_param_2016_Greenland_P3.xls';
day_seg = '20160509_10';
frms = 1;
passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

param_fn = 'rds_param_2018_Greenland_P3.xls';
day_seg = '20180418_05';
frms = 1:2;
passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

param.combine_passes.pass_name = pass_name;
param.combine_passes.passes = passes;
param.combine_passes.master_pass_idx = master_pass_idx;
param.combine_passes.start = start;
param.combine_passes.stop = stop;
param.combine_passes.dist_min = dist_min;
param.day_seg = passes(master_pass_idx).day_seg;

%% Automated section
% =========================================================================

%read radar names
global gRadar;

% Input checking
if exist('param_override','var')
  param_override = merge_structs(gRadar,param_override);
else
  param_override = gRadar;
end

% multipass.combine_passes(param, param_override);
multipass.combine_passes;
