%% User Settings
% =========================================================================

param_override = [];
param = [];
passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});

% master_pass_idx = struct('day_seg','20140505_01','frms',17:18,'param_fn','rds_param_2014_Greenland_P3.xls','in_path','CSARP_post/standard');
% master_pass_idx = struct('day_seg','20170331_01','frms',16:17,'param_fn','rds_param_2017_Greenland_P3.xls','in_path','CSARP_post/standard');

%% 2002 (not processed yet)
% param_fn = 'rds_param_2002_Greenland_P3.xls';
% pass_name = sprintf('2002_Greenland_P3');
% passes(end+1) = struct('frm','20020528_06_005','wf_adc',[1 1],'param_fn',param_fn);
% passes(end+1) = struct('frm','20020528_06_006','wf_adc',[1 1],'param_fn',param_fn);
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

%% Petermann Line 1 2014, 2018
pass_name = sprintf('Petermann_line1_2014_2018');
dist_min = 500;
master_pass_idx = 1;
start = struct('lat',80.503546,'lon',-60.034880);
stop = struct('lat',80.964739,'lon',-61.618307);

param_fn = 'rds_param_2014_Greenland_P3.xls';
day_seg = '20140505_01';
frms = 17:18;
passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

param_fn = 'rds_param_2018_Greenland_P3.xls';
day_seg = '20180405_01';
frms = 15:16;
passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');

%% loading
% fn_dir_2011 = 'X:\ct_data\rds\2011_Greenland_P3\CSARP_post\CSARP_csarp-combined\';
% addpath(fn_dir_2011);
% fn_dir_2014 = 'X:\ct_data\rds\2014_Greenland_P3\CSARP_post\CSARP_standard\';
% addpath(fn_dir_2014);
% fn_dir_2017 = 'X:\ct_data\rds\2017_Greenland_P3\CSARP_post\CSARP_standard\';
% addpath(fn_dir_2017);

%2011
% for passes_idx = 1:length(param_fn)
%   if passes.param_fn = 'rds_param_2011_Greenland_P3.xls'
%     load(passes.frm)
%   end 
%   mdata = load(passes.frm(1:3));
% end
% 
% %2011
% for passes_ = 'rds_param_2011_Greenland_P3.xls'
%   fn_dir_2011 = 'X:\ct_data\rds\2011_Greenland_P3\CSARP_post\CSARP_csarp-combined\';
%   mdata = load(passes.frm(1:3));
% end
% 
% %%
% %2014
% for passes_idx = 1:length(passes)
%  mdata = load(passes.frm);
% end
% %2017
% for passes_idx = 1:length(passes)
% mdata = load(passes.frm);
% 
% mdata = load('passes.param_fn');
% 
%   mdata=load('Data_20170413_01_056')
%   dist_min = 2000;
% 
%   mdata=load('Data_20170413_01_056')
%  fcs = mdata.param_combine.combine.fcs{1}{1};
% 
% end
% 
% param.dist_min = dist_min;

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
