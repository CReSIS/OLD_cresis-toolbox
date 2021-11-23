% script run_combine_passes
%
% Script for running combine_passes.m
%
% Authors: Cody Barnett, Bailey Miller, John Paden
%
% See also: multipass.combine_passes.m, multipass.run_combine_passes.m,
% multipass.multipass.m, multipass.run_multipass.m

%% User Settings
% =========================================================================

param_override = [];
param = [];
passes = [];

example_str = 'egig_2011_2012_2014_2018_allwf';

if strcmpi(example_str,'Thwaites_201902_201912_202001')
  %% Thwaites Line 1 20190201_01, 20191225_01, 20200127_01
  pass_name = sprintf('Thwaites_201902_201912_202001');
  dist_min = 300;
  master_pass_idx = 1;
  start = struct('lat',-75.137955,'lon',-105.538244);
  stop = struct('lat',-75.259933,'lon',-105.418265);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
  param_fn = 'accum_param_2018_Antarctica_TObas.xls';
  day_seg = '20190201_01';
  frms = 33:34;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'accum_param_2019_Antarctica_TObas.xls';
  day_seg = '20191225_01';
  frms = 20:22;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'accum_param_2019_Antarctica_TObas.xls';
  day_seg = '20200127_01';
  frms = 33:34;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
end

if strcmpi(example_str,'Petermann_line1_2011_2014_2015_2017_2018_2019')
  %% Petermann Line 1 2011, 2014, 2015, 2017, 2018, 2019
  pass_name = sprintf('Petermann_line1_2011_2014_2015_2017_2018_2019');
  dist_min = 500;
  master_pass_idx = 2;
  start = struct('lat',80.499370,'lon',-60.013440);
  stop = struct('lat',80.964739,'lon',-61.618307);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
  param_fn = 'rds_param_2011_Greenland_P3.xls';
  day_seg = '20110507_02';
  frms = 9:11;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  day_seg = '20140505_01';
  frms = 17:18;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2015_Greenland_C130.xls';
  day_seg = '20150505_02';
  frms = 16:17;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2017_Greenland_P3.xls';
  day_seg = '20170331_01';
  frms = 16:17;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2018_Greenland_P3.xls';
  day_seg = '20180405_01';
  frms = 15:16;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
%   param_fn = 'rds_param_2019_Greenland_P3.xls';
%   day_seg = '20190417_01';
%   frms = 15:16;
%   passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
end

if strcmpi(example_str,'Petermann_line2_2013_2014_2017')
  %% Petermann Line 2 2013, 2014, 2017
  pass_name = sprintf('Petermann_line2_2013_2014_2017');
  dist_min = 500;
  master_pass_idx = 2;
  start = struct('lat',80.517191,'lon',-59.844969);
  stop = struct('lat',80.903411,'lon',-61.332159);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
  param_fn = 'rds_param_2013_Greenland_P3.xls';
  day_seg = '20130420_02';
  frms = 4:5;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  day_seg = '20140512_01';
  frms = 16:17;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2017_Greenland_P3.xls';
  day_seg = '20170331_01';
  frms = 22:23;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post\standard'); % Folder is not in standard
end

if strcmpi(example_str,'Petermann_line4_2010_2011_2013_2014_2017')
  %% Petermann Line 4 2011, 2010, 2011, 2013, 2014, 2017
  pass_name = sprintf('Petermann_line4_2010_2011_2013_2014_2017');
  dist_min = 500;
  master_pass_idx = 4;
  start = struct('lat',80.5295278,'lon',-59.567146);
  stop = struct('lat',80.995874,'lon',-61.357055);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
  param_fn = 'rds_param_2010_Greenland_DC8.xls';
  day_seg = '20100324_01';
  frms = 25:27;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2011_Greenland_P3.xls';
  day_seg = '20110507_02';
  frms = 13:15;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2013_Greenland_P3.xls';
  day_seg = '20130420_02';
  frms = 7:8;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  day_seg = '20140512_01';
  frms = 13:14;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2017_Greenland_P3.xls';
  day_seg = '20170331_01';
  frms = 13:14;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
end

if strcmpi(example_str,'79N_line1_2010_2014_2016_2018')
  %% 79N Line 1 2010, 2014, 2016, 2018
  pass_name = sprintf('79N_line1_2010_2014_2016_2018');
  dist_min = 16000;
  master_pass_idx = 2;
  start = struct('lat',79.346109,'lon',-22.575643);
  stop = struct('lat',79.559002,'lon',-19.329911);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
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
end

if strcmpi(example_str,'Ryder_line1_2011_2013_2015_2019')
  %% Ryder line1 2011,2013,2015,2019
  pass_name = sprintf('Ryder_line1_2011_2013_2015_2019');
  dist_min = 16000;
  master_pass_idx = 3; % Use 2015 since it runs parallel to glacier
  start = struct('lat',81.537723,'lon',-50.392724);
  stop = struct('lat',81.800,'lon',-50.669);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
  param_fn = 'rds_param_2011_Greenland_P3.xls';
  day_seg = '20110509_01';
  frms = 19:20;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2013_Greenland_P3.xls';
  day_seg = '20130426_01';
  frms = 17:18;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2015_Greenland_C130.xls';
  day_seg = '20150506_02';
  frms = 19:20;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2019_Greenland_P3.xls';
  day_seg = '20190423_01';
  frms = 15:16;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
end

if strcmpi(example_str,'Steensby_line1_2011_2013_2015_2019')
  %% Steensby line1 2011,2013,2015,2019
  pass_name = sprintf('Steensby_line1_2011_2013_2015_2019');
  dist_min = 16000;
  master_pass_idx = 1;
  start = struct('lat',81.457622,'lon',-54.244113);
  stop = struct('lat',81.687467,'lon',-54.537472);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
  param_fn = 'rds_param_2011_Greenland_P3.xls';
  day_seg = '20110509_01';
  frms = 16;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
  param_fn = 'rds_param_2013_Greenland_P3.xls';
  day_seg = '20130426_01';
  frms = 15:16;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2015_Greenland_C130.xls';
  day_seg = '20150506_02';
  frms = 16:17;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2019_Greenland_P3.xls';
  day_seg = '20190423_01';
  frms = 13:14;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
end

if strcmpi(example_str,'Humboldt_line1_2012_2013_2014_2017')
  %% Humboldt line1 2012,2013,2014,2017
  pass_name = sprintf('Humboldt_line1_2012_2013_2014_2017');
  dist_min = 16000;
  master_pass_idx = 2;
  start = struct('lat',79.801630,'lon',-64.155016);
  stop = struct('lat',79.798388,'lon',-64.752165);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
  
  param_fn = 'rds_param_2012_Greenland_P3.xls';
  day_seg = '20120511_01';
  frms = 39;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2013_Greenland_P3.xls';
  day_seg = '20130420_01';
  frms = 35;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  day_seg = '20140520_03';
  frms = 31;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2017_Greenland_P3.xls';
  day_seg = '20170417_01';
  frms = 40;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
end

if strcmpi(example_str,'ZI_line1_2010_2014_2016_2017_2018_2019')
  %% Z.I. line1 2010,2014,2016,2017,2018, 2019
  pass_name = sprintf('ZI_line1_2010_2014_2016_2017_2018_2019');
  dist_min = 16000;
  master_pass_idx = 2;
  start = struct('lat',78.939138,'lon',-21.354236);
  stop = struct('lat',78.927951,'lon',-19.105975);
  input_type = 'echo';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{});
 
  param_fn = 'rds_param_2010_Greenland_P3.xls';
  day_seg = '20100525_04';
  frms = 09;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  day_seg = '20140508_01';
  frms = 46:47;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2016_Greenland_P3.xls';
  day_seg = '20160509_07';
  frms = 10;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2017_Greenland_P3.xls';
  day_seg = '20170403_02';
  frms = 4:5;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','CSARP_post/standard');
  
  param_fn = 'rds_param_2018_Greenland_P3.xls';
  day_seg = '20180418_04';
  frms = 18:19;
  passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
  
%   param_fn = 'rds_param_2019_Greenland_P3.xls';
%   day_seg = '20190405_01';
%   frms = 22:23;
%   passes(end+1) = struct('day_seg',day_seg,'frms',frms,'param_fn',param_fn,'in_path','standard');
end

if strcmpi(example_str,'camp_century_2014_2_weeks')
  %% Camp Century: 2014 Greenland P3 2 weeks
  pass_name = 'camp_century_2014_2_weeks';
  dist_min = 300;
  master_pass_idx = 8;
  start = struct('lat', 77.10,'lon', -62.3);
  stop = struct('lat', 77.13, 'lon', -61.9);
  input_type = 'sar';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{},'imgs',[]);
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140429_01','frms',67,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140515_02','frms',4,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
end

if strcmpi(example_str,'camp_century_2014_same_day')
  %% Camp Century: 2014 Greenland P3 same day
  pass_name = 'camp_century_2014_same_day';
  dist_min = 300;
  master_pass_idx = 8;
  start = struct('lat', 77.10,'lon', -62.3);
  stop = struct('lat', 77.13, 'lon', -61.9);
  input_type = 'sar';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{},'imgs',[]);
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140429_01','frms',[67 5],'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
end

if strcmpi(example_str,'camp_century_2011_2012_2013_2014')
  %% Camp Century: 2011, 2012, 2013, 2014 Greenland P3
  pass_name = sprintf('camp_century_2011_2012_2013_2014');
  dist_min = 300;
  master_pass_idx = 15+15+7+8;
  start = struct('lat', 77.10,'lon', -62.3);
  stop = struct('lat', 77.13, 'lon', -61.9);
  input_type = 'sar';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{},'imgs',[]);
  
  param_fn = 'rds_param_2011_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20110502_02','frms',32,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc]}});
  end
  param_fn = 'rds_param_2012_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20120516_01','frms',89,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc]}});
  end
  param_fn = 'rds_param_2013_Greenland_P3.xls';
  for adc = 1:7
    passes(end+1) = struct('day_seg','20130419_01','frms',4,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc]}});
  end
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140429_01','frms',67,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140429_01','frms',5,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140515_02','frms',4,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
end

if strcmpi(example_str,'summit_2012_2014_allwf')
  %% Summit Camp: 2012-2014
  pass_name = sprintf('summit_2012_2014_allwf');
  dist_min = 300;
  master_pass_idx = 8;
  start = struct('lat', 72.646,'lon', -37.898);
  stop = struct('lat', 72.791, 'lon', -38.461);
  input_type = 'sar';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{},'imgs',[]);
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140502_01','frms',41,'param_fn',param_fn, ...
      'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
  
  param_fn = 'rds_param_2012_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20120330_03','frms',8,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc]}});
  end
end

if strcmpi(example_str,'egig_2011_2012_2014_2018_allwf')
  %% EGIG line: 2011-2018
  pass_name = sprintf('egig_2011_2012_2014_2018_allwf');
  start = struct('lat', 71.147248,'lon', -36.921416);
  stop = struct('lat', 71.179219, 'lon', -36.442813);
  dist_min = 300;
  
  master_pass_idx = 8;
  input_type = 'sar';
  passes = struct('day_seg',{},'frms',{},'param_fn',{},'in_path',{},'imgs',[]);

  % Initial line:
  % 2011: 20110426_11_005: Raw data on tape
  % 2012: 20120411_02_009 and 010: Raw data on tape
  % 2014: 20140410_01_057: Raw data on tape
  % 2017: 20170506_01_057: Raw data on Indiana University tape
  % 2018: 20180501_01_052: Raw data at Indiana University DC
  %
  % Additional data to capture physical signals of interest
  % 2011: 20110426_11_006 and 007
  % 2012: 20120411_02_007 and 008
  % 2014: 20140410_01_058 and 059
  % 2017: 20170506_01_058 and 059
  % 2018: 20180501_01_051 and 053, 054, and 055
  
  param_fn = 'rds_param_2014_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20140410_01','frms',57,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  end
  
  param_fn = 'rds_param_2011_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20110426_11','frms',5,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc]}});
  end
  
  param_fn = 'rds_param_2012_Greenland_P3.xls';
  for adc = 2:16
    passes(end+1) = struct('day_seg','20120411_02','frms',[9 10],'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc]}});
  end
  
  % param_fn = 'rds_param_2017_Greenland_P3.xls';
  % for adc = 2:16
  %   passes(end+1) = struct('day_seg','20170506_01','frms',57,'param_fn',param_fn,'in_path','','imgs',{{[1 adc], [2 adc], [3 adc]}});
  % end
  
  param_fn = 'rds_param_2018_Greenland_P3.xls';
  for adc = [1:4,6:16]
    passes(end+1) = struct('day_seg','20180501_01','frms',[52],'param_fn',param_fn, ...
      'in_path','','imgs',{{[1 adc], [3 adc], [5 adc]}});
  end
  for adc = [1:4,6:16]
    passes(end+1) = struct('day_seg','20180501_01','frms',[52],'param_fn',param_fn, ...
      'in_path','','imgs',{{[2 adc], [4 adc], [6 adc]}});
  end

end

%% Automated section
% =========================================================================
%Load param variable
param.combine_passes.input_type = input_type;
param.combine_passes.pass_name = pass_name;
param.combine_passes.passes = passes;
param.combine_passes.master_pass_idx = master_pass_idx;
param.combine_passes.start = start;
param.combine_passes.stop = stop;
param.combine_passes.dist_min = dist_min;
param.day_seg = passes(master_pass_idx).day_seg;

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
