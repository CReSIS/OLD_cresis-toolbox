
% Example to run a specific segment and frame by overriding parameter spreadsheet values
params = ct_set_params(params,['cmd.' cmd_method],0);
% -------------------------------------------------------------------------
% Eqip Line 1
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140414_02');
% params = ct_set_params(params,'cmd.frms',[18],'day_seg','20140414_02'); % 13
% Eqip Line 2
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20130406_01');
% params = ct_set_params(params,'cmd.frms',[19],'day_seg','20130406_01'); % 19
% -------------------------------------------------------------------------
% Equalization
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140325_07');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20110506_02');
% params = ct_set_params(params,'cmd.frms',[18],'day_seg','20110506_02'); % 18
% -------------------------------------------------------------------------
% Cody Petermann Track 1
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20110507_02');
% params = ct_set_params(params,'cmd.frms',[8:16],'day_seg','20110507_02'); % 9 10 11
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140505_01');
% params = ct_set_params(params,'cmd.frms',[16:19],'day_seg','20140505_01'); % 17 18
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20180405_01');
% params = ct_set_params(params,'cmd.frms',[14:17],'day_seg','20180405_01'); % 15 16
% Cody Petermann Track 2
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20130420_02');
% params = ct_set_params(params,'cmd.frms',[3:9],'day_seg','20130420_02'); % 4 5
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140512_01');
% params = ct_set_params(params,'cmd.frms',[12:18],'day_seg','20140512_01'); % 16 17
% Cody Petermann Track 3: All 2010
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20100324_01');
% params = ct_set_params(params,'cmd.frms',[24:30],'day_seg','20100324_01'); % 28 29
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20100420_02');
% params = ct_set_params(params,'cmd.frms',[6:9],'day_seg','20100420_02'); % 7 8
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20100420_03');
% params = ct_set_params(params,'cmd.frms',[3:6],'day_seg','20100420_03'); % 8 9
% Cody Petermann Track 4
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20100324_01');
% params = ct_set_params(params,'cmd.frms',[24:27],'day_seg','20100324_01'); % 25 26: See Track 3
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20110507_02');
% params = ct_set_params(params,'cmd.frms',[13:16],'day_seg','20110507_02'); % 13 14 15: See Track 1
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20130420_02');
% params = ct_set_params(params,'cmd.frms',[6:9],'day_seg','20130420_02'); % 7 8: See Track 2
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140512_01');
% params = ct_set_params(params,'cmd.frms',[12:15],'day_seg','20140512_01'); % 13 14: See Track 2
% Cody 79N
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20100525_04');
% params = ct_set_params(params,'cmd.frms',[10:14],'day_seg','20100525_04'); % 11 12 13
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140429_01');
% params = ct_set_params(params,'cmd.frms',[42:45],'day_seg','20140429_01'); % 43 44
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20160509_10');
% params = ct_set_params(params,'cmd.frms',[1],'day_seg','20160509_10'); % 1
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20180418_05');
% params = ct_set_params(params,'cmd.frms',[1:3],'day_seg','20180418_05'); % 1 2
% -------------------------------------------------------------------------
% CAA
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140325_05');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140325_06');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140325_07');
% params = ct_set_params(params,'cmd.frms',[4 5],'day_seg','20140325_07');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140401_03');
% params = ct_set_params(params,'cmd.frms',[1 2],'day_seg','20140401_03');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140506_01');
% params = ct_set_params(params,'cmd.frms',[3:4],'day_seg','20140506_01');
% -------------------------------------------------------------------------
% Multipass Camp Century
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140429_01');
% params = ct_set_params(params,'cmd.frms',[67]);

% -------------------------------------------------------------------------
% 2018 Antarctica Ground
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181219_01');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181217');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181219');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181220');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181221');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181222');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181223');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181224');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181225');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181226');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181227');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181228');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181229');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181231');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190102');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190103');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190104');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190105');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190106');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190107');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190108');

% params = ct_set_params(params,'radar.wfs(1).deconv.en',true);
% params = ct_set_params(params,'radar.wfs(2).deconv.en',true);

for param_idx = 1:length(params)
  param = params(param_idx);
  if strcmpi(cmd_method,'generic')
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
  elseif ~param.cmd.(cmd_method)
    continue;
  end
  for wf = 1:length(params(param_idx).radar.wfs)
    params(param_idx).radar.wfs(wf).deconv.en = 0;
    if strcmpi(params(param_idx).season_name,'2010_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2010_Greenland_DC8')
    elseif strcmpi(params(param_idx).season_name,'2011_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2012_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2013_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2014_Greenland_P3')
      if strcmp(param.day_seg,'20140429_01')
      else
        params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
        params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold';
      end
    elseif strcmpi(params(param_idx).season_name,'2016_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
    else
      params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
      params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold';
    end
  end
  for wf = 1:length(params(param_idx).radar.wfs)
    if strcmpi(params(param_idx).season_name,'2010_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2010_Greenland_DC8')
    elseif strcmpi(params(param_idx).season_name,'2011_Greenland_P3')
      params(param_idx).radar.wfs(wf).Tadc_adjust = -7.8534e-07; % 854.656e-9 shift
      %       rx_paths = params(param_idx).radar.wfs(wf).rx_paths;
      %       params(param_idx).radar.wfs(wf).chan_equal_deg = zeros(1,max(rx_paths));
      %       params(param_idx).radar.wfs(wf).chan_equal_dB = zeros(1,max(rx_paths));
      %       params(param_idx).radar.wfs(wf).Tsys = zeros(1,max(rx_paths));
      params(param_idx).radar.wfs(wf).chan_equal_deg = [-33.1 53.9 0 -175.8 116.8 84.5 96.3 -107.9 173.4 63.2 -28.2 -15.2 52.3 -178.4 73];
      params(param_idx).radar.wfs(wf).chan_equal_dB = [2.4 0.8 0 -3.6 0.1 0.5 -0.4 3.2 0 1.4 1.9 3.6 0.9 -2.5 2.4];
      params(param_idx).radar.wfs(wf).Tsys = [0 1.79 0.4 2.97 1.94 1.95 1.9 -37.59 -38.09 -34.75 -36.44 -35.97 -34.49 -37.66 -34.95]/1e9;
      params = ct_set_params(params,'qlook.out_path','qlook_equal');
    elseif strcmpi(params(param_idx).season_name,'2012_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2013_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2014_Greenland_P3')
      params(param_idx).radar.wfs(wf).Tadc_adjust = -0.00000164;
      params(param_idx).radar.wfs(wf).chan_equal_deg = [-96.4 142.3 0 -146.9 173.2 41.2 148.2 -168 -136.2 174.3 -169.6 -48.2 -131.4 -119.2 -119.2];
      params(param_idx).radar.wfs(wf).chan_equal_dB = [2.9 1.3 0 -2.5 0.5 0.1 1.4 2.5 1.5 -0.3 2 2.8 1.8 1.1 2.2];
      params(param_idx).radar.wfs(wf).Tsys = [3.73 3.09 0 5.59 3.05 0.92 2.28 -19.24 -22.03 -26.23 -29.8 -28.74 -25.71 -22.8 -19.42]/1e9;
      params = ct_set_params(params,'qlook.out_path','qlook_equal');
    elseif strcmpi(params(param_idx).season_name,'2016_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
      params(param_idx).radar.wfs(wf).Tadc_adjust = -0.00000164;
      params(param_idx).radar.wfs(wf).Tsys = [0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9;
      params(param_idx).radar.wfs(wf).chan_equal_dB = [6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2];
      params(param_idx).radar.wfs(wf).chan_equal_deg = [-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6];
      params = ct_set_params(params,'qlook.out_path','qlook_equal');
    elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
    end
  end
end

if 1
  % Standard
  params = ct_set_params(params,'array.tomo_en',false);
  params = ct_set_params(params,'array.out_path','standard');
  if strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
    params = ct_set_params(params,'array.imgs',{[ones(1,7); 6:12].', [2*ones(1,7); 6:12].', [3*ones(1,7); 6:12].'});
  elseif strcmpi(params(param_idx).season_name,'2010_Greenland_DC8')
  elseif strcmpi(params(param_idx).season_name,'2010_Greenland_P3')
  elseif strcmpi(params(param_idx).season_name,'2013_Greenland_P3')
  elseif strcmpi(params(param_idx).season_name,'2016_Greenland_P3')
  elseif any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2014_Greenland_P3'}))
    params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
  elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
  else
    keyboard
  end
  params = ct_set_params(params,'array.Nsv',1);
  params = ct_set_params(params,'array.bin_rng',[0]);
  params = ct_set_params(params,'array.line_rng',[-5:5]);
  params = ct_set_params(params,'array.Nsrc',1);
elseif 1
  % MVDR
  params = ct_set_params(params,'array.tomo_en',false);
  params = ct_set_params(params,'array.method','mvdr');
  if 1
    % 15 Elements
    params = ct_set_params(params,'array.out_path','mvdr');
    if strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
      params = ct_set_params(params,'array.imgs',{[ones(1,15); [1:4,6:16]].', [2*ones(1,15); [1:4,6:16]].', [3*ones(1,15); [1:4,6:16]].'});
    elseif any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2014_Greenland_P3'}))
      params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].', [3*ones(1,15); 2:16].'});
    else
      keyboard
    end
  else
    % 7 Elements
    params = ct_set_params(params,'array.out_path','mvdr');
    if any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2014_Greenland_P3'}))
      params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
    elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
    else
      keyboard
    end
  end
  params = ct_set_params(params,'array.Nsv',1);
  params = ct_set_params(params,'array.DCM.bin_rng',[-3:3]);
  params = ct_set_params(params,'array.DCM.line_rng',[-30:30]);
  params = ct_set_params(params,'array.bin_rng',[0]);
  params = ct_set_params(params,'array.line_rng',[-5:5]);
  params = ct_set_params(params,'array.Nsrc',1);
elseif 0
  % MUSIC
  params = ct_set_params(params,'array.tomo_en',true);
  params = ct_set_params(params,'array.out_path','music3D');
  params = ct_set_params(params,'array.method','music');
  if strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
    params = ct_set_params(params,'array.imgs',{[ones(1,15); [1:4,6:16]].', [2*ones(1,15); [1:4,6:16]].', [3*ones(1,15); [1:4,6:16]].'});
  elseif any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2014_Greenland_P3'}))
    params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].', [3*ones(1,15); 2:16].'});
  elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
  else
    keyboard
  end
  params = ct_set_params(params,'array.Nsv',128);
  params = ct_set_params(params,'array.bin_rng',[-1:1]);
  params = ct_set_params(params,'array.line_rng',[-10:10]);
  params = ct_set_params(params,'array.Nsrc',3);
elseif 0
  % GEONULL
  params = ct_set_params(params,'array.tomo_en',true);
  params = ct_set_params(params,'array.method','geonull');
  params = ct_set_params(params,'array.surf_layer.source','surfData');
  params = ct_set_params(params,'array.surf_layer.name','top twtt');
  if any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2014_Greenland_P3'}))
    params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
  elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
  else
    keyboard
  end
  params = ct_set_params(params,'array.Nsv',1);
  params = ct_set_params(params,'array.bin_rng',[0]);
  params = ct_set_params(params,'array.line_rng',[-5:5]);
  params = ct_set_params(params,'array.Nsrc',1);
elseif 0
  % GSLC
  params = ct_set_params(params,'array.tomo_en',true);
  params = ct_set_params(params,'array.out_path','gslc');
  params = ct_set_params(params,'array.method','music');
  if any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2014_Greenland_P3'}))
    params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
  elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
  else
    keyboard
  end
  params = ct_set_params(params,'array.Nsv',1);
  params = ct_set_params(params,'array.DCM.bin_rng',[-3:3]);
  params = ct_set_params(params,'array.DCM.line_rng',[-30:30]);
  params = ct_set_params(params,'array.bin_rng',[0]);
  params = ct_set_params(params,'array.line_rng',[-5:5]);
  params = ct_set_params(params,'array.Nsrc',1);
end
params = ct_set_params(params,'array.dline',6);

