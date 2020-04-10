% rds_settings
%
% Script that should be called from run scripts to do standard processing of RDS data.
% Has many seasons and different projects setup.

%% cmd
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
% Multipass Summit
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20120330_03');
% params = ct_set_params(params,'cmd.frms',[8 9]);
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20140502_01');
% params = ct_set_params(params,'cmd.frms',[41 42]);

% -------------------------------------------------------------------------
% 2018 Antarctica Ground
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181224_03');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181224_02');
params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181224_03');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181224_04');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181217');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181219');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181220');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181221');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181222');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181223');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181224');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181225');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181226');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181227');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181228');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181229');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20181231');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190102');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190103');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190104');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190105');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190106');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190107');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190108');
% params = ct_set_params(params,['cmd.' cmd_method],0,'cmd.notes','do not process');
% -------------------------------------------------------------------------
% 2019 Antarctica Ground
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20190925_04');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20200107_01');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20191230');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20191231');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20200101');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20200104');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20200105');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20200106');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20200107');
% params = ct_set_params(params,['cmd.' cmd_method],1,'day_seg','20200108');


for param_idx = 1:length(params)
  param = params(param_idx);
  if strcmpi(cmd_method,'generic')
    if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
      continue;
    end
  elseif ~param.cmd.(cmd_method)
    continue;
  end
  
  %% qlook
  params = ct_set_params(params,'qlook.out_path','qlook');
  if strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
    %params(param_idx).qlook.out_path = 'qlook_test';
    params(param_idx).qlook.out_path = 'qlook_test_adcs5678';
    %params(param_idx).qlook.out_path = 'qlook';
    params(param_idx).qlook.motion_comp = false;
    if 0
      adcs = [1:4]; Nchan = length(adcs);
      params(param_idx).qlook.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].'};
    else
      adcs = [5:8]; Nchan = length(adcs);
      params(param_idx).qlook.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].'};
    end
  elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
    %params(param_idx).qlook.out_path = 'qlook_test';
    params(param_idx).qlook.out_path = 'qlook';
    records = records_load(params(param_idx),'gps_source');
    if ~isempty(regexpi(records.gps_source,'cresis'))
      params(param_idx).qlook.motion_comp = true;
    else
      params(param_idx).qlook.motion_comp = false;
    end
    params(param_idx).qlook.motion_comp = false;
    adcs = [1:6]; Nchan = length(adcs);
    if length(params(param_idx).radar.wfs) == 3
      params(param_idx).qlook.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan); adcs].'};
    else
      params(param_idx).qlook.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan) 4*ones(1,Nchan); adcs adcs].'};
    end
  end
  
  %% sar
  if strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
    params = ct_set_params(params,'sar.out_path','sar');
    params = ct_set_params(params,'sar.sigma_x',1);
    params(param_idx).sar.mocomp.en = true;
  elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
    params = ct_set_params(params,'sar.out_path','sar');
    %params = ct_set_params(params,'sar.out_path','sar_tukey');
    params = ct_set_params(params,'sar.sigma_x',1);
    records = records_load(params(param_idx),'gps_source');
    if ~isempty(regexpi(records.gps_source,'cresis'))
      params(param_idx).sar.mocomp.en = true;
    else
      params(param_idx).sar.mocomp.en = false;
    end
    params(param_idx).sar.mocomp.en = false;
  end
  
  %% general
  if strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
  elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
    %     params = ct_set_params(params,'records.data_map',{[2 0 1 1;2 1 1 2;5 0 2 1;5 1 2 2;8 0 3 1;8 1 3 2],[2 0 1 3;2 1 1 4;5 0 2 3;5 1 2 4;8 0 3 3;8 1 3 4],[2 0 1 5;2 1 1 6;5 0 2 5;5 1 2 6;8 0 3 5;8 1 3 6],[2 0 1 7;2 1 1 8;5 0 2 7;5 1 2 8;8 0 3 7;8 1 3 8]});
    %     params = ct_set_params(params,'records.data_map',{[2 0 1 1; 2 1 1 2; 5 0 2 1; 5 1 2 2; 8 0 3 1; 8 1 3 2; 11 0 4 1; 11 1 4 2],[2 0 ...
    %       1 3; 2 1 1 4; 5 0 2 3; 5 1 2 4; 8 0 3 3; 8 1 3 4; 11 0 4 3; 11 1 4 4],[2 0 1 5; 2 1 1 6; 5 0 2 5; 5 1 2 6; 8 0 3 ...
    %       5; 8 1 3 6; 11 0 4 5; 11 1 4 6],[2 0 1 7; 2 1 1 8; 5 0 2 7; 5 1 2 8; 8 0 3 7; 8 1 3 8; 11 0 4 7; 11 1 4 8]}, ...
    %       'day_seg','20200101_02|20200101_03|20200101_04|20200104_01|20200104_02|20200104_03|20200104_04|20200104_05|20200105_01|20200105_02|20200105_03|20200105_04|20200106_01|20200106_02|20200107_01|20200107_02|20200107_03|20200108_01|20200108_02');
  end
  
  %% radar.wfs
  for wf = 1:length(params(param_idx).radar.wfs)
    params(param_idx).radar.wfs(wf).deconv.en = 0;
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
    elseif strcmpi(params(param_idx).season_name,'2012_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2013_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2014_Greenland_P3')
      if any(strcmp(param.day_seg,{'20140429_01','20140502_01'}))
      else
        params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
        params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold';
      end
      params(param_idx).radar.wfs(wf).Tadc_adjust = -0.00000164;
      params(param_idx).radar.wfs(wf).chan_equal_deg = [-96.4 142.3 0 -146.9 173.2 41.2 148.2 -168 -136.2 174.3 -169.6 -48.2 -131.4 -119.2 -119.2];
      params(param_idx).radar.wfs(wf).chan_equal_dB = [2.9 1.3 0 -2.5 0.5 0.1 1.4 2.5 1.5 -0.3 2 2.8 1.8 1.1 2.2];
      params(param_idx).radar.wfs(wf).Tsys = [3.73 3.09 0 5.59 3.05 0.92 2.28 -19.24 -22.03 -26.23 -29.8 -28.74 -25.71 -22.8 -19.42]/1e9;
    elseif strcmpi(params(param_idx).season_name,'2016_Greenland_P3')
    elseif strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
      params(param_idx).radar.wfs(wf).Tadc_adjust = -0.00000164;
      params(param_idx).radar.wfs(wf).Tsys = [0.46 -4.66 0.14 -1.77 0 -2.63 -3.38 -69.66 -75.57 -75.45 -80.42 -80.49 -75.71 -77.69 -70.53]/1e9;
      params(param_idx).radar.wfs(wf).chan_equal_dB = [6.8 -0.6 3 0.1 0 3.5 3.9 7 3.3 4.8 6.1 6.2 4.6 3.1 6.2];
      params(param_idx).radar.wfs(wf).chan_equal_deg = [-166.2 -142.7 177 -95.9 0 -25.9 -86.5 -27.4 128.1 41.6 -46.8 43 90.7 121.3 31.6];
    elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
      params(param_idx).radar.wfs(wf).Tadc_adjust = -0.259e-6;
      if wf == 1
        params(param_idx).radar.wfs(wf).chan_equal_dB = [-0.9 0.7 0 0.5 1.6 -0.7 -0.8 0.4];
        params(param_idx).radar.wfs(wf).chan_equal_deg = [-172.7 -173.5 -174.8 -169.3 -51.1 0 -26 -34.5];
      elseif wf == 2
        params(param_idx).radar.wfs(wf).chan_equal_dB = [-0.9 0.7 0 0.5 1.6 -0.7 -0.8 0.4];
        params(param_idx).radar.wfs(wf).chan_equal_deg = [-172.7 -173.5 -174.8 -169.3 -51.1 0 -26 -34.5];
      end
    elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
      if strcmpi(cmd_method,'sar')
        if strcmp(params(param_idx).sar.out_path,'sar')
          params(param_idx).radar.wfs(wf).ft_wind = @hanning;
          params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
          params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold';
        elseif strcmp(params(param_idx).sar.out_path,'sar_tukey')
          params(param_idx).radar.wfs(wf).ft_wind = @(N) tukeywin_trim(N,0.5);
          params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
          params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold_tukey';
        end
      elseif strcmpi(cmd_method,'qlook')
        if strcmp(params(param_idx).qlook.out_path,'qlook')
          params(param_idx).radar.wfs(wf).ft_wind = @(N) tukeywin_trim(N,0.5);
          params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
          params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold_tukey';
        elseif strcmp(params(param_idx).qlook.out_path,'qlook_test')
          params(param_idx).radar.wfs(wf).ft_wind = @(N) tukeywin_trim(N,0.5);
        end
      end
      params(param_idx).radar.wfs(wf).Tadc_adjust = -0.25e-6;
      if wf  < 3
        params(param_idx).radar.wfs(wf).chan_equal_dB = [10 6 10 9.6 -11.4 -0.8 0 9];
        params(param_idx).radar.wfs(wf).chan_equal_deg = [-53.4 -60.5 116.8 134 118.4 64.7 0 -57];
      elseif wf == 3
        params(param_idx).radar.wfs(wf).chan_equal_dB = [10.2 6.4 17.8 12.6 15.3 11.7 0 9.2];
        params(param_idx).radar.wfs(wf).chan_equal_deg = [-52.7 -59 117.5 134.2 92.6 73.2 0 -57.6];
      elseif wf == 4
        params(param_idx).radar.wfs(wf).chan_equal_dB = [10.6 6.8 20.2 16.3 20.2 6.1 1.4 10.6];
        params(param_idx).radar.wfs(wf).chan_equal_deg = [-47.1 -53.7 127.2 147.1 111.9 122.6 11.2 -35.2];
      end
    else
      params(param_idx).radar.wfs(wf).coh_noise_method = 'analysis';
      params(param_idx).radar.wfs(wf).coh_noise_arg.fn = 'analysis_threshold';
    end
  end
  
  %% analysis coh_noise cmd
  if isfield(params(param_idx),'analysis') && strcmp(params(param_idx).analysis.cmd{1}.method,'coh_noise')
    for wf = 1:length(params(param_idx).radar.wfs)
      params(param_idx).radar.wfs(wf).bad_value = NaN;
    end
    
    % radar.wfs
    if isfield(param_override,'analysis') && isfield(param_override.analysis,'out_path')
      for wf = 1:length(params(param_idx).radar.wfs)
        if ~isempty(regexp(param_override.analysis.out_path,'tukey'))
          params(param_idx).radar.wfs(wf).ft_wind = @(N) tukeywin_trim(N,0.5);
          if ~isempty(regexp(param_override.analysis.out_path,'threshold'))
            params = ct_set_params(params,'analysis.cmd{1}.threshold','analysis_tukey');
          end
        else
          if ~isempty(regexp(param_override.analysis.out_path,'threshold'))
            params = ct_set_params(params,'analysis.cmd{1}.threshold','');
          end
        end
      end
    end
    
    if isfield(param_override,'collate_coh_noise')
      if strcmp(param_override.collate_coh_noise.in_path,'analysis_tukey')
        param_override.collate_coh_noise.debug_out_dir = 'collate_coh_noise_tukey';
      elseif strcmp(param_override.collate_coh_noise.in_path,'analysis')
        param_override.collate_coh_noise.debug_out_dir = 'collate_coh_noise';
      elseif strcmp(param_override.collate_coh_noise.in_path,'analysis_threshold')
        param_override.collate_coh_noise.debug_out_dir = 'collate_coh_noise_threshold';
      elseif strcmp(param_override.collate_coh_noise.in_path,'analysis_threshold_tukey')
        param_override.collate_coh_noise.debug_out_dir = 'collate_coh_noise_threshold_tukey';
      end
      for img = 1:length(params(param_idx).analysis.imgs)
        for wf_adc = 1:size(params(param_idx).analysis.imgs{img},1)
          wf = params(param_idx).analysis.imgs{img}(wf_adc,1);
          adc = params(param_idx).analysis.imgs{img}(wf_adc,2);
          Tpd = params(param_idx).radar.wfs(wf).Tpd;
          if strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
            params(param_idx).collate_coh_noise.method{img} = 'firdec';
            params(param_idx).collate_coh_noise.firdec_fs{img} = 1/7.5;
            params(param_idx).collate_coh_noise.firdec_fcutoff{img} = @(t)(t<Tpd+1.2e-6)*1/120+(t>=Tpd+1.2e-6)*1/120;
            if wf == 1
              params(param_idx).collate_coh_noise.threshold_eval{img}{wf_adc} ...
                = 'threshold = max_filt1(threshold+6,11);';
            elseif wf == 2
              params(param_idx).collate_coh_noise.threshold_eval{img}{wf_adc} ...
                = 'threshold = max_filt1(threshold+6,11);';
            elseif any(wf == [3 4])
              if any(adc == [1 2 3 4 5 6 7 8])
                params(param_idx).collate_coh_noise.threshold_eval{img}{wf_adc} ...
                  = 'threshold = max(nanmin(threshold(time>Tpd+1.2e-6 & time<time(end)-Tpd))+6*ones(size(threshold)),max_filt1((10*log10(abs(dft_noise(:,1)).^2)+9).*(time<Tpd+0e-6) + (threshold+9).*(time>=Tpd+0e-6 & time<=Tpd+13e-6)-1e6*(time>(Tpd+13e-6)),21));';
              end
            else
              params(param_idx).collate_coh_noise.threshold_eval{img}{wf_adc} = 'threshold = -inf(size(threshold));';
            end
          else
            % Coherent noise estimate by finding DC of the entire segment
            % (as opposed to using a firdec low pass filter):
            params(param_idx).collate_coh_noise.method{img} = 'dft';
            params(param_idx).collate_coh_noise.dft_corr_time{img} = inf;
            
            % For time >= Tpd+1.2e-6, the minimum threshold used that is
            % estimated by the collate process is used in the region of time >=
            % Tpd+1.2e-6 and time(end)-Tpd. The first part of the time gate is
            % ignored because the transmission is occuring and the last part is
            % ignored due to pulse compression roll off.  A guard of 6 dB is
            % added to this threshold to ensure random noise fluctuation is
            % accounted for. There is some risk that this data dependent minimum
            % threshold value is too high and a fixed threshold should be used
            % instead of "nanmin(threshold(time>Tpd+1.2e-6 & time<time(end)-Tpd))"
            %
            % For time < Tpd+1.2e-6, the max_filt1 filtered version of the
            % minimum of the (DC coherent noise + 15 dB) and the (threshold + 6
            % dB) waveforms are used. max_filt1 is used because the nulls
            % sometimes shift around with time. Threshold+6 dB is ideal, but is
            % occasionally too high due to signal interference so the DC average
            % is also used.
            params(param_idx).collate_coh_noise.threshold_eval{img} = 'threshold = max(nanmin(threshold(time>Tpd+1.2e-6 & time<time(end)-Tpd))+6*ones(size(threshold)),max_filt1(min(threshold+6,10*log10(abs(noise_dft(:,1)).^2)+15)-1e6*(time>(Tpd+1.2e-6)),5));';
          end
        end
      end
    end
  end
  
  %% analysis waveform cmd (equalization)
  if isfield(params(param_idx),'analysis') && strcmp(params(param_idx).analysis.cmd{1}.method,'waveform') && isfield(param_override,'analysis')
    if strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
      if strcmp(param_override.analysis.out_path,'analysis_equal_001')
        params(param_idx).analysis.imgs = {[1*ones([8 1]),(1:8).']};
        params(param_idx).analysis.cmd{1}.start_time = struct('name','equal_001','eval',struct('cmd','s=s-0.75e-6;'));
      elseif strcmp(param_override.analysis.out_path,'analysis_equal_002')
        params(param_idx).analysis.imgs = {[2*ones([8 1]),(1:8).']};
        params(param_idx).analysis.cmd{1}.start_time = struct('name','equal_002','eval',struct('cmd','s=s-0.75e-6;'));
      end
    end
    if strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
      if strcmp(param_override.analysis.out_path,'analysis_equal_002')
        params(param_idx).analysis.imgs = {[2*ones([8 1]),(1:8).']};
        params(param_idx).analysis.cmd{1}.start_time = struct('name','equal_002','eval',struct('cmd','s=s-1.5e-6;'));
      elseif strcmp(param_override.analysis.out_path,'analysis_equal_003')
        params(param_idx).analysis.imgs = {[[3*ones([8 1]),(1:8).']; [4*ones([8 1]),(1:8).']]};
        params(param_idx).analysis.cmd{1}.start_time = struct('name','equal_003','eval',struct('cmd','s=s-1.5e-6;'));
      end
    end
  end
  
  if isfield(param_override,'collate_equal')
    if strcmpi(params(param_idx).season_name,'2014_Greenland_P3') &&  strcmpi(params.day_seg,'20140325_07')
      param_override.collate_equal.rlines = [14000:18500];
      param_override.collate_equal.ref = 3;
      param_override.collate_equal.debug_plots = {'before_comp','after_comp','surf','final','comp_image'};
      param_override.collate_equal.retrack_en = true;
      
    elseif strcmpi(params(param_idx).season_name,'2011_Greenland_P3') && strcmpi(params.day_seg,'20110506_02')
      if 1
        param_override.collate_equal.img_lists = {[2]};
        param_override.collate_equal.rlines = [170000:180000];
      else
        param_override.collate_equal.img_lists = {[1]};
        param_override.collate_equal.rlines = [80000:117500];
      end
      param_override.collate_equal.wf_adcs = {{[2:16]}};
      param_override.collate_equal.ref = 3;
      param_override.collate_equal.debug_plots = {'before_comp','after_comp','surf','final','comp_image'};
      
      param_override.collate_equal.retrack_en = true;
      
    elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground') && strcmpi(params.day_seg,'20181224_03')
      if strcmp(param_override.collate_equal.in_path,'analysis_equal_001')
        params(param_idx).analysis.imgs = {[1*ones([8 1]),(1:8).']};
        param_override.collate_equal.rlines = []; % wf == 1, equal_001 layer
        param_override.collate_equal.debug_out_dir = 'collate_equal_001';
        param_override.collate_equal.ref = 6;
      elseif strcmp(param_override.collate_equal.in_path,'analysis_equal_002')
        params(param_idx).analysis.imgs = {[2*ones([8 1]),(1:8).']};
        param_override.collate_equal.rlines = []; % wf == 2 equal_002 layer
        param_override.collate_equal.debug_out_dir = 'collate_equal_002';
        param_override.collate_equal.ref = 6;
      end
      %param_override.collate_equal.debug_plots = {'visible','before_comp','after_comp','surf','final','comp_image'};
      param_override.collate_equal.debug_plots = {'before_comp','after_comp','surf','final','comp_image'};
      param_override.collate_equal.retrack_en = false;
      
    elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground') && strcmpi(params.day_seg,'20200107_01')
      if strcmp(param_override.collate_equal.in_path,'analysis_equal_002')
        params(param_idx).analysis.imgs = {[2*ones([8 1]),(1:8).']};
        param_override.collate_equal.rlines = []; % wf == 2 equal_002 layer
        param_override.collate_equal.debug_out_dir = 'collate_equal_002';
      elseif strcmp(param_override.collate_equal.in_path,'analysis_equal_003')
        params(param_idx).analysis.imgs = {[[3*ones([8 1]),(1:8).']; [4*ones([8 1]),(1:8).']]};
        param_override.collate_equal.rlines = []; % wf == 3,4, equal_003 layer
        param_override.collate_equal.debug_out_dir = 'collate_equal_003';
      end
      param_override.collate_equal.debug_plots = {'visible','before_comp','after_comp','surf','final','comp_image'};
      %       param_override.collate_equal.debug_plots = {'before_comp','after_comp','surf','final','comp_image'};
      param_override.collate_equal.retrack_en = false;
      
    end
    
  end
  
end

%% array
if isfield(param_override,'array') && isfield(param_override.array,'out_path')
  for param_idx = 1:length(params)
    param = params(param_idx);
    if strcmpi(cmd_method,'generic')
      if ~isfield(param.cmd,'generic') || iscell(param.cmd.generic) || ischar(param.cmd.generic) || ~param.cmd.generic
        continue;
      end
    elseif ~param.cmd.(cmd_method)
      continue;
    end
    if strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
      params(param_idx).array.in_path = 'sar';
    elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
      if strcmpi(param_override.array.out_path,'standard')
        params(param_idx).array.in_path = 'sar';
      elseif strcmpi(param_override.array.out_path,'standard_tukey')
        params(param_idx).array.in_path = 'sar_tukey';
      elseif strcmpi(param_override.array.out_path,'mvdr')
        params(param_idx).array.in_path = 'sar';
      elseif strcmpi(param_override.array.out_path,'music3D')
        params(param_idx).array.in_path = 'sar';
      end
    end
    params(param_idx).array.ft_over_sample = 2;
    
    if ~isempty(regexp(param_override.array.out_path,'standard'))
      % Standard
      params(param_idx).array.tomo_en = false;
      params(param_idx).array.method = 'standard';
      if strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 6:12].', [2*ones(1,7); 6:12].', [3*ones(1,7); 6:12].'});
      elseif strcmpi(params(param_idx).season_name,'2010_Greenland_DC8')
      elseif strcmpi(params(param_idx).season_name,'2010_Greenland_P3')
      elseif strcmpi(params(param_idx).season_name,'2011_Greenland_P3')
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].'});
      elseif strcmpi(params(param_idx).season_name,'2012_Greenland_P3')
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].'});
      elseif strcmpi(params(param_idx).season_name,'2013_Greenland_P3')
      elseif strcmpi(params(param_idx).season_name,'2016_Greenland_P3')
      elseif any(strcmpi(params(param_idx).season_name,{'2014_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
      elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
        params(param_idx).array.imgs = {{[ones(1,4); 1:4].',[ones(1,4); 5:8].'},{[2*ones(1,4); 1:4].',[2*ones(1,4); 5:8].'}};
      elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
        adcs = [1:6]; Nchan = length(adcs);
        if length(params(param_idx).radar.wfs) == 3
          params(param_idx).array.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan); adcs].'};
        else
          params(param_idx).array.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan) 4*ones(1,Nchan); adcs adcs].'};
        end
      else
        keyboard
      end
      params(param_idx).array.Nsv = 1;
      params(param_idx).array.bin_rng = [0];
      params(param_idx).array.line_rng = [-5:5];
      params(param_idx).array.Nsrc = 1;
    elseif ~isempty(regexp(param_override.array.out_path,'mvdr'))
      % MVDR
      params(param_idx).array.tomo_en = false;
      params(param_idx).array.method = 'mvdr';
      if strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
        adcs = [1:6]; Nchan = length(adcs);
        if length(params(param_idx).radar.wfs) == 3
          params(param_idx).array.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan); adcs].'};
        else
          params(param_idx).array.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan) 4*ones(1,Nchan); adcs adcs].'};
        end
      else
        if 1
          % 15 Elements
          params = ct_set_params(params,'array.out_path','mvdr');
          if strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
            params = ct_set_params(params,'array.imgs',{[ones(1,15); [1:4,6:16]].', [2*ones(1,15); [1:4,6:16]].', [3*ones(1,15); [1:4,6:16]].'});
          elseif any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2012_Greenland_P3'}))
            params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].'});
          elseif any(strcmpi(params(param_idx).season_name,{'2014_Greenland_P3'}))
            params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].', [3*ones(1,15); 2:16].'});
          else
            keyboard
          end
        else
          % 7 Elements
          params = ct_set_params(params,'array.out_path','mvdr');
          if any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2012_Greenland_P3','2014_Greenland_P3'}))
            params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
          elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
          elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
          else
            keyboard
          end
        end
      end
      params(param_idx).array.Nsv = 1;
      params(param_idx).array.bin_rng = [0];
      params(param_idx).array.line_rng = [-5:5];
      params(param_idx).array.Nsrc = 1;
      params(param_idx).array.DCM.bin_rng = [-3:3];
      params(param_idx).array.DCM.line_rng = [-30:30];
    elseif ~isempty(regexp(param_override.array.out_path,'music3D'))
      % MUSIC
      params(param_idx).array.tomo_en = true;
      params(param_idx).array.method = 'music';
      if strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
        params = ct_set_params(params,'array.imgs',{[ones(1,15); [1:4,6:16]].', [2*ones(1,15); [1:4,6:16]].', [3*ones(1,15); [1:4,6:16]].'});
        params = ct_set_params(params,'array.Nsv',128);
        params = ct_set_params(params,'array.bin_rng',[-1:1]);
        params = ct_set_params(params,'array.line_rng',[-10:10]);
        params = ct_set_params(params,'array.Nsrc',3);
      elseif any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2012_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].'});
        params = ct_set_params(params,'array.Nsv',128);
        params = ct_set_params(params,'array.bin_rng',[-1:1]);
        params = ct_set_params(params,'array.line_rng',[-10:10]);
        params = ct_set_params(params,'array.Nsrc',3);
      elseif any(strcmpi(params(param_idx).season_name,{'2014_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].', [3*ones(1,15); 2:16].'});
        params = ct_set_params(params,'array.Nsv',128);
        params = ct_set_params(params,'array.bin_rng',[-1:1]);
        params = ct_set_params(params,'array.line_rng',[-10:10]);
        params = ct_set_params(params,'array.Nsrc',3);
      elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
      elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
        adcs = [1:6]; Nchan = length(adcs);
        if length(params(param_idx).radar.wfs) == 3
          params(param_idx).array.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan); adcs].'};
        else
          params(param_idx).array.imgs = {[ones(1,Nchan); adcs].', [2*ones(1,Nchan); adcs].', [3*ones(1,Nchan) 4*ones(1,Nchan); adcs adcs].'};
        end
        params(param_idx).array.Nsv = 64;
        params(param_idx).array.bin_rng = [-1:1];
        params(param_idx).array.line_rng = [-10:10];
        params(param_idx).array.Nsrc = 2;
      else
        keyboard
      end
    elseif 1
      % GEONULL
      params = ct_set_params(params,'array.tomo_en',true);
      params = ct_set_params(params,'array.method','geonull_cal');
      params = ct_set_params(params,'array.surf_layer.source','surfData');
      params = ct_set_params(params,'array.surf_layer.name','top twtt');
      if any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2012_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].'});
      elseif any(strcmpi(params(param_idx).season_name,{'2014_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
      elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
      elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
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
      if any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2012_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].'});
      elseif any(strcmpi(params(param_idx).season_name,{'2014_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
      elseif strcmpi(params(param_idx).season_name,'2018_Antarctica_Ground')
      elseif strcmpi(params(param_idx).season_name,'2019_Antarctica_Ground')
      else
        keyboard
      end
      params = ct_set_params(params,'array.Nsv',1);
      params = ct_set_params(params,'array.DCM.bin_rng',[-3:3]);
      params = ct_set_params(params,'array.DCM.line_rng',[-30:30]);
      params = ct_set_params(params,'array.bin_rng',[0]);
      params = ct_set_params(params,'array.line_rng',[-5:5]);
      params = ct_set_params(params,'array.Nsrc',1);
    elseif 1
      % SNAPSHOT
      params = ct_set_params(params,'array.tomo_en',true);
      params = ct_set_params(params,'array.in_path','sar_air');
      params = ct_set_params(params,'array.out_path','snapshot');
      params = ct_set_params(params,'array.method','snapshot');
      params = ct_set_params(params,'array.surf_layer.source','surfData');
      params = ct_set_params(params,'array.surf_layer.name','top twtt');
      if strcmpi(params(param_idx).season_name,'2018_Greenland_P3')
        params = ct_set_params(params,'array.imgs',{[ones(1,15); [1:4,6:16]].', [2*ones(1,15); [1:4,6:16]].', [3*ones(1,15); [1:4,6:16]].'});
      elseif any(strcmpi(params(param_idx).season_name,{'2011_Greenland_P3','2012_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].'});
      elseif any(strcmpi(params(param_idx).season_name,{'2014_Greenland_P3'}))
        params = ct_set_params(params,'array.imgs',{[ones(1,15); 2:16].', [2*ones(1,15); 2:16].', [3*ones(1,15); 2:16].'});
      else
        keyboard
      end
      %   params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
      params = ct_set_params(params,'array.imgs',{[ones(1,7); 2:8].', [2*ones(1,7); 2:8].', [3*ones(1,7); 2:8].'});
      params = ct_set_params(params,'array.dline',1);
    end
    params = ct_set_params(params,'array.dline',6);
  end
end
