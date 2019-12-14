clearvars -except gRadar
close all
%File names for the different years
fn2011_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/rds_thule_2011_2014_wf2_insar.mat';
fn2011_2012 = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_2011_2012_wf2_insar3.mat';
fn2012_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/rds_thule_2012_2014_wf2_insar.mat';
fn2012_2013 = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_2012_2013_wf2_insar3.mat';
fn2013_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/rds_thule_2013_2014_wf2_insar.mat';
fn2014_2014 = '/cresis/snfs1/dataproducts/ct_data/rds/2014_Greenland_P3/CSARP_insar/rds_thule_2014_1Day_wf2_insar3.mat';

fns = {fn2011_2014, fn2011_2012, fn2012_2014, fn2012_2013, fn2013_2014,fn2014_2014};
leg = {'2011 to 2014', '2011 to 2012', '2012 to 2014', '2012 to 2013', '2013 to 2014','2014 to 2014'};
saven = {'PhaseSlope_Data_1Chunk'};
%Settings 
set1 = {1:15, 1:15, 1:15, 1:7, 1:15,1:15};
set2 = {16:30, 16:30, 16:30, 8:22, 16:22, 16:30};
if 0 %7 Element Setting
  set1 = {5:11, 5:11, 5:11, 5:11, 5:11};
  set2 = {20:26, 20:26, 16:22, 20:26, 20:26};
  saven  = cellfun(@(c)[c '_7Elements'],saven,'UniformOutput',false);
end
rbins = [240:420];
filt_x = 3; filt_y = 31;
filt_med = 5;
filt_psd = 301; %Also called snapshots
PSDthresh = 500;
for fn_id = 1:length(fns)
  if ~exist('tmp','var') || length(tmp)<fn_id
    tmp{fn_id} = load(fns{fn_id});
  end
  %Get the sampling frequency
  fs = tmp{fn_id}.ref.wfs(1).fs;
  %Calculate difference in days
  gpsdiff = (tmp{fn_id}.pass(set1{fn_id}(1)).gps_time(1)-tmp{fn_id}.pass(set2{fn_id}(1)).gps_time(1));
  daydiff(fn_id) = gpsdiff/(60*60*24); %seconds to days

  insar1 = mean(tmp{fn_id}.data(rbins,:,set1{fn_id}),3);
  insar2 = mean(tmp{fn_id}.data(rbins,:,set2{fn_id}),3);
  insar_data = insar2 .* conj(insar1) ./ (abs(insar1).*abs(insar2));

  %Block averaging (smooths image)
%   insar_data_filtrow = fir_dec(transpose(insar_data),ones(1,filt_y)/filt_y,1);
%   insar_data_filt = fir_dec(transpose(insar_data_filtrow),ones(1,filt_x)/filt_x);
  insar_data_filt = fir_dec(fir_dec(insar_data,ones(1,filt_y)/filt_y,1).',ones(1,filt_x)/filt_x).';
  %Normalization (makes phase begin at 0)
  if ~strcmp(leg{fn_id},'2014')
    insar_data_filt = bsxfun(@times,insar_data_filt, ...
      exp(-1i*angle(insar_data_filt(round(size(insar_data_filt,1)/2),:))));
  end  
  %Do the fft in chunks along the short time axis
  numchunks = 1;
  numsamps = size(insar_data,1);
  sampchunk = floor(numsamps/(numchunks-floor(numchunks/2)));
  %Generate number of sample points from fft
  fftpts = sampchunk*100;
  for id_ch = 1:numchunks
    startind = (sampchunk*(id_ch-1)+1)-sampchunk/2*(id_ch-1);
    endind = (sampchunk*(id_ch))-sampchunk/2*(id_ch-1);
    insar_chunk = insar_data(startind:endind,:);
    %Spectral estimation
    PSDinsar_raw = fft(insar_chunk,fftpts,1);
    PSDinsar = abs(PSDinsar_raw).^2;
    %Do median filtering
    PSDinsar_medfilt = medfilt1(PSDinsar,filt_med);
    %Incoherent averaging
    PSDinsar_avg = fir_dec(PSDinsar_medfilt,ones(1,filt_psd)/filt_psd);
    %Find the peak
    [PSDmax, PSDmaxind]  = max(PSDinsar_avg,[],1);
    %Do thresholding to grab only high-powered peaks
    PSDmaxind(PSDmax<PSDthresh) = nan;
    %Do median filtering of indexes
    PSDmaxind_medfilt = medfilt1(PSDmaxind,filt_med);
    
    %Check PSD abs^2
    if 0
      colchk = 3500;
      figure(100+id_ch);
      plot(1:length(PSDinsar(:,colchk)),PSDinsar(:,colchk))
      hold on
      if fn_id==length(fns)
        hold off
        legend(leg)
        grid on
        title(sprintf('|PSD(%d)|^2 Chunk %d',colchk,id_ch))
      end
    end
    
    %Determine the phase for each max
    %Modulate phase to compensate for phase wrapping
    phaseraw = 2*pi*(PSDmaxind/fftpts);
    phaseraw(phaseraw>pi) = phaseraw(phaseraw>pi)-2*pi;
    phaseout{fn_id}(id_ch,:) = phaseraw;
    
    %Check phaseout
    if 0
      figure(300+id_ch)
      hold on
      plot(1:length(PSDmaxind),phaseout{fn_id}(id_ch,:))
      if fn_id == length(fns)
        hold off
        legend(leg)
        grid on
        title(sprintf('Phase Out Chunk %d',id_ch))
      end
    end
  end
end

%Convert phase to physical length
eps_ice = 3.15; wfuse = 2;
kappa = 1/sqrt(eps_ice);
cls = physconst('LightSpeed');
centfreq = (tmp{1}.ref.wfs(wfuse).f0+tmp{1}.ref.wfs(wfuse).f1)/2;
BW = tmp{1}.ref.wfs(wfuse).f1-tmp{1}.ref.wfs(wfuse).f0;
efflambda = centfreq/(kappa*cls);
alongx={};
for pid = 1:length(fns)
  disp_off{pid} = (phaseout{pid}/(4*pi))*efflambda; %m
  %Convert geodetic to along track
  alongx{pid} = geodetic_to_along_track(tmp{pid}.ref.lat,tmp{pid}.ref.lon,tmp{pid}.ref.elev)/1000;
end

%Convert rbins to depth
  %Grab the rbin corresponding to the middle of the overlapping segment
  for id_ch = 1:numchunks
    startind = (sampchunk*(id_ch-1)+1)-sampchunk/2*(id_ch-1);
    endind = (sampchunk*(id_ch))-sampchunk/2*(id_ch-1);
    middleind(id_ch) = mean([startind endind])+rbins(1);
    %Grab the associated time value
    timemid(id_ch) = tmp{1}.ref.wfs(wfuse).time(floor(middleind(id_ch)));
  end
  %Convert time to depth
  binsize = (cls/(2*BW*sqrt(eps_ice)));
  crossy = middleind*binsize;

%% Save everything
save([saven{1} '.mat'],'disp_off','phaseout','middleind','alongx','leg',...
  'crossy', 'numchunks', 'BW', 'eps_ice', 'fns','daydiff')
% PlotPhaseSlope_1Chunk