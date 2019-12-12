clearvars -except gRadar tmp
close all
%Settings
set1 = {1:15, 1:15, 1:15};
set2 = {16:30, 16:30, 16:22};
rbins = [240:420];
numchunks = 3;
%File names for the different years
fn2011 = '/cresis/snfs1/dataproducts/ct_data/rds/2011_Greenland_P3/CSARP_insar/rds_thule_2011_2014_insar.mat';
fn2012 = '/cresis/snfs1/dataproducts/ct_data/rds/2012_Greenland_P3/CSARP_insar/rds_thule_2012_2014_wf1_insar.mat';
fn2013 = '/cresis/snfs1/dataproducts/ct_data/rds/2013_Greenland_P3/CSARP_insar/rds_thule_2013_2014_wf2_insar.mat';
fns = {fn2011, fn2012, fn2013};
leg = {'2011', '2012', '2013'};

for fn_id = 1%:length(fns)
  %% Load the data
  if ~exist('tmp','var') || length(tmp)<fn_id
    tmp{fn_id} = load(fns{fn_id});
  end
  %Grab the sets of interest
  insar1 = mean(tmp{fn_id}.data(rbins,:,set1{fn_id}),3);
  insar2 = mean(tmp{fn_id}.data(rbins,:,set2{fn_id}),3);
  insar_data = insar2 .* conj(insar1) ./ (abs(insar1).*abs(insar2));
  
  %% Plot the initial plots
  insar_data_filt = fir_dec(fir_dec(insar_data,ones(1,31)/31,1).',ones(1,3)/3).';
  %Generate the title string
  title_str = sprintf('%s to 2014',leg{fn_id});
  %Plot the interferogram
  fignum = fn_id + 30;
  figure(fignum); clf;
  subplot(2,1,1)
  imagesc(hsv_plot(insar_data_filt,-10));
  colormap(hsv(256))
  h_colorbar = colorbar;
  caxis([-pi pi])
  set(get(h_colorbar,'ylabel'),'string','Angle (radians)');
  grid on;
  xlabel('Range line');
  ylabel('Range bin');
  title(title_str);
  %Plot the coherence
  subplot(2,1,2)
  imagesc(abs(insar_data_filt));
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Coherence');
  grid on;
  xlabel('Range line');
  ylabel('Range bin');
  title(title_str);
  
  %% Do the cross correlation of the two data sets
  %Construct snap shot matrix 
  numchunks = 200;
  for id_ch = 1:numchunks
    if 0 %overlapping chunks
      snapx = floor(size(insar1,2)/(numchunks-floor(numchunks/2)));
      sid = (snapx*(id_ch-1)+1)-snapx/2*(id_ch-1);
      eid = (snapx*(id_ch))-snapx/2*(id_ch-1);
    elseif 1 %Non-overlapping chunks
      snapx = floor(size(insar1,2)/numchunks);
      sid = (snapx*(id_ch-1)+1);
      eid = (snapx*(id_ch));
    end
    snap1 = insar1(:,sid:eid); snap2 = insar2(:,sid:eid);
    crosssnap = snap1*snap2';
    %Do the fft of the cross correlation matrix
    R_f_shift=abs(fftshift(fft2(crosssnap)));
    %Find the maximum correlation for each chunk
    [maxval(id_ch), mid]=max(R_f_shift(:));
    [yoff(id_ch),xoff(id_ch)]=ind2sub(size(R_f_shift),mid);
%     [maxval(id_ch), mid]=max(crosssnap(:));
%     [yoff(id_ch),xoff(id_ch)]=ind2sub(size(crosssnap),mid);
  end
  %Do a weighted average based on maximum value
  avgy = mean(yoff.*abs(maxval))/mean(abs(maxval));
  avgx = mean(xoff.*abs(maxval))/mean(abs(maxval));
  offsetinx = round((avgy-avgx))
%   offsetinx = round((avgy-avgx)*snapx/2)
  %Adjust insar1 and insar 2
  insar1adj = insar1(:,offsetinx+1:end); insar2adj=insar2(:,1:end-offsetinx);
  %Plot result
  insar_dataadj = insar2adj .* conj(insar1adj) ./ (abs(insar1adj).*abs(insar2adj));
  insar_dataadj_filt = fir_dec(fir_dec(insar_dataadj,ones(1,31)/31,1).',ones(1,3)/3).';
  %Generate the title string
  title_str = sprintf('%s to 2014 (Coregistered)',leg{fn_id});
  %Plot the interferogram
  fignum = fn_id + 40;
  figure(fignum); clf;
  subplot(2,1,1)
  imagesc(hsv_plot(insar_dataadj_filt,-10));
  colormap(hsv(256))
  h_colorbar = colorbar;
  caxis([-pi pi])
  set(get(h_colorbar,'ylabel'),'string','Angle (radians)');
  grid on;
  xlabel('Range line');
  ylabel('Range bin');
  title(title_str);
  %Plot the coherence
  subplot(2,1,2)
  imagesc(abs(insar_dataadj_filt));
  h_colorbar = colorbar;
  set(get(h_colorbar,'ylabel'),'string','Coherence');
  grid on;
  xlabel('Range line');
  ylabel('Range bin');
  title(title_str);
  %% Redo to confirm result
  %Construct snap shot matrix
  clearvars snap1 snap2 crosssnap
  numchunks = 200;
  for id_ch = 1:numchunks
    if 0 %overlapping chunks
      snapx = floor(size(insar1adj,2)/(numchunks-floor(numchunks/2)));
      sid = (snapx*(id_ch-1)+1)-snapx/2*(id_ch-1);
      eid = (snapx*(id_ch))-snapx/2*(id_ch-1);
    elseif 1 %Non-overlapping chunks
      snapx = floor(size(insar1adj,2)/numchunks);
      sid = (snapx*(id_ch-1)+1);
      eid = (snapx*(id_ch));
    end
    snap1 = insar1adj(:,sid:eid); snap2 = insar2adj(:,sid:eid);
    crosssnap = snap1*snap2';
    %Do the fft of the cross correlation matrix
    R_f_shift=abs(fftshift(fft2(crosssnap)));
    %Find the maximum correlation for each chunk
    [maxval(id_ch), mid]=max(R_f_shift(:));
    [yoffadj(id_ch),xoffadj(id_ch)]=ind2sub(size(R_f_shift),mid);
%     [maxval(id_ch), mid]=max(crosssnap(:));
%     [yoffadj(id_ch),xoffadj(id_ch)]=ind2sub(size(crosssnap),mid);
  end
  %Do a weighted average based on maximum value
  avgyadj = mean(yoffadj.*abs(maxval))/mean(abs(maxval));
  avgxadj = mean(xoffadj.*abs(maxval))/mean(abs(maxval));
  offsetinxadj = round((avgyadj-avgxadj))
%   offsetinxadj = round((avgyadj-avgxadj)*snapx/2)
end



%% Use cross correlation to determine offset in x and y indices

%% Plot resulting interferogram
