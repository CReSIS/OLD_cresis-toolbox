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
%Set tile width and length
winrow = 150; wincol = 450;
mastrow = 75; mastcol = 225; % These two define the maximum window size
slavrow = mastrow; slavcol = mastcol;
rowoverlap = 50; coloverlap = 25;
%Set rlines
rlines = 1:900;
for fn_id = 1%:length(fns)
  %% Load the data
  if ~exist('tmp','var') || length(tmp)<fn_id
    tmp{fn_id} = load(fns{fn_id});
  end
  %Grab the two sets of data
  slaveinsar = mean(tmp{fn_id}.data(rbins,rlines,set1{fn_id}),3);
  masterinsar = mean(tmp{fn_id}.data(rbins,rlines,set2{fn_id}),3);
  if 0 %Test slave with shifted master
    rowshift = 3; colshift = 100;
    fprintf('\n Slave to Master Comparison \n\tRowOff: %d\n\tColOff: %d\n',rowshift,colshift);
    slaveinsar = zeros(size(masterinsar));
    slaveinsar(rowshift+1:end,colshift+1:end)=masterinsar(1:end-rowshift,1:end-colshift);
  end
  
  if 1
    figure(1+fn_id); clf
    subplot(2,1,1)
    imagesc(lp(masterinsar(:,1:4*mastcol)));
    hold on
    title('Master')
    subplot(2,1,2)
    imagesc(lp(slaveinsar(:,1:4*mastcol)));
    hold on
    title('Slave')
  end
    
  %Sweep through the slave tiles
  convout = {}; rowout = []; colout = []; convMARout =[];
  colidx = 0; colind_end = 0;
  while colind_end <size(masterinsar,2)
    %Determine the column indices for the master image
    colind_str = colind_end-coloverlap;
    if colind_str<1;
      colind_str = 1;
    end
    colind_end = colind_str+mastcol;
    if colind_end >=size(masterinsar,2)
      colind_end = size(masterinsar,2);
    end
    colinds_mast = colind_str:colind_end;
    %Determine the column indices for the slave image
    colind_strsl = (min(colinds_mast)+floor((mastcol-slavcol)/2));
    colind_endsl = (max(colinds_mast)-floor((mastcol-slavcol)/2));
    colinds_slav = colind_strsl:colind_endsl;
    rowidx = 0; rowind_end = 0;
    while rowind_end<size(masterinsar,1)
      %Determine the row indices for the master image
      rowind_str = rowind_end-rowoverlap;
      if rowind_str<1
        rowind_str = 1;
      end
      rowind_end = rowind_str+mastrow;
      if rowind_end >=size(masterinsar,1)
        rowind_end = size(masterinsar,1);
      end
      rowinds_mast = rowind_str:rowind_end;
      %Determine row indices for the slave image
      rowind_strsl = (min(rowinds_mast)+floor((mastrow-slavrow)/2));
      rowind_endsl = (max(rowinds_mast)-floor((mastrow-slavrow)/2));
      rowinds_slav = rowind_strsl:rowind_endsl;
      if 1
        figure(1+fn_id)
        subplot(2,1,1)
        rowplot = [rowind_str, rowind_str, rowind_end, rowind_end, rowind_str];
        colplot = [colind_str, colind_end, colind_end, colind_str, colind_str];
        plot(colplot, rowplot,'b','Linewidth',2)
        subplot(2,1,2)
        rowplot = [rowind_strsl, rowind_strsl, rowind_endsl, rowind_endsl, rowind_strsl];
        colplot = [colind_strsl, colind_endsl, colind_endsl, colind_strsl, colind_strsl];
        plot(colplot, rowplot,'r','Linewidth',2)
        drawnow
      end 
      %% Do Cross-Correlation in frequency Domain
      %Generate master tile
      mainsar = masterinsar(rowinds_mast,colinds_mast);
      madata = zeros(winrow,wincol);
%       macols = floor((wincol-mastcol)/2):floor((wincol-mastcol)/2)+size(mainsar,2)-1;
%       marows= floor((winrow-mastrow)/2):floor((winrow-mastrow)/2)+size(mainsar,1)-1;
      macols = 1:size(mainsar,2); marows = 1:size(mainsar,1);
      madata(marows,macols) = mainsar;
      %Generate slave tile
      slinsar = slaveinsar(rowinds_slav,colinds_slav);
      %Pad slave tile with zeros to match master tile size
      sldata = zeros(size(madata));
%       slcols = floor((wincol-slavcol)/2):floor((wincol-slavcol)/2)+size(slinsar,2)-1;
%       slrows = floor((winrow-slavrow)/2):floor((winrow-slavrow)/2)+size(slinsar,1)-1;
      slcols = 1:size(slinsar,2); slrows = 1:size(slinsar,1);
      sldata(slrows,slcols) = slinsar;
      %Do fft of slave tile and master tile
      slfft = []; mafft = [];
      slfft = fft2(sldata); mafft = fft2(madata);
      %Do convolution (multiplication in freq domain)
      fftmult = [];
      fftmult = mafft.*conj(slfft); 
      %Take back to time domain (oversample ~x20)
      oversample_factor = 20; insarconv = [];
      insarconv = abs(ifft2(fftmult,size(fftmult,1)*oversample_factor,size(fftmult,2)*oversample_factor));
      %Find the peak of the convolution
      [peak_row, peak_col] = ind2sub(size(insarconv),find(insarconv==max(max(insarconv)),1));
      %Using region around peak, do quadratic interpolation to determine
        %"true" peak
        regionrows = 20; regioncols = 75; interp_interval = 0.25;
        peak_strow =  peak_row-floor(regionrows/2);
        peak_endrow = peak_row+floor(regionrows/2);
        if peak_endrow>size(insarconv,1)
          peak_endrow = size(insarconv,1);
        elseif peak_strow<=1
          peak_strow = 1;
        end
        prowinds = peak_strow:peak_endrow;
        prowinterp = peak_strow:interp_interval:peak_endrow;
        
        peak_stcol = peak_col-floor(regioncols/2);
        peak_endcol = peak_col+floor(regioncols/2);
        if peak_endcol>size(insarconv,2)
          peak_endcol = size(insarconv,2);
        elseif peak_stcol<=1
          peak_stcol = 1;
        end
        pcolinds = peak_stcol:peak_endcol;
        pcolinterp = peak_stcol:interp_interval:peak_endcol;
      peak_region = insarconv(prowinds,pcolinds);
      [colinterp, rowinterp] = meshgrid(pcolinds,prowinds);
      [colq, rowq] = meshgrid(pcolinterp,prowinterp);
      peak_region_interp = interp2(colinterp, rowinterp, peak_region, colq, rowq,'spline');
      [peak_interp_row, peak_interp_col] = ind2sub(size(peak_region_interp),find(peak_region_interp==max(max(peak_region_interp)),1));
      col_offset = colq(peak_interp_row, peak_interp_col);
      row_offset = rowq(peak_interp_row, peak_interp_col);
%       col_offset = peak_col; row_offset = peak_row;
      %% Plot the results
      if 0
        %Check the reverse fft result
        sl_ifft = ifft2(slfft,size(slfft,1)*oversample_factor,size(slfft,2)*oversample_factor);
        ma_ifft = ifft2(mafft,size(mafft,1)*oversample_factor,size(mafft,2)*oversample_factor);
        %Plot it to check
        figure(100+fn_id); clf
        subplot(2,1,1)
        imagesc(lp(sl_ifft));
        colorbar
%         caxis([-150 -100])
        title('ifft');
        subplot(2,1,2)
        imagesc(lp(sldata))
        colorbar
        caxis([-150 -100])
        title('Original')
        
        figure(200+fn_id); clf
        subplot(2,1,1)
        imagesc(lp(ma_ifft));
        colorbar
        caxis([-200 -150])
        title('ifft Master');
        subplot(2,1,2)
        imagesc(lp(madata))
        colorbar
        caxis([-150 -100])
        hold on
        slcoloff = floor((wincol-slavcol)/2); slrowoff = floor((winrow-slavrow)/2);
        rowplot = [slrowoff, slrowoff, size(madata,1)-slrowoff, size(madata,1)-slrowoff, slrowoff];
        colplot = [slcoloff, size(madata,2)-slcoloff, size(madata,2)-slcoloff, slcoloff, slcoloff];
        plot( colplot,rowplot,'r','LineWidth',2)
        hold off
        title('Original Master')
      end
      if 0
        figure(10+fn_id); clf
        imagesc(lp(insarconv'));
        colorbar
        hold on
        plot(row_offset,col_offset,'rd','MarkerFaceColor','r','Markersize',10);
        hold off
        title({'Convolution',sprintf('Max @ Row: %.3f of %.0f and Column: %.3f of %.0f',row_offset,max(prowinds),col_offset,max(pcolinds))});
        legend('Maximum Location')
        set(gca,'YDir','Normal')
      end
      %Save the results
      convout{rowidx+1,colidx+1}=insarconv;  convMARout(rowidx+1,colidx+1)=max(max(insarconv))/mean(mean(insarconv));
      rows = linspace(0,size(madata,1),size(insarconv,1));
      cols = linspace(0,size(madata,2),size(insarconv,2));
      rowout(rowidx+1,colidx+1)=size(madata,1)-rows(round(row_offset)); 
      colout(rowidx+1,colidx+1)=size(madata,2)-cols(round(col_offset));
      %Update start row
      rowidx = rowidx+1;
    end
    %Update start col
    colidx=colidx+1;
  end
  %% Process the data
  %Remove results with low MAR (keeping results above some threshold)
  thresh = 35;
  logicvec = convMARout<=thresh;
  colout(logicvec) = nan; rowout(logicvec) = nan;
  Out(fn_id).row = nanmean(nanmean(rowout)); Out(fn_id).col = nanmean(nanmean(colout));
  
  keyboard
end



%% Use cross correlation to determine offset in x and y indices

%% Plot resulting interferogram
