function Data = detrending(data)

xx=double(data);
xx1=10*log10(xx(:,:));
[M N]=size(xx1);
xy=zeros(M,N);
%xy_before_denoise=zeros(M,N);
xinput=zeros(2^nextpow2(M),1);
win=151;
for ii=1:N;  
  [detrended_data, x_trend] = detrending_method(xx1(:,ii), win, 2);       
  %xy_before_denoise(:,ii)=detrended_data;  
  
  xinput(1:M,1)=detrended_data;
  
  [sigDEN,wDEC] = func_denoise_sw1d(xinput);
   xy(:,ii)=sigDEN(1:M);
 end

Data=xy;
return