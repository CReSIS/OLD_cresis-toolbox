% This function is used to detect and cancel burst noise 
% (i.e.those thin vertical stripes in echogram )---qishi
function [data_out]=icards_burst_noise_detection(data_in)
time_start=tic;

test_data=data_in;
filter_length=51;%FIR and median filter length
noise_threshold=7;%threshold in dB
consecutive_ones=5;%determine the shortest length of noise could be regarded as burst noise
num_col_fix=4;%number of columns to be averaged to get a value to replace the noise points(to the oneside)
%====================find burst noise======================================
test_data_power=real(test_data).^2+imag(test_data).^2; %get signal power     
test_data_dec=fir_dec(test_data_power,ones(1,filter_length)/filter_length);%FIR filter
test_data_med=medfilt1(double(test_data_dec),filter_length,[],2);%median filter
mask=(lp(test_data_power)>lp(test_data_med)+noise_threshold);%possible burst noise data points


mask_extend=[zeros(consecutive_ones-1,size(mask,2));...
    mask;zeros(consecutive_ones-1,size(mask,2))];%extend original mask up and down with zeros for shift accumulation
mask_add=zeros(size(mask));%initial the accumulated data
for ii=1:2*consecutive_ones-1% accumulation of shifted mask
  mask_add=mask_add+mask_extend(ii:ii+size(mask,1)-1,:);
end
[idx_row,idx_col,~]=find(mask_add>=consecutive_ones);%find the index of those consequtive ones
mask_final=zeros(size(mask));% initial the final mask
idx_row_u=unique(idx_row);%those row index where consequtive ones exist
for ii=1:length(idx_row_u)
  [idx_row_col,~]=find(idx_row==idx_row_u(ii));%find the corresponding column index of each row index
  mask_final(idx_row_u(ii),idx_col(idx_row_col)')=1;
end
%====================fix burst noise=======================================
test_data_after=test_data;
[~,idx_col,~]=find(mask_final==1);
idx_col=unique(idx_col);
% for ii=1:length(idx_col)
%   idx_row=find(mask_final(:,idx_col(ii))==1);
%   if idx_col(ii)-num_col_fix<=0 %noise points lie in the column which is too near to the first column
%     test_data_after(idx_row,idx_col(ii))=mean(test_data(idx_row,idx_col(ii)+1:idx_col(ii)+num_col_fix*2),2);
%   elseif idx_col(ii)+num_col_fix>size(test_data,2)%noise points lie in the column which is too near to the last column
%     test_data_after(idx_row,idx_col(ii))=mean(test_data(idx_row,idx_col(ii)-num_col_fix*2:idx_col(ii)-1),2);
%   else 
%     test_data_after(idx_row,idx_col(ii))=mean(test_data(idx_row,[idx_col(ii)-num_col_fix:idx_col(ii)-1,[],idx_col(ii)+1:idx_col(ii)+num_col_fix]),2);
%   end
% end
test_data_after=test_data.*(~mask_final);
data_out=test_data_after;
time_spent=toc(time_start);
fprintf('%f s is used to detect burst noise\n',time_spent);
end


