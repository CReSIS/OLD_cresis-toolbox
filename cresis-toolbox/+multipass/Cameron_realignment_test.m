%% Correlation Window approach
maxlag = 250; % minimum bin offset used by cameron
corr_window = 4* max_lag; 
% width of correlation window, real corr_window is +1

% make nan arrays the length of each profile for surface and bed
% surface
AT_data.window_S_11 = NaN(1,length(pass(1).layers(1).layer_elev));
AT_data.window_S_14 = NaN(1,length(pass(2).layers(1).layer_elev));
AT_data.window_S_18 = NaN(1,length(pass(3).layers(1).layer_elev));
% Bed
AT_data.window_B_11 = NaN(1,length(pass(1).layers(2).layer_elev));
AT_data.window_B_14 = NaN(1,length(pass(2).layers(2).layer_elev));
AT_data.window_B_18 = NaN(1,length(pass(3).layers(2).layer_elev));

% temporary correlation window file
temp_corr = NaN(length(AT_data.elevS.P2011),2*max_lag + 1);
xcorr = [];
temp_shift = [];

%%
[temp_corr,lags] = xcorr(tmp_signal_1, tmp_signal_2, maxlag,'coeff');

%% 2011-2014 surface adjustment
for cidx = ((corr_window/2)+1):(length(AT_data.window_S_11) - (corr_window/2) - 1)
  % create tmp_signals 1 and 2
  tmp_signal_1 = pass(1).layers(1).layer_elev(cidx - (corr_window/2):cidx+(corr_window/2));
  tmp_signal_2 = AT_data.elevS.P2014(cidx - (corr_window/2):cidx+(corr_window/2)); 
  
  
  [temp_corr(cidx,:),temp_lag] = xcorr(tmp_signal_1 - nanmean(tmp_signal_1),tmp_signal_2 - nanmean(tmp_signal_2), maxlag, 'coeff');
  vel_pick_window = (-10:10) + max_lag + 1 + AT_data.shift_data_14;
 
  
  if sum(temp_corr(cidx, max_lag + 1:end) == max(temp_corr(cidx, maxlag + 1:end))) ==1
    temp_shift = temp_lag(temp_corr(cidx, vel_pick_window) == max(temp_corr(cidx, vel_pick_window))) + vel_pick_window(1) -1;
    if temp_shift > vel_pick_window(1) - max_lag - 1
      AT_data.window_S_14(cidx) = temp_shift;
    end
  end
end
%% 
%2011-2014 bed adjustment
for cidx = ((corr_window/2)+1):(length(AT_data.window_B_11) - (corr_window/2) - 1)
  % create tmp_signals 1 and 2
  tmp_signal_1 = AT_data.elevB.P2011(cidx - (corr_window/2):cidx+(corr_window/2));
  tmp_signal_2 = AT_data.elevB.P2014(cidx - (corr_window/2):cidx+(corr_window/2)); 
  
  [temp_corr(cidx,:),temp_lag] = xcorr(tmp_signal_1 - mean(tmp_signal_1),(tmp_signal_2 - mean(tmp_signal_2)), maxlag, 'coeff');
  vel_pick_window = (-10:10) + max_lag + 1 + AT_data.shift_data_14;
  

  if sum(temp_corr(cidx, max_lag + 1:end) == max(temp_corr(cidx, maxlag + 1:end))) ==1
    temp_shift = temp_lag(temp_corr(cidx, vel_pick_window) == max(temp_corr(cidx, vel_pick_window))) + vel_pick_window(1) -1;
    if temp_shift > vel_pick_window(1) - max_lag - 1
      AT_data.window_B_14(cidx) = temp_shift;
    end
  end
end


